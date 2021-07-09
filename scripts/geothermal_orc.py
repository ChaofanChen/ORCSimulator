from tespy.connections import Connection, Bus, Ref
from tespy.networks import Network
from tespy.components import (
    HeatExchanger, Pump, Turbine, Source, Sink, CycleCloser, Splitter,
    Merge, Condenser, Drum, Valve, HeatExchangerSimple, Compressor
)
from tespy.tools.helpers import UserDefinedEquation, TESPyNetworkError
from tespy.tools.fluid_properties import h_mix_pT, T_mix_ph

import CoolProp as CP
from CoolProp.CoolProp import PropsSI as PSI

from collections import OrderedDict
import numpy as np
import pygmo as pg
import pandas as pd

import os


def desuperheat(ude):
    ttd_min = ude.params['ttd_min']
    return (
        ude.conns[0].h.val_SI - ude.conns[1].h.val_SI -
        ude.params['distance'] * (
            ude.conns[0].h.val_SI - h_mix_pT(
                ude.conns[1].get_flow(),
                T_mix_ph(ude.conns[2].get_flow()) + ttd_min)))


def desuperheat_deriv(ude):
    ude.jacobian[ude.conns[0]][2] = 1 - ude.params['distance']
    ude.jacobian[ude.conns[1]][1] = ude.numeric_deriv('p', 1)
    ude.jacobian[ude.conns[1]][2] = -1
    ude.jacobian[ude.conns[2]][1] = ude.numeric_deriv('p', 2)
    ude.jacobian[ude.conns[2]][2] = ude.numeric_deriv('h', 2)
    return ude.jacobian


class PowerPlant():

    def __init__(self, working_fluid):
        """Set up model."""
        self.working_fluid = working_fluid
        fluids = ['water', self.working_fluid, 'air']
        self.nw = Network(fluids=fluids)
        self.nw.set_attr(p_unit='bar', T_unit='C', h_unit='kJ / kg')

        # geo parameters

        self.geo_mass_flow = 200
        geo_steam_share = 0.1
        self.T_brine_in = 140

        # ambient parameters

        self.T_amb = 5
        self.p_amb = 0.6

        # main components

        geo_steam = Source('geosteam source')
        geo_brine = Source('geobrine source')
        geo_reinjection = Sink('re-injection')

        air_in = Source('air source')
        air_out = Sink('air sink')
        air_fan = Compressor('air fan')
        air_cond = Condenser('condenser')

        orc_cc = CycleCloser('orc cycle closer')

        evap_splitter = Splitter('splitter evaporation')
        evap_merge = Merge('merge evaporation')
        evap_steam = Condenser('geosteam evaporator')
        evap_brine = HeatExchanger('geobrine evaporator')
        dr = Drum('drum')
        geo_merge = Merge('merge brine')

        pre = HeatExchanger('preheater')
        feed_working_fluid_pump = Pump('feed pump')

        tur = Turbine('turbine')

        ihe = HeatExchanger('internal heat exchanger')

        # busses
        net_power = Bus('net power output')
        net_power.add_comps(
            {'comp': tur, 'char': 0.97},
            {'comp': feed_working_fluid_pump, 'char': 0.97, 'base': 'bus'},
            {'comp': air_fan, 'char': 0.97, 'base': 'bus'}
        )

        ORC_power_bus = Bus('cycle gross power output')
        ORC_power_bus.add_comps(
            {'comp': tur}, {'comp': feed_working_fluid_pump}
        )

        geothermal_bus = Bus('thermal input')
        geothermal_bus.add_comps(
            {'comp': pre, 'char': -1}, {'comp': evap_brine, 'char': -1},
            {'comp': evap_steam, 'char': -1}
        )

        self.nw.add_busses(net_power, ORC_power_bus, geothermal_bus)

        # turbine to condenser
        c1 = Connection(orc_cc, 'out1', tur, 'in1', label='1')
        c2 = Connection(tur, 'out1', ihe, 'in1', label='2')
        c3 = Connection(ihe, 'out1', air_cond, 'in1', label='3')
        self.nw.add_conns(c1, c2, c3)

        # condenser to steam generator
        c4 = Connection(air_cond, 'out1', feed_working_fluid_pump, 'in1', label='4')
        c5 = Connection(feed_working_fluid_pump, 'out1', ihe, 'in2', label='5')
        self.nw.add_conns(c4, c5)

        # steam generator
        c6 = Connection(ihe, 'out2', pre, 'in2', label='6')
        c7 = Connection(pre, 'out2', dr, 'in1', label='7')
        c8 = Connection(dr, 'out1', evap_splitter, 'in1', label='8')
        c9 = Connection(evap_splitter, 'out2', evap_steam, 'in2', label='9')
        c10 = Connection(evap_steam, 'out2', evap_merge, 'in2', label='10')
        c11 = Connection(evap_splitter, 'out1', evap_brine, 'in2', label='11')
        c12 = Connection(evap_brine, 'out2', evap_merge, 'in1', label='12')
        c13 = Connection(evap_merge, 'out1', dr, 'in2', label='13')
        c0 = Connection(dr, 'out2', orc_cc, 'in1', label='0')
        self.nw.add_conns(c6, c7, c8, c11, c9, c12, c10, c13, c0)

        # condenser cold side
        c20 = Connection(air_in, 'out1', air_fan, 'in1', label='20')
        c21 = Connection(air_fan, 'out1', air_cond, 'in2', label='21')
        c22 = Connection(air_cond, 'out2', air_out, 'in1', label='22')
        self.nw.add_conns(c20, c21, c22)

        # geo source
        c30 = Connection(geo_steam, 'out1', evap_steam, 'in1', label='30')
        c31 = Connection(evap_steam, 'out1',  geo_merge, 'in1', label='31')
        c32 = Connection(geo_brine, 'out1', geo_merge, 'in2', label='32')
        c33 = Connection(geo_merge, 'out1', evap_brine, 'in1', label='33')
        self.nw.add_conns(c30, c31, c32, c33)

        c34 = Connection(evap_brine, 'out1', pre, 'in1', label='34')
        c35 = Connection(pre, 'out1', geo_reinjection, 'in1', label='35')
        self.nw.add_conns(c34, c35)

        # generate a set of stable starting values of every working fluid
        # fluid settings
        c6.set_attr(fluid={self.working_fluid: 1.0, 'air': 0.0, 'water': 0.0})
        c20.set_attr(fluid={self.working_fluid: 0.0, 'air': 1.0, 'water': 0.0})
        c30.set_attr(fluid={self.working_fluid: 0.0, 'air': 0.0, 'water': 1.0})
        c32.set_attr(fluid={self.working_fluid: 0.0, 'air': 0.0, 'water': 1.0})

        # connection parameters
        p0 = PSI('P', 'T', self.T_brine_in + 273.15, 'Q', 1, self.working_fluid)
        c1.set_attr(p0=p0 / 1e5)
        ws_stable_h0 = (
            PSI('H', 'T', self.T_amb + 273.15, 'Q', 1, self.working_fluid) +
            0.5 * (
                PSI('H', 'T', self.T_brine_in + 273.15, 'Q', 1, self.working_fluid) -
                PSI('H', 'T', self.T_amb + 273.15, 'Q', 1, self.working_fluid)
            )
        ) / 1e3
        c2.set_attr(h=ws_stable_h0)
        p0 = PSI('P', 'T', self.T_amb + 273.15, 'Q', 1, self.working_fluid)
        c3.set_attr(Td_bp=5, design=['Td_bp'], p0=p0 / 1e5)
        c5.set_attr(h=Ref(c4, 1, 1))

        # steam generator
        c30.set_attr(
            m=self.geo_mass_flow * geo_steam_share,
            T=self.T_brine_in, x=1, p0=5)
        c32.set_attr(
            m=self.geo_mass_flow * (1 - geo_steam_share),
            T=self.T_brine_in, x=0)

        c13.set_attr()
        c12.set_attr(x=0.5)
        c10.set_attr(x=0.5, design=['x'])
        c34.set_attr(h=Ref(c33, 1, -50))

        c7.set_attr(Td_bp=-2)

        # main condenser
        c20.set_attr(p=self.p_amb, T=self.T_amb)
        c22.set_attr(T=self.T_amb + 15, p=self.p_amb)

        # component parameters
        # condensing
        ihe.set_attr(pr1=0.98, pr2=0.98)
        air_cond.set_attr(pr1=1, pr2=0.995, ttd_u=10)
        air_fan.set_attr(eta_s=0.6)

        # steam generator
        evap_brine.set_attr(pr1=0.98, ttd_l=8)
        pre.set_attr(pr1=0.98, pr2=0.98)

        self.nw.set_attr(iterinfo=False)
        self.nw.solve('design')
        self.nw.save('stable_' + self.working_fluid)

        # specify actual parameters
        tur.set_attr(eta_s=0.9)
        feed_working_fluid_pump.set_attr(eta_s=0.75)
        c2.set_attr(h=None)
        c5.set_attr(h=None)
        c34.set_attr(h=None, T=Ref(c33, 1, -10))

        self.nw.solve('design')
        c22.set_attr(T=None)
        c3.set_attr(Td_bp=None)

        self.ude_IHE_size = UserDefinedEquation(
            label='ihe deshuperheat ratio',
            func=desuperheat, deriv=desuperheat_deriv,
            latex={
                'equation':
                    r'0 = h_3 - h_2 - x_\mathrm{IHE} \cdot \left(h_3 -'
                    r'h\left(p_2, T_5 + \Delta T_\mathrm{t,u,min} \right)'
                    r'\right)'},
            conns=[
                self.nw.get_conn('2'),
                self.nw.get_conn('3'),
                self.nw.get_conn('5')],
            params={'distance': 0.0, 'ttd_min': 2}
        )
        if self.nw.lin_dep or self.nw.res[-1] > 1e-3:
            msg = 'No stable solution found.'
            raise TESPyNetworkError(msg)
        print(
            'Generated stable starting values for working fluid ' +
            self.working_fluid + '.')

    def run_simulation(
            self, p_before_tur=None, Q_ihe=None, Q_brine_ev=None,
            T_before_tur=None, T_reinjection=None, brine_evap_Td=None,
            dT_air=None, IHE_sizing=None, geo_steam_share=None):
        """Run simulation on specified parameter set."""

        self.nw.get_comp('internal heat exchanger').set_attr(Q=Q_ihe)
        self.nw.get_conn('1').set_attr(p=p_before_tur, T=T_before_tur)
        self.nw.get_conn('35').set_attr(T=T_reinjection)
        self.nw.get_comp('geobrine evaporator').set_attr(Q=Q_brine_ev)

        if geo_steam_share is not None:
            self.nw.get_conn('30').set_attr(
                m=self.geo_mass_flow * geo_steam_share)
            self.nw.get_conn('32').set_attr(
                m=self.geo_mass_flow * (1 - geo_steam_share))

        if brine_evap_Td is not None:
            self.nw.get_conn('34').set_attr(
                T=Ref(self.nw.get_conn('33'), 1, brine_evap_Td))
        else:
            self.nw.get_conn('34').set_attr(T=None)

        if dT_air is not None:
            self.nw.get_conn('22').set_attr(T=Ref(self.nw.get_conn('21'), 1, dT_air))
        else:
            self.nw.get_conn('22').set_attr(T=None)

        if IHE_sizing is None:
            if self.ude_IHE_size in self.nw.user_defined_eq.values():
                self.nw.del_ude(self.ude_IHE_size)
            self.nw.get_comp('internal heat exchanger').set_attr(pr1=0.98, pr2=0.98)
        else:
            if self.ude_IHE_size not in self.nw.user_defined_eq.values():
                self.nw.add_ude(self.ude_IHE_size)
            self.ude_IHE_size.params['distance'] = IHE_sizing
            if IHE_sizing == 0:
                self.nw.get_comp('internal heat exchanger').set_attr(pr1=1, pr2=1)
            else:
                self.nw.get_comp('internal heat exchanger').set_attr(pr1=0.98, pr2=0.98)

        try:
            self.nw.solve('design')
#            self.nw.print_results()
        except ValueError:
            self.nw.res = [1]
            pass

    def check_simulation(self, value):
        """Check if simulation converged."""
        if self.nw.lin_dep or self.nw.res[-1] > 1e-3:
            self.nw.solve(
                'design', init_path='stable_' + self.working_fluid,
                init_only=True)
            return np.nan
        else:
            for cp in self.nw.comps['object']:
                if isinstance(cp, HeatExchanger):
                    if cp.Q.val > 0:
                        print(cp.label)
                        return np.nan
                    elif cp.kA.val <= 0 or (np.isnan(cp.kA.val) and cp.Q.val != 0):
                        print(cp.label)
                        return np.nan
        return value

    def get_power(self):
        """Calculate ORC gross power (main cycle only)."""
        return self.check_simulation(self.nw.busses['cycle gross power output'].P.val)

    def get_net_power(self):
        """Calculate net power."""
        return self.check_simulation(self.nw.busses['net power output'].P.val)

    def get_thermal_efficiency(self):
        """Calculate thermal efficiency."""
        return self.check_simulation(
            -self.nw.busses['cycle gross power output'].P.val /
            self.nw.busses['thermal input'].P.val)

    def get_net_efficiency(self):
        """Calculate net efficiency."""
        return self.check_simulation(
            -self.nw.busses['net power output'].P.val /
            self.nw.busses['thermal input'].P.val)

    def get_geosteam_share(self):
        """Return a geosteam share."""
        return self.check_simulation(
            self.nw.get_conn('geosteam').m.val_SI / self.geo_mass_flow)

    def get_connection_param(self, conn, param):
        """Return a connection parameter."""
        return self.check_simulation(
            self.nw.get_conn(conn).get_attr(param).val)

    def get_component_param(self, comp, param):
        """Return a component parameter."""
        return self.check_simulation(
            self.nw.get_comp(comp).get_attr(param).val)

    def get_misc_param(self, param):
        """Get non component or connection parameters."""
        if param == 'gross power output':
            return self.get_power()
        elif param == 'net power output':
            return self.get_net_power()
        elif param == 'thermal efficiency':
            return self.get_thermal_efficiency()
        elif param == 'net efficiency':
            return self.get_net_efficiency()
        elif param == 'IHE sizing factor':
            return self.ude_IHE_size.params['distance']

    def get_objective_func(self, objective):
        """Return corresponding objective function."""
        if objective == 'net power output':
            return self.get_net_power
        elif objective == 'gross power output':
            return self.get_power
        else:
            msg = (
                'Please specify valid objective function: "net power output" '
                'or "gross power output".')
            raise ValueError(msg)


class _MultivariateOptimizationProblem():

    def fitness(self, x):
        """Fitness function."""
        inputs = {}
        i = 0
        for val in x:
            inputs.update({list(self.params_to_opt.keys())[i]: x[i]})
            i += 1

        self.model.run_simulation(**inputs, **self.boundary_conditions)
        f1 = self.objective()
        ci1 = 70 - self.model.get_connection_param('35', 'T')
        return [f1, ci1]

    def get_nobj(self):
        """Return number of objectives."""
        return 1

    # equality constraints
    def get_nec(self):
        return 0

    # inequality constraints
    def get_nic(self):
        return 1

    def get_bounds(self):
        """Return bounds of decision variables."""
        return (
            [self.params_to_opt[param]['min']
             for param in self.params_to_opt.keys()],
            [self.params_to_opt[param]['max']
             for param in self.params_to_opt.keys()])


def _process_generation_data(
        pop, gen, individuals, params_list, objectives_list,
        constraint_list):

    individual = 0
    for x in pop.get_x():
        for i in range(len(x)):
            individuals.loc[[(gen, individual)], params_list[i]] = x[i]
        # individuals.loc[(gen, individual), params_list] = x
        individual += 1

    individual = 0
    for objective in pop.get_f():
        individuals.loc[
            [(gen, individual)], objectives_list + constraint_list] = objective
        individual += 1

    individuals['valid'] = individuals[constraint_list] < 0

    return individuals


def multivariate_optimization(**input_data):

    fluid_list = input_data['working_fluid_list']

    result = {}

    opt_results = pd.DataFrame()

    for fluid in fluid_list:
        print('+' * 75)
        print('Working fluid:', fluid)
        print('+' * 75)

        optimize = _MultivariateOptimizationProblem()
        optimize.model = PowerPlant(working_fluid=fluid)

        optimize.params_to_opt = OrderedDict(input_data['variables'])
        optimize.boundary_conditions = input_data['boundary_conditions']

        optimize.objective = optimize.model.get_objective_func(
            input_data['objective'])
        objectives_list = [input_data['objective']]
        constraint_list = ['constraints']
        params_list = list(optimize.params_to_opt.keys())

        # run optimization problem
        prob = pg.problem(optimize)
        num_gen = input_data['num_gen']
        num_ind = input_data['num_ind']
        pop = pg.population(prob, size=num_ind)
        algo = pg.algorithm(pg.ihs(gen=num_gen))

        individuals = pd.DataFrame(
            columns=params_list + objectives_list + constraint_list,
            index=[(gen, ind) for gen in range(num_gen)
                   for ind in range(num_ind)])

        gen = 0
        for gen in range(num_gen - 1):
            individuals = _process_generation_data(
                pop, gen, individuals, params_list, objectives_list,
                constraint_list)
            print()
            print('Evolution: {}'.format(gen))
            for i in range(len(objectives_list)):
                print(objectives_list[i] +
                      ': {}'.format(round(-pop.champion_f[i]/1e6, 4)))
            for i in range(len(params_list)):
                print(params_list[i] +
                      ': {}'.format(round(pop.champion_x[i], 4)))
            pop = algo.evolve(pop)

        gen += 1
        individuals = _process_generation_data(
            pop, gen, individuals, params_list, objectives_list,
            constraint_list)

        print()

        optimize.fitness(pop.champion_x)

        for conn, param_list in input_data['result_data']['connections'].items():
            for param in param_list:
                opt_results.loc[fluid, param + '_' + conn] = optimize.model.get_connection_param(conn, param)
        for comp, param_list in input_data['result_data']['components'].items():
            for param in param_list:
                opt_results.loc[fluid, param + '_' + comp] = optimize.model.get_component_param(comp, param)

        for param in input_data['result_data']['misc']:
            opt_results.loc[fluid, param] = optimize.model.get_misc_param(param)

        print('Final evolution: {}'.format(gen))
        for i in range(len(objectives_list)):
            print(objectives_list[i] +
                  ': {}'.format(round(pop.champion_f[i], 4)))
        for i in range(len(params_list)):
            print(params_list[i] +
                  ': {}'.format(round(pop.champion_x[i], 4)))

        print()

        result[fluid] = individuals

    if input_data['save_result']:
        if not os.path.isdir(input_data['scenario']):
            os.mkdir('./' + input_data['scenario'])
        for fluid in fluid_list:
            result[fluid].to_csv(input_data['scenario'] + '/' + fluid + '.csv')

        opt_results.to_csv(input_data['scenario'] + '/champion_data.csv')

    return result, opt_results


def _golden_ratio_search(
        function, get_param, a, b, tol=1e-5, direction='min',
        param_to_opt=None, func_params={}):
    """Golden ratio search"""

    invphi = (5 ** 0.5 - 1) / 2  # 1 / phi
    invphi2 = (3 - 5 ** 0.5) / 2  # 1 / phi^2

    (a, b) = (min(a, b), max(a, b))
    h = b - a
    if h <= tol:
        return np.array([a, b])

    # Required steps to achieve tolerance
    n = int(np.ceil(np.log(tol / h) / np.log(invphi)))

    c = a + invphi2 * h
    d = a + invphi * h

    if direction == 'min':
        func_params[param_to_opt] = c
        function(**func_params)
        yc = get_param()
        func_params[param_to_opt] = d
        function(**func_params)
        yd = get_param()
    else:
        func_params[param_to_opt] = c
        function(**func_params)
        yc = 1 / get_param()
        func_params[param_to_opt] = d
        function(**func_params)
        yd = 1 / get_param()

    for k in range(n-1):
        if yc < yd:
            b = d
            d = c
            yd = yc
            h = invphi * h
            c = a + invphi2 * h
            func_params[param_to_opt] = c
            if direction == 'min':
                function(**func_params)
                yc = get_param()
            else:
                function(**func_params)
                yc = 1 / get_param()
        else:
            a = c
            c = d
            yc = yd
            h = invphi * h
            d = a + invphi * h
            func_params[param_to_opt] = d
            if direction == 'min':
                function(**func_params)
                yd = get_param()
            else:
                function(**func_params)
                yd = 1 / get_param()

    if yc < yd:
        return np.array([a, d])
    else:
        return np.array([c, b])


def single_optimization(**input_data):

    fluid_list = input_data['working_fluid_list']
    result = {}

    for key, value in input_data['variables'].items():
        result[key] = pd.DataFrame(index=fluid_list)
        for fluid in fluid_list:
            print('+' * 75)
            ORC = PowerPlant(working_fluid=fluid)
            print('Working fluid:', fluid)

            p_opt = _golden_ratio_search(
                ORC.run_simulation,
                ORC.get_objective_func(input_data['objective']),
                a=value['min'], b=value['max'], tol=value['tol'],
                param_to_opt=key,
                func_params=input_data['boundary_conditions'])

            for conn, param_list in input_data['result_data']['connections'].items():
                for param in param_list:
                    result[key].loc[fluid, param + '_' + conn] = ORC.get_connection_param(conn, param)
            for comp, param_list in input_data['result_data']['components'].items():
                for param in param_list:
                    result[key].loc[fluid, param + '_' + comp] = ORC.get_component_param(comp, param)

            for param in input_data['result_data']['misc']:
                result[key].loc[fluid, param] = ORC.get_misc_param(param)

    if input_data['save_result']:
        if not os.path.isdir(input_data['scenario']):
            os.mkdir('./' + input_data['scenario'])
        for key in input_data['variables'].keys():
            result[key].to_csv(input_data['scenario'] + '/' + key + '.csv')

    return result


def single_parameter_influence(**input_data):
    """Evaluate the parameter influence of a single parameter."""

    fluid_list = input_data['working_fluid_list']

    result = {}

    for fluid in fluid_list:
        ORC = PowerPlant(working_fluid=fluid)
        print('Working fluid:', fluid)

        result[fluid] = pd.DataFrame()

        for key, value in input_data['variables'].items():
            for x in np.linspace(value['min'], value['max'], value['num']):
                ORC.run_simulation(**input_data['boundary_conditions'], **{key: x})

                for conn, param_list in input_data['result_data']['connections'].items():
                    for param in param_list:
                        result[fluid].loc[x, param + '_' + conn] = ORC.get_connection_param(conn, param)
                for comp, param_list in input_data['result_data']['components'].items():
                    for param in param_list:
                        result[fluid].loc[x, param + '_' + comp] = ORC.get_component_param(comp, param)

                for param in input_data['result_data']['misc']:
                    result[fluid].loc[x, param] = ORC.get_misc_param(param)

    if input_data['save_result']:
        if not os.path.isdir(input_data['scenario']):
            os.mkdir('./' + input_data['scenario'])
        for fluid in fluid_list:
            result[fluid].to_csv(input_data['scenario'] + '/' + fluid + '.csv')

    return result
