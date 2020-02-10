# -*- coding: utf-8 -*-
"""
Created on Tue Apr 30 15:01:36 2019

@author: Chaofan
"""
from wtforms import Form, FloatField, validators, SelectField

working_fluid_type = [('Isobutane', 'Isobutane'), ('Isopentane', 'Isopentane'), ('R245fa', 'R245fa')]

class InputForm(Form):
    fluid_type = SelectField(
            u'Type of working fluid', choices = working_fluid_type,
            validators=[validators.InputRequired()])
    Q_w = FloatField(
        label='Mass flow rate of the hot water (kg/s)', default=30.0,
        validators=[validators.InputRequired()])
    Q_s = FloatField(
        label='Mass flow rate of the steam (kg/s)', default=1.0,
        validators=[validators.InputRequired()])
    T_b_p = FloatField(
        label='Brine temperature of the product well (K)', default=383.15,
        validators=[validators.InputRequired()])
    T_b_i = FloatField(
        label='Brine temperature of the injection well (K)', default=328.15,
        validators=[validators.InputRequired()])
    p_b = FloatField(
        label='Pressure of the brine (Pa)', default=2e5,
        validators=[validators.InputRequired()])
    x_c = FloatField(
        label='Dryness fraction (steam mass fraction)', default=0,
        validators=[validators.InputRequired()])
    T_env = FloatField(
        label='Environment temperature (K)', default=293.15,
        validators=[validators.InputRequired()])
    eta = FloatField(
        label='Heat transfer efficiency between brine and working fluid', default=1,
        validators=[validators.InputRequired()])

