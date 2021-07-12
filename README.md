# ORCSimulator

## Installation and Usage

After downloading install the requirements within a fresh virtual environment:

```sh
python -m pip install -r requirements.txt
```

Three different tasks can be accomplished with the software:

- single parameter influence analysis (`from geothermal_orc import single_parameter_influence`),
- single parameter optimization (`from geothermal_orc import single_optimization`) and
- multivariate parameter optimization (`from geothermal_orc import multivariate_optimization`)

Import the respective functionality and run it with the .json control file.
The result date are passed to the function call. Example applications can be
found within the scripts folder:

- single.py
- parameter_influence.py
- multivariate.py

Sample input data can be found in the same folder.

## Topology and System Design

![Alt flow diagram of the geothermal ORC](./flowdiagram.svg)
*Flow diagram of the geothermal ORC plant (a) and Ts-diagram of the process (b)*

The table below highlights the design parameters of the system.

| Item                | Parameter            | Symbol                  | Value | Unit  |
|:--------------------|:---------------------|:------------------------|:------|:------|
| Geothermal resource | Steam temperature    | *T*<sub>gs</sub>    | 140   |  °C   |
|                     | Steam mass flow rate | *ṁ*<sub>gs</sub>    | 180   |  kg/s |
|                     | Brine temperature    | *T*<sub>gb</sub>    | 140   |  °C   |
|                     | Brine mass flow rate | *ṁ*<sub>gb</sub>    | 20    |  kg/s |
|                     | Steam mass fraction  | *x*                     | 0.1   | \-    |
|                     | Brine/steam pressure | *p*<sub>geo</sub> | 3.615 | bar   |
| Ambient condition   | Average temperature  | *T*<sub>am</sub>    | 5     |  °C   |
|                     | Average pressure     | *p*<sub>am</sub>    | 0.6   |  bar  |

| Location             | Parameter                             | Symbol                 | Value | Unit |
|:---------------------|:--------------------------------------|:-----------------------|------:|:----:|
| Turbine              | Isentropic efficiency                 | *η*<sub>s, t</sub>     |    90 |  %   |
| Feed pump            | Isentropic efficiency                 | *η*<sub>s, fp</sub>    |    75 |  %   |
| Air fan              | Isentropic efficiency                 | *η*<sub>s, af</sub>    |    60 |  %   |
| Main condenser       | Upper terminal temperature difference | *ΔT*<sub>t, u</sub>  |    10 |  °C  |
|                      | Pressure ratio on hot side            | *pr*<sub>1</sub>     |     1 |  \-  |
|                      | Pressure ratio on cold side           | *pr*<sub>2</sub>     | 0.995 |  \-  |
| Geo-steam evaporator | Pressure ratio on hot side            | *pr*<sub>1</sub>     |     1 |  \-  |
|                      | Pressure ratio on cold side           | *pr*<sub>2</sub>     |     1 |  \-  |
| Geo-brine evaporator | Pinch point temperature difference    | *ΔT*<sub>pp</sub>    |     8 |  °C  |
|                      | Pressure ratio on hot side            | *pr*<sub>1</sub>     |  0.98 |  \-  |
|                      | Pressure ratio on cold side           | *pr*<sub>2</sub>     |     1 |  \-  |
| IHE & preheater      | Pressure ratio on hot side            | *pr*<sub>1</sub>     |  0.98 |  \-  |
|                      | Pressure ratio on cold side           | *pr*<sub>2</sub>     |  0.98 |  \-  |
| Preheater outlet     | Approach point temperature difference | *ΔT*<sub>ap</sub>    |     2 |  °C  |
| Generator            | Efficiency                            | *η*<sub>el, mech</sub> |    97 |  %   |
| Motors               | Efficiency                            | *η*<sub>el, mech</sub> |    97 |  %   |

## Configuring the input file

It is possible to choose the decision variables for the optimization as well as
the objective of the optimization. Additionally, other boundary conditions can
be specified. Apart from the geosteam share three variables or boundary
conditions must be specified in total. Available parameters are listed below:

| Parameter     | Meaning                                                            |
|---------------|--------------------------------------------------------------------|
| p_before_tur  | pressure at connection 1                                           |
| T_before_tur  | temperature at connection 1                                        |
| T_reinjection | temperature at connection 35                                       |
| brine_evap_Td | temperature change from connection 33 to 34 (negative)             |
| Q_brine_ev    | heat transferred by brine evaporator                               |
| dT_air        | temperature change from connection 21 to 22                        |
| IHE_sizing    | IHE sizing factor: 0 = IHE non existent; 1 = maximum heat transfer |
| Q_ihe         | heat transferred by IHE                                            |

### Boundary conditions

- It is **mandatory** to specify the geosteam share.
- For single parameter investigation specify two additional parameters
- For multivariate investigation specify less than two additional parameters

### Variables

- Specify variable(s) to investigate with upper and lower limit as well as
  tolerance in case of single optimization

### Results

- Specify components, connections and misc for retrieving the respective data
  in the results DataFrames.

Parameters for misc:

| Parameter          | Meaning                                                              |
|--------------------|----------------------------------------------------------------------|
| gross power output | power output of the rankine cycle only                               |
| thermal efficiency | efficiency considering gross power output                            |
| net power output   | power output including the power required for the condensator's fans |
| net efficiency     | efficiency considering net power output                              |
| IHE sizing factor  | IHE sizing factor (see above)                                        |

### Possible objectives

Choose from:

- "gross power output"
- "net power output"

## Reference

An archived version of this repository can be found at zenodo:
https://zenodo.org/record/SOMERECORD. For more information also see the
respective publication. A link will be added here once published.

## License

Copyright (c) 2021 Francesco Witte, Chaofan Chen

Permission is hereby granted, free of charge, to any person obtaining a copy of
this software and associated documentation files (the "Software"), to deal in
the Software without restriction, including without limitation the rights to
use, copy, modify, merge, publish, distribute, sublicense, and/or sell copies of
the Software, and to permit persons to whom the Software is furnished to do so,
subject to the following conditions:

The above copyright notice and this permission notice shall be included in all
copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS
FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR
COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER
IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN
CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.
