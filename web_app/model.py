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
        label='Mass flow rate of the hot water (kg/s)', default=190.56,
        validators=[validators.InputRequired()])
    Q_s = FloatField(
        label='Mass flow rate of the steam (kg/s)', default=20.28,
        validators=[validators.InputRequired()])
    T_b_p = FloatField(
        label='Brine temperature of the production well (C)', default=146.6,
        validators=[validators.InputRequired()])
    T_b_i = FloatField(
        label='Brine temperature of the injection well (C)', default=69.1,
        validators=[validators.InputRequired()])
    T_env = FloatField(
        label='Environment temperature (C)', default=-4.7,
        validators=[validators.InputRequired()])
    p_env = FloatField(
        label='Environment pressure (bar)', default=0.61,
        validators=[validators.InputRequired()])
