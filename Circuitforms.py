__author__ = 'moore'

from flask_wtf import Form
from wtforms import RadioField, FloatField, SelectField


class Baseform(Form):
    vSource = FloatField('Vs(V)', default=0)
    vType = SelectField('', choices=[(1, 'u(t - t0)') , (2, 'u(-t - t0)'), (3, '(u(t) - u(t-t0)')])
    t_0 = FloatField('t0', default=0)
    resistance = FloatField('R', default=1)
    rOrder = RadioField('', choices=[(1, '&Omega;'), (1000, 'k&Omega;'), (1000000, 'M&Omega;')], default=1, coerce=int)


class RCForm(Baseform):
    capacitance = FloatField('C', default=1)
    cOrder = RadioField('', choices=[(1, 'F'), (0.001, 'mF'), (0.000001, 'uF'), (0.000000001, 'nF')], default=1,
                        coerce=float)


class LRForm(Baseform):
    inductance = FloatField('L', default=1)
    lOrder = RadioField('', choices=[(1, 'H'), (0.001, 'mH')], default=1, coerce=float)


class LRCForm(RCForm):
    inductance = FloatField('L', default=1)
    lOrder = RadioField('', choices=[(1, 'H'), (1000, 'KH'), (1000000, 'MH')], default=1)