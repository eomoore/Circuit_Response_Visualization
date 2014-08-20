__author__ = 'moore'

from flask_wtf import Form
from wtforms import RadioField, FloatField, SelectField
from wtforms.validators import DataRequired


class Baseform(Form):
    vType = SelectField('Source Selector', choices=[(1, 'Vs * u(t - t0)'), (2, 'Vs * u(-t - t0)'),
                                                    (3, 'A * u(-t) + B * u(t)'), (4, 'Vs * (u(t) - u(t-t0))')])
    vSource = FloatField('Vs(V)', default=0, validators=[DataRequired()])
    a_value = FloatField('A(V)', default=0, id='a', validators=[DataRequired()])
    b_value = FloatField('B(V)', default=0, id='b', validators=[DataRequired()])
    t_0 = FloatField('t0', default=0)
    resistance = FloatField('R', default=1, validators=[DataRequired()])
    rOrder = RadioField('', choices=[(1, '&Omega;'), (1000, 'k&Omega;'), (1000000, 'M&Omega;')], default=1, coerce=int)


class RCForm(Baseform):
    capacitance = FloatField('C', default=1, validators=[DataRequired()])
    cOrder = RadioField('', choices=[(1, 'F'), (0.001, 'mF'), (0.000001, 'uF'), (0.000000001, 'nF')], default=1,
                        coerce=float)


class LRForm(Baseform):
    inductance = FloatField('L', default=1, validators=[DataRequired()])
    lOrder = RadioField('', choices=[(1, 'H'), (0.001, 'mH')], default=1, coerce=float)


class LRCForm(RCForm):
    inductance = FloatField('L', default=1, validators=[DataRequired()])
    lOrder = RadioField('', choices=[(1, 'H'), (1000, 'KH'), (1000000, 'MH')], default=1)