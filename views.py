from flask import render_template, request
import numpy as np
from scipy.integrate import odeint, simps

import matplotlib

matplotlib.use('Agg')
import matplotlib.pyplot as plt

plt.ioff()
import mpld3
# from mpld3 import plugins


# Import all the forms
from Circuitforms import *

form = 0


def index():
    return render_template('index.html')


def rc():
    global form
    form = RCForm()

    if request.method == 'GET':
        return render_template('plot.html', form=form, type='rc', title="RC Circuit")


def lr():
    global form
    form = LRForm()

    if request.method == 'GET':
        return render_template('plot.html', form=form, type='lr', title="LR Circuit")


def lrc(value):
    if value is 'series':
        lrc_series()
    elif value is 'parallel':
        lrc_parallel()
    else:
        return render_template('404.html'), 404


def lrc_series():
    global form

    form = LRCForm()

    if request.method == 'GET':
        return render_template('lrcplot.html', form=form, type='lrcseries', title="Series LRC Circuit")

    elif request.method == 'POST':
        return render_template('lrcplot.html', form=form, type='lrcseries', title="Series LRC Circuit")


def lrc_parallel():
    global form

    form = LRCForm()

    if request.method == 'GET':
        return render_template('lrcplot.html', form=form, type='lrcparallel', title="Parallel LRC Circuit")

    elif request.method == 'POST':
        return render_template('lrcplot.html', form=form, type='lrcparallel', title="Parallel LRC Circuit")


def plot1(value1=None, value2=None):
    global form

    if value1 == 'rc':
        form = RCForm()
        if request.method == 'POST':
            v_s, v_t, time_0, r, r_order, c, c_order, v_a, v_b = get_values(form)
            v, i, t, tau, i_tau, v_tau = solver(v_source=v_s, v_type=v_t, t_0=time_0, res=r * r_order, cap=c * c_order,
                                                circuit_type=0, v_a=v_a, v_b=v_b)

            if value2 == 'annotate':
                return plot(v, i, t, time_0, tau, i_tau, v_tau, 'rc', v_t, annotate=1)

            else:
                return plot(v, i, t, time_0, tau, i_tau, v_tau, 'rc', v_t, annotate=0)

    elif value1 == 'lr':
        form = LRForm()
        if request.method == 'POST':
            (v_s, v_t, time_0, r, r_order, ind, ind_order, v_a, v_b) = get_values(form)

            v, i, t, tau, i_tau, v_tau = solver(v_source=v_s, v_type=v_t, t_0=time_0, res=r * r_order,
                                                ind=ind * ind_order, circuit_type=1, v_a=v_a, v_b=v_b)

            if value2 == 'annotate':
                return plot(v, i, t, time_0, tau, i_tau, v_tau, 'lr', v_t, annotate=1)

            else:
                return plot(v, i, t, time_0, tau, i_tau, v_tau, 'lr', v_t, annotate=0)

    elif value1 == 'lrcseries':
        form = LRCForm()
        if request.method == 'POST':
            (v_s, v_t, r, r_order, c, c_order, ind, ind_order) = get_values1(form)

            v_c, v_r, v_l, i, t = solver(v_source=v_s, v_type=v_t, res=r * r_order, cap=c*c_order,
                                         ind=ind * ind_order, circuit_type=2)

            if value2 == 'annotate':
                return lrcplot(v_c, v_r, v_l, i, t)

            else:
                return lrcplot(v_c, v_r, v_l, i, t)

    elif value1 == 'lrcparallel':
        form = LRCForm()
        if request.method == 'POST':
            (v_s, v_t, r, r_order, c, c_order, ind, ind_order) = get_values1(form)

            v, i_c, i_r, i_l, t = solver(v_source=v_s, v_type=v_t, res=r * r_order, cap=c*c_order,
                                         ind=ind * ind_order, circuit_type=3)

            if value2 == 'annotate':
                return lrcplot1(v, i_c, i_r, i_l, t)

            else:
                return lrcplot1(v, i_c, i_r, i_l, t)

    else:
        return render_template('404.html'), 404


def lrcplot(v_c, v_r, v_l, i, t):
    plt.clf()

    # Create a new figure
    fig = plt.figure(1)

    # Create a new subplot from a 2x1 grid
    a = plt.subplot(211)

    # Plot the Voltage
    a.plot(t, v_c, linewidth=2.5)
    a.plot(t, v_r, linewidth=2.5)
    a.plot(t, v_l, linewidth=2.5)

    # Add a label for the y-axis
    plt.ylabel('Voltage (v)')

    # Create a new subplot from a 2x1 grid
    plt.subplot(212)

    # Plot the Current
    plt.plot(t, i, linewidth=2.5)

    # Add a label for the y-axis
    plt.ylabel('Current (A)')
    # Add a label for the y-axis
    plt.xlabel('Time (s)')

    return mpld3.fig_to_html(fig)


def lrcplot1(v, i_c, i_r, i_l, t):
    plt.clf()

    # Create a new figure
    fig = plt.figure(1)

    # Create a new subplot from a 2x1 grid
    a = plt.subplot(211)

    # Plot the Voltage
    a.plot(t, v, linewidth=2.5)

    # Add a label for the y-axis
    plt.ylabel('Voltage (v)')

    # Create a new subplot from a 2x1 grid
    plt.subplot(212)

    # Plot the Current
    plt.plot(t, i_c, linewidth=2.5)
    plt.plot(t, i_r, linewidth=2.5)
    plt.plot(t, i_l, linewidth=2.5)

    # Add a label for the y-axis
    plt.ylabel('Current (A)')
    # Add a label for the y-axis
    plt.xlabel('Time (s)')

    return mpld3.fig_to_html(fig)


def heaviside(x):
    x = np.array(x)
    if x.shape != ():
        y = np.zeros(x.shape)
        y[x < 0.0] = 0
        y[x > 0.0] = 1
        y[x == 0.0] = 0.5
    else:  # special case for 0d array (a number)
        if x > 0:
            y = 1
        elif x == 0:
            y = 0.5
        else:
            y = 0
    return y


def heaviside1(x, value):
    x = np.array(x)
    if x.shape != ():
        y = np.zeros(x.shape)
        for ind, temp in np.ndenumerate(x):
            if x[ind] >= value:
                y[ind] = value
            else:
                y[ind] = temp

    else:  # special case for 0d array (a number)
        if x < 0:
            y = 1
        elif x == 0:
            y = 0.5
        else:
            y = 0
    return y


def heaviside2(x, v):
    x = np.array(x)
    if x.shape != ():
        y = np.zeros(x.shape)
        for ind, temp in np.ndenumerate(x):
            if x[ind] <= v:
                y[ind] = v
            else:
                y[ind] = temp

    else:  # special case for 0d array (a number)
        if x < 0:
            y = 1
        elif x == 0:
            y = 0.5
        else:
            y = 0
    return y


def solver(circuit_type=0, v_source=0, v_type=1, t_0=0, res=0, cap=0, ind=0, v_a=0, v_b=0):
    """



    :type v_source: object
    :param circuit_type:
        The type of circuit
        0 = RC circuit
        1 = LR circuit
        2 = Series LRC Circuit
        3 = Parallel LRC Circuit
    :param v_source:
        The magnitude of the voltage
    :param v_type:
        The type of voltage source
        '1' = V * u(t)
        '2' = V * u(-t)
        '3' = V * u(t)-u(t-t0)
    :param t_0:
        Time delay
    :param res:
        Resistance magnitude
    :param cap:
        Capacitance magnitude
    :param ind:
        Inductance magnitude
    :rtype : v, i, t, tau, i_tau, v_tau
            v is the voltage vector
            i is the current vector
            t is the time vector
            tau is the value of tau
            i_tau is the value of the current at t = tau
            v_tau is the value of the voltage at t = tau
    """
    tau, v_tau, i_tau = 0, 0, 0

    if circuit_type == 0:                                                   # RC Circuit
        tau = res * cap
        t = np.linspace(-tau + t_0, (7 * tau) + t_0, 500)
        v = np.zeros(t.shape)
        i = np.zeros(t.shape)
        if v_type == '1':                                                   # V * u(t)
            v = v_source * (1 - (np.exp(-(t - t_0) / tau)))
            v *= heaviside(v)

            i = np.gradient(v) * cap
            v_tau = v_source * (1 - (np.exp(-tau / tau)))
            i_tau = (v_source - v_tau) / res

        elif v_type == '2':                                                 # V * u(-t)
            v = v_source * (np.exp(-(t - t_0) / tau))
            v_temp = v_source * (np.exp(-0 / tau))
            v = heaviside1(v, v_temp)
            i = np.gradient(v) * cap
            v_tau = v_source * (np.exp(-tau / tau))
            i_tau = (0 - v_tau) / res

        elif v_type == '3':                                                 # A*u(-t) + B*u(t)
            v = v_b + ((v_a - v_b)*(np.exp(-(t - t_0) / tau)))

            if v_b > v_a:
                for x in np.nditer(v, op_flags=['readwrite']):
                    if x < v_a:
                        x[...] = v_a

            elif v_b < v_a:
                for x in np.nditer(v, op_flags=['readwrite']):
                    if x > v_a:
                        x[...] = v_a

            i = np.gradient(v) * cap
            # v_tau = v_b + ((v_a - v_b)*(np.exp(-(tau - t_0) / tau)))
            # i_tau = (v_source - v_tau)/res

        elif v_type == '4':                                                 # V * u(t) - u(t-t0)
            t = np.linspace(-2 * tau, (7 * tau) + t_0, 500)
            v1 = v_source * (1 - (np.exp(-t / tau)))
            v2 = v_source * (1 - (np.exp(-(t - t_0) / tau)))
            v1 *= heaviside(v1)
            v2 *= heaviside(v2)
            v = v1 - v2
            i = np.gradient(v) * cap
            v_tau = v_source * (1 - (np.exp(-tau / tau)))
            i_tau = (v_source - v_tau) / res

    elif circuit_type == 1:                                                 # LR Circuit
        tau = ind / res
        t = np.linspace(-tau + t_0, (7 * tau) + t_0, 500)
        v = np.zeros(t.shape)
        i = np.zeros(t.shape)
        if v_type == '1':  # V * u(t)
            i = (v_source / res) * (1 - (np.exp(-(t - t_0) / tau)))
            i *= heaviside(i)
            v = np.gradient(i) * ind * -1
            i_tau = (v_source / res) * (1 - (np.exp(-tau / tau)))
            v_tau = v_source - (i_tau * res)

        elif v_type == '2':                                                 # V * u(-t)
            i = (v_source / res) * (np.exp(-(t - t_0) / tau))
            i_temp = (v_source / res) * (np.exp(-0 / tau))
            i = heaviside1(i, i_temp)
            v = np.gradient(i) * ind * -1
            i_tau = (v_source / res) * (np.exp(-tau / tau))
            v_tau = 0 - (i_tau * res)

        elif v_type == '3':                                                 # A*u(-t) + B*u(t)
            i = (v_b / res) + ((v_a - v_b) / res) * (np.exp(-(t - t_0) / tau))

            if v_b > v_a:
                for x in np.nditer(i, op_flags=['readwrite']):
                    if x < (v_a/res):
                        x[...] = v_a/res

            elif v_b < v_a:
                for x in np.nditer(v, op_flags=['readwrite']):
                    if x > (v_a/res):
                        x[...] = v_a/res

            v = np.gradient(i) * ind * -1
            # i_tau = (v_b / res) + ((v_a - v_b) / res) * (np.exp(-(tau) / tau))
            # v_tau = v_source - (i_tau * res)

        elif v_type == '4':                                                 # V * u(t) - u(t-t0)
            t = np.linspace(-2 * tau, (7 * tau) + t_0, 100)
            i1 = (v_source / res) * (1 - (np.exp(-t / tau)))
            i1 *= heaviside(i1)
            i2 = (v_source / res) * (1 - (np.exp(-(t - t_0) / tau)))
            i2 *= heaviside(i2)
            i = i1 - i2
            v = np.gradient(i) * ind * -1
            i_tau = (v_source / res) * (1 - (np.exp(-tau / tau)))
            v_tau = v_source - (i_tau * res)

    elif circuit_type == 2:                                                 # Series LRC Circuit
        v0 = v_source
        i0 = 0
        i0_prime = (-1/ind) * (res*i0 + v0)

        def deriv(i, t):
            return np.array([i[1], (v_source/ind) - (res/ind)*i[0] - i[1]/(ind*cap)])

        t = np.linspace(0, 10, 1001)
        i_init = np.array([i0, i0_prime], dtype=np.float)

        q, i = odeint(deriv, i_init, t).T

        v_c = q/cap
        v_r = res * i
        v_l = np.gradient(i) * ind

        return v_c, v_r, v_l, i, t

    elif circuit_type == 3:                                                 # Parallel LRC Circuit
        a = -1/(res * cap)
        b = -1/(ind * cap)
        v0 = 0
        i0 = v_source
        v0_prime = (-1/cap) * (i0 + (v0/res))

        def deriv(v, t):
            return np.array([v[1], (v_source/cap)-v[0]/(ind * cap) + -v[1]/(res * cap)])

        t = np.linspace(0.0, 10.0, 1001)
        v_init = np.array([v0, v0_prime], dtype=np.float)

        x, v = odeint(deriv, v_init, t).T

        i_l = x/ind
        i_r = v/res
        i_c = np.gradient(v) * cap

        return v, i_c, i_r, i_l, t

    return v, i, t, tau, i_tau, v_tau


def get_values(form):
    try:
        vs = form.vSource.data
        vt = form.vType.data
        t0 = form.t_0.data
        r = form.resistance.data
        rorderofmag = form.rOrder.data
        c = form.capacitance.data
        corderofmag = form.cOrder.data
        v_a = form.a_value.data
        v_b = form.b_value.data

        return vs, vt, t0, r, rorderofmag, c, corderofmag, v_a, v_b

    except AttributeError:
        vs = form.vSource.data
        vt = form.vType.data
        t0 = form.t_0.data
        r = form.resistance.data
        rorderofmag = form.rOrder.data
        l = form.inductance.data
        lorderofmag = form.lOrder.data
        v_a = form.a_value.data
        v_b = form.b_value.data

        return vs, vt, t0, r, rorderofmag, l, lorderofmag, v_a, v_b


def get_values1(form):
    try:
        vs = form.vSource.data
        vt = form.vType.data
        r = form.resistance.data
        rorderofmag = form.rOrder.data
        c = form.capacitance.data
        corderofmag = form.cOrder.data
        l = form.inductance.data
        lorderofmag = form.lOrder.data

        return vs, vt, r, rorderofmag, c, corderofmag, l, lorderofmag

    except AttributeError:
        pass


def plot(v, i, t, t_0, tau, i_tau, v_tau, circ_type, v_type, annotate=0):
    """



    :param v:
    :param i:
    :param t:
    :param tau:
    :param i_tau:
    :param v_tau:
    :param circ_type:
    :param annotate:
    :param v_type:
    :rtype : object
    """
    plt.clf()

    # Create a new figure
    fig = plt.figure(1)

    # Create a new subplot from a 2x1 grid
    a = plt.subplot(211)

    # Plot the Voltage
    a.plot(t, v, linewidth=2.5)

    # Add a label for the y-axis
    plt.ylabel('Voltage (v)')

    if annotate and ((v_type == '1' and circ_type == 'rc') or (v_type == '2' and circ_type == 'lr')):
        plt.plot([0 + t_0, tau + t_0], [v.min(), v.max()], color='red', linewidth=1, linestyle='--')
        plt.annotate('tau = %s' % tau, xy=(tau, v_tau), xycoords='data',
                     xytext=(+tau, +v_tau), textcoords='offset points', fontsize=16)

    elif annotate and ((v_type == '2' and circ_type == 'rc') or (v_type == '1' and circ_type == 'lr')):
        plt.plot([0 + t_0, tau + t_0], [v.max(), v.min()], color='red', linewidth=1, linestyle='--')
        plt.annotate('tau = %s' % tau, xy=(tau, v_tau), xycoords='data',
                     xytext=(+tau, +v_tau), textcoords='offset points', fontsize=16)

    elif annotate and ((v_type == '4' and circ_type == 'rc') or (v_type == '4' and circ_type == 'lr')):
        plt.plot([0, 0], [v.min(), v.max()], color='red', linewidth=1, linestyle='--')
        plt.plot([0, t_0], [v.max(), v.max()], color='red', linewidth=1, linestyle='--')
        plt.plot([t_0, t_0], [v.min(), v.max()], color='red', linewidth=1, linestyle='--')

    # Create a new subplot from a 2x1 grid
    plt.subplot(212)

    # Plot the Current
    plt.plot(t, i, linewidth=2.5)

    # Add a label for the y-axis
    plt.ylabel('Current (A)')
    # Add a label for the y-axis
    plt.xlabel('Time (s)')

    if annotate and ((v_type == '1' and circ_type == 'rc') or (v_type == '2' and circ_type == 'lr')):
        plt.plot([0 + t_0, tau + t_0], [i.max(), i.min()], color='red', linewidth=1, linestyle='--')
        plt.annotate('tau = %ss' % tau, xy=(tau, i_tau), xycoords='data',
                     xytext=(+0, +0), textcoords='offset points', fontsize=16)

    elif annotate and ((v_type == '2' and circ_type == 'rc') or (v_type == '1' and circ_type == 'lr')):
        plt.plot([0 + t_0, tau + t_0], [i.min(), i.max()], color='red', linewidth=1, linestyle='--')
        plt.annotate('tau = %ss ' % tau, xy=(tau, v_tau), xycoords='data',
                     xytext=(+0, +0), textcoords='offset points', fontsize=16)

    # elif annotate and ((v_type == '4' and circ_type == 'rc') or (v_type == '4' and circ_type == 'lr')):
    #     plt.plot([0, 0], [i.min(), i.max()], color='red', linewidth=1, linestyle='--')
    #     plt.plot([0, t_0], [i.max(), i.max()], color='red', linewidth=1, linestyle='--')
    #     plt.plot([t_0, t_0], [i.min(), i.max()], color='red', linewidth=1, linestyle='--')

    return mpld3.fig_to_html(fig)