from flask import render_template, request, send_file
import numpy as np
import cStringIO
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
        # return render_template('rcplot.html', form=form, title="RC Circuit")


def plot(value1=None, value2 = None):
    """

    :param value1:
    :param value2:
    :return:
    """
    global form


    if value1 == 'rc':
        form = RCForm()

        if request.method == 'POST':
            (v_s, v_t, time_0, r, r_order, c, c_order) = get_values(form)
            v, i, t, tau, i_tau, v_tau = solver(v_source=v_s, v_type=v_t, t_0=time_0, res=r * r_order, cap=c * c_order,
                                                circuit_type=0)

            if value2 == 'annotate':
                return plot(v, i, t, time_0, tau, i_tau, v_tau, 'rc', v_t, annotate=1)

            else:
                return plot(v, i, t, time_0, tau, i_tau, v_tau, 'rc', v_t, annotate=0)

    elif value1 == 'lr':
        form = LRForm()
        if request.method == 'POST':
            (v_s, v_t, time_0, r, r_order, ind, ind_order) = get_values(form)

            v, i, t, tau, i_tau, v_tau = solver(v_source=v_s, v_type=v_t, t_0=time_0, res=r * r_order, ind=ind * ind_order,
                                                circuit_type=1)

            if value2 == 'annotate':
                return plot(v, i, t, time_0, tau, i_tau, v_tau, 'lr', v_t, annotate=1)

            else:
                return plot(v, i, t, time_0, tau, i_tau, v_tau, 'lr', v_t, annotate=0)

    else:
        return render_template('404.html'), 404


def rcplot(value=0):
    global form
    form = RCForm()
    # if request.method == 'POST' and form.validate():
    # if form.validate_on_submit():
    if request.method == 'POST':
        (v_s, v_t, time_0, r, r_order, c, c_order) = get_values(form)
        v, i, t, tau, i_tau, v_tau = solver(v_source=v_s, v_type=v_t, t_0=time_0, res=r * r_order, cap=c * c_order,
                                            circuit_type=0)

        if value == 'annotate':
            return plot(v, i, t, time_0, tau, i_tau, v_tau, 'rc', v_t, annotate=1)

        else:
            return plot(v, i, t, time_0, tau, i_tau, v_tau, 'rc', v_t, annotate=0)


def lr():
    """

    

    :rtype : object
    :return: 
    """
    global form
    form = LRForm()

    if request.method == 'GET':
        return render_template('lrplot.html', form=form, title="LR Circuit")


def lrplot(value=0):
    global form
    form = LRForm(request.form)
    if request.method == 'POST':
        (v_s, v_t, time_0, r, r_order, ind, ind_order) = get_values(form)

        v, i, t, tau, i_tau, v_tau = solver(v_source=v_s, v_type=v_t, t_0=time_0, res=r * r_order, ind=ind * ind_order,
                                            circuit_type=1)

        if value == 'annotate':
            return plot(v, i, t, time_0, tau, i_tau, v_tau, 'lr', v_t, annotate=1)

        else:
            return plot(v, i, t, time_0, tau, i_tau, v_tau, 'lr', v_t, annotate=0)


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
        return render_template('plot.html', form=form, title="LRC Circuit")

    elif request.method == 'POST':
        return render_template('plot.html', form=form, title="LRC Circuit")


def lrc_parallel():
    global form

    form = LRCForm()

    if request.method == 'GET':
        return render_template('plot.html', form=form, title="LRC Circuit")

    elif request.method == 'POST':
        return render_template('plot.html', form=form, title="LRC Circuit")


def lrcplot():
    global form

    # if request.method == 'POST' and rcform.validate():
    if request.method == 'GET' or request.method == 'POST':
        # if rcform.validate_on_submit():
        # if request.method == 'POST':
        Vs = form.vSource.data
        r = form.resistance.data
        rOrderOfMag = form.rOrder.data
        l = form.inductance.data
        lOrderOfMag = form.lOrder.data
        tau = r * rOrderOfMag * l * lOrderOfMag
        # tau = r * l

        # Generate the plot
        x = np.linspace(0, 7 * tau)
        Vc = Vs * (1 - (np.exp(-x / tau)))
        i = (Vs - Vc) / r

        plt.clf()
        plt.figure(1)
        plt.subplot(211)
        plt.plot(x, Vc)
        plt.ylabel('Voltage (v)')

        plt.subplot(212)
        plt.plot(x, i)

        plt.ylabel('Current (A)')
        plt.xlabel('Time (s)')

        f = cStringIO.StringIO()
        plt.savefig(f, format='png')

        # Serve up the data
        f.seek(0)

    return send_file(f, mimetype='image/png')


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


def heaviside1(x, v):

    x = np.array(x)
    if x.shape != ():

        y = np.zeros(x.shape)
        for ind, temp in np.ndenumerate(x):
            if x[ind] >= v:
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


def solver(circuit_type=0, v_source=0, v_type=1, t_0=0, res=0, cap=0, ind=0):
    """


    :param circuit_type:
        The type of circuit
        0 = RC circuit
        1 = LR circuit
    :param v_source:
        The type of voltage source
        '1' = V * u(t)
        '2' = V * u(-t)
        '3' = V * u(t)-u(t-t0)
    :param v_type:
    :param t_0:
    :param res:
    :param cap:
    :param ind:
    :rtype : v, i, t, tau, i_tau, v_tau
            v is the voltage vector
            i is the current vector
            t is the time vector
            tau is the value of tau
            i_tau is the value of the current at t = tau
            v_tau is the value of the voltage at t = tau
    """
    v, i, t, tau, v_tau, i_tau = 0, 0, 0, 0, 0, 0
    if circuit_type == 0:                                       # RC Circuit
        tau = res * cap
        t = np.linspace(-tau + t_0, (7 * tau) + t_0, 100)
        if v_type == '1':                                       # V * u(t)
            v = v_source * (1 - (np.exp(-(t - t_0) / tau)))
            v *= heaviside(v)
            i = (v_source - v) / res
            v_tau = v_source * (1 - (np.exp(-tau / tau)))
            i_tau = (v_source - v_tau) / res

        elif v_type == '2':                                     # V * u(-t)
            v = v_source * (np.exp(-(t - t_0) / tau))
            v_temp = v_source * (np.exp(-0 / tau))
            v = heaviside1(v, v_temp)
            i = (0 - v) / res
            v_tau = v_source * (np.exp(-tau / tau))
            i_tau = (0 - v_tau) / res

        elif v_type == '3':                                     # V * u(t) - u(t-t0)
            t = np.linspace(-2*tau , (7 * tau) + t_0, 100)
            v1 =  v_source * (1 - (np.exp(-t / tau)))
            v2 =  v_source * (1 - (np.exp(-(t-t_0) / tau)))
            v1 = heaviside(v1) * v1
            v2 = heaviside(v2) * v2
            v = v1 - v2
            i = (v_source - v) / res
            v_tau = v_source * (1 - (np.exp(-tau / tau)))
            i_tau = (v_source - v_tau) / res


    elif circuit_type == 1:                                     # LR Circuit
        tau = ind / res
        t = np.linspace(-tau + t_0, (7 * tau) + t_0)

        if v_type == '1':                                       # V * u(t)
            i = (v_source / res) * (1 - (np.exp(-(t-t_0) / tau)))
            i *= heaviside(i)
            v = v_source - (i * res)
            i_tau = (v_source / res) * (1 - (np.exp(-tau / tau)))
            v_tau = v_source - (i_tau * res)

        elif v_type == '2':                                     # V * u(-t)
            i = (v_source / res) * (np.exp(-(t-t_0) / tau))
            i_temp = (v_source / res) * (np.exp(-0 / tau))
            i = heaviside1(i, i_temp)
            v = 0 - (i * res)
            i_tau = (v_source / res) * (np.exp(-tau / tau))
            v_tau = 0 - (i_tau * res)

        elif v_type == '3':
            t = np.linspace(-2*tau , (7 * tau) + t_0, 100)
            i1 = (v_source / res) * (1 - (np.exp(-t / tau)))
            i1 *= heaviside(i1)
            i2 = (v_source / res) * (1 - (np.exp(-(t-t_0) / tau)))
            i2 *= heaviside(i2)
            i = i1 - i2
            v = v_source - (i * res)
            i_tau = (v_source / res) * (1 - (np.exp(-tau / tau)))
            v_tau = v_source - (i_tau * res)

    elif circuit_type == 2:
        pass
    elif circuit_type == 3:
        pass

    return v, i, t, tau, i_tau, v_tau


def get_values(form):
    try:
        Vs = form.vSource.data
        Vt = form.vType.data
        t0 = form.t_0.data
        r = form.resistance.data
        rOrderOfMag = form.rOrder.data
        c = form.capacitance.data
        cOrderOfMag = form.cOrder.data

        return Vs, Vt, t0, r, rOrderOfMag, c, cOrderOfMag

    except AttributeError:
        Vs = form.vSource.data
        Vt = form.vType.data
        t0 = form.t_0.data
        r = form.resistance.data
        rOrderOfMag = form.rOrder.data
        l = form.inductance.data
        lOrderOfMag = form.lOrder.data

        return Vs, Vt, t0, r, rOrderOfMag, l, lOrderOfMag


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
    a.plot(t, v)

    # Add a label for the y-axis
    plt.ylabel('Voltage (v)')

    if annotate and ((v_type == '1' and circ_type == 'rc') or (v_type == '2' and circ_type == 'lr')):
        plt.plot([0 + t_0, tau+t_0], [v.min(), v.max()], color='red', linewidth=1, linestyle='--')
        plt.annotate('tau = %s' % tau, xy=(tau, v_tau), xycoords='data',
                     xytext=(+tau, +v_tau), textcoords='offset points', fontsize=16,
                     arrowprops=dict(arrowstyle="->", connectionstyle="arc3,rad=.2"))

    elif annotate and ((v_type == '2' and circ_type == 'rc') or (v_type == '1' and circ_type == 'lr')):
        plt.plot([0+t_0, tau+t_0], [v.max(), v.min()], color='red', linewidth=1, linestyle='--')
        plt.annotate('tau = %s' % tau, xy=(tau, v_tau), xycoords='data',
                     xytext=(+tau, +v_tau), textcoords='offset points', fontsize=16,
                     arrowprops=dict(arrowstyle="->", connectionstyle="arc3,rad=.2"))

    elif annotate and ((v_type == '3' and circ_type == 'rc') or (v_type == '3' and circ_type == 'lr')):
        plt.plot([0,0], [v.min(), v.max()], color='red', linewidth=1, linestyle='--')
        plt.plot([0,t_0], [v.max(), v.max()], color='red', linewidth=1, linestyle='--')
        plt.plot([t_0, t_0], [v.min(), v.max()], color='red', linewidth=1, linestyle='--')


        # Create a new subplot from a 2x1 grid
    plt.subplot(212)

    # Plot the Current
    plt.plot(t, i)

    # Add a label for the y-axis
    plt.ylabel('Current (A)')
    # Add a label for the y-axis
    plt.xlabel('Time (s)')

    if annotate and ((v_type == '1' and circ_type == 'rc') or (v_type == '2' and circ_type == 'lr')):
        plt.plot([0+t_0, tau+t_0], [i.max(), i.min()], color='red', linewidth=1, linestyle='--')
        plt.annotate('tau = %s' % tau, xy=(tau, i_tau), xycoords='data',
                     xytext=(+0, +0), textcoords='offset points', fontsize=16,
                     arrowprops=dict(arrowstyle="->", connectionstyle="arc3,rad=.2"))

    elif annotate and ((v_type == '2' and circ_type == 'rc') or (v_type == '1' and circ_type == 'lr')):
        plt.plot([0+t_0, tau+t_0], [i.min(), i.max()], color='red', linewidth=1, linestyle='--')
        plt.annotate('tau = %s' % tau, xy=(tau, v_tau), xycoords='data',
                     xytext=(+0, +0), textcoords='offset points', fontsize=16,
                     arrowprops=dict(arrowstyle="->", connectionstyle="arc3,rad=.2"))

    elif annotate and ((v_type == '3' and circ_type == 'rc') or (v_type == '3' and circ_type == 'lr')):
        plt.plot([0,0], [i.min(), i.max()], color='red', linewidth=1, linestyle='--')
        plt.plot([0,t_0], [i.max(), i.max()], color='red', linewidth=1, linestyle='--')
        plt.plot([t_0, t_0], [i.min(), i.max()], color='red', linewidth=1, linestyle='--')

    return mpld3.fig_to_html(fig)
