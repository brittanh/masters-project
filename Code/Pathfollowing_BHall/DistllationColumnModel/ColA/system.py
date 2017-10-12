#!/opt/local/bin/python
# -*- encoding: ascii -*-
"""
    @purpose: Solving the
    @author: Brittany Hall
    @date: 07.10.2017
    @version: 0.1
    @updates:
"""
from casadi import *
import scipy
from col_cstr_LV import *

def system(t, x, u, T):
    global uc
    uc = u
    ode15s =  scipy.integrate.ode(col_cstr_LV)
    ode15s.set_integrator('vode', method='bdf', order=15, nsteps=3000)
    ode15s.set_initial_value(x, t)
    x_out = ode15s.integrate(t+T)
    lengthx = shape(x_out)
    y = transpose(x_out[lengthx[0],:])
    return y
