#!/opt/local/bin/python
# -*- encoding: ascii -*-
"""
    @purpose: Solving the
    @author: Brittany Hall
    @date: 07.10.2017
    @version: 0.1
    @updates:
"""
#from casadi import *
import scipy
from col_cstr_model import *
from numpy import floor, zeros

def system(t, x, u, T):
    uc = u
    
    def _col_cstr_LV(t, X, uc):
        """
            Inputs- uc[0] - reflux LT
            uc[1] - boilup VB
            These inputs are set by altering col_LV.py file.
            Outputs-x - liquid coposition and hold up for stages 1 to NT
        """
        NT = params['dist']['NT']          #Number of stages in distillation column
        LT = uc[0]
        print LT
        raw_input()
        VB = uc[1]
        F = uc[2]
        D = uc[3]
        B = uc[4]
        F_0 = params['dist']['F_0']                                  #CSTR feedrate
        zF = params['dist']['zF'][0]                          #Feed composition [A]
        
        #Collecting inputs and disturbances
        u_all = append(LT,VB)
        u_all = append(u_all, D)
        u_all = append(u_all, B)
        u_all = append(u_all, F)
        u_all = append(u_vall, zF)
        u_all = append(u_all, F_0)
        
        xprime = col_cstr_model(t,X,u_all)
        return xprime
    #Specifying the integrator to use
    ode15s =  scipy.integrate.ode(_col_cstr_LV)
    ode15s.set_integrator('vode', method='bdf', order =15)
    
    #Set time range
    t_start = t
    t_final = t+T
    delta_t = 1
    num_steps = int(floor((t_final-t_start)/delta_t) + 1)
    
    #Set initial conditions
    x_init = x
    ode15s.set_initial_value(x_init, t_start)
    
    #Create vectors to store trajectories
    t = zeros((num_steps,1))
    x_out = zeros((num_steps,len(x)))
    t[0] = t_start
    x_out[0,:] = x_init
    
    #Integrate the ODE across each delta_t timestep
    k = 1
    while ode15s.successful() and k<num_steps:
        ode15s.integrate(ode15s.t + delta_t)
        #Store values for later
        t[k] = ode15s.t
        x_out[k] = ode15s.y[0]
        k += 1
    
    raw_input()
    lengthx = shape(x_out)
    y = transpose(x_out[lengthx[0],:])

    return y
