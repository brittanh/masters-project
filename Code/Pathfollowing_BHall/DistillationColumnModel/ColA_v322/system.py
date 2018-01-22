#!/opt/local/bin/python
# -*- encoding: ascii -*-
"""
    @purpose: Solving the
    @author: Brittany Hall
    @date: 07.10.2017
    @version: 0.1
    @updates:
"""
import scipy
from col_cstr_model import col_cstr_model
from numpy import floor, zeros, append, shape, transpose
from params import *

def system(t, x, u, T):
    global uc
    uc = u
    def _col_cstr_LV(t, X):
        global uc
        """
            Inputs- uc[0] - reflux LT
                    uc[1] - boilup VB
                    uc[2] - feed F
                    uc[3] - distillate D
                    uc[4] - bottoms B
                    
            Outputs-x - liquid coposition and hold up for stages 1 to NT
        """
        NT = params['dist']['NT']          #Number of stages in column
        LT = uc[0]
        VB = uc[1]
        F = uc[2]
        D = uc[3]
        B = uc[4]
        F_0 = params['dist']['F_0']                       #CSTR feedrate
        zF = params['dist']['zF'][0]               #Feed composition [A]
        
        #Collecting inputs and disturbances
        u_all = append(LT,VB)
        u_all = append(u_all, D)
        u_all = append(u_all, B)
        u_all = append(u_all, F)
        u_all = append(u_all, zF)
        u_all = append(u_all, F_0)

        xprime = col_cstr_model(t,X,u_all)
        return xprime
    
    #Specifying the integrator to use
    ode15s =  scipy.integrate.ode(_col_cstr_LV)
    ode15s.set_integrator('vode', method='bdf', with_jacobian=False)
    
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
    while ode15s.successful() and k < num_steps:
        ode15s.integrate(ode15s.t + delta_t)
        
        #Store values for later
        t[k] = ode15s.t + delta_t
        x_out[k] = ode15s.y
        k += 1
    
    lengthx = shape(x_out)
    y = transpose(x_out[lengthx[0]-1,:])

    return y
