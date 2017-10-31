#!/opt/local/bin/python
# -*- encoding: ascii -*-
"""
    @purpose: Simulation of column and CSTR with LV-configuration.
    Function calls the model col_model and includes control of the condenser
    and reboiler using two P-controllers with the LV-configuration.
    @author: Brittany Hall
    @date: 05.10.2017
    @version: 0.1
    @updates:
"""
from numpy import *
from col_cstr_model import *
from params import *
from system import system

def col_cstr_LV(t, X):
    """
        Inputs- uc[0] - reflux LT
                uc[1] - boilup VB
                These inputs are set by altering col_LV.py file.
        Outputs-x - liquid coposition and hold up for stages 1 to NT
    """
    NT = params['dist']['NT']          #Number of stages in distillation column
    print system()
    raw_input()
    LT = system.uc[0]
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
