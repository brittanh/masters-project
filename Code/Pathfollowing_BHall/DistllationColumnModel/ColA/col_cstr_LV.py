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

def col_cstr_LV(t, X):
    """
        Inputs- uc[0] - reflux LT
                uc[1] - boilup VB
                These inputs are set by altering col_LV.py file.
        Outputs-x - liquid coposition and hold up for stages 1 to NT
    """
    NT = 41                            #Number of stages in distillation column

    LT = uc[0]
    VB = uc[1]
    F = uc[2]
    D = uc[3]
    B = uc[4]
    F_0 = 0.3                                                    #CSTR feedrate
    zF = 1.0                                              #Feed composition [A]

    #Collecting inputs and disturbances
    u_all = append(LT,VB)
    u_all = append(u_all, D)
    u_all = append(u_all, B)
    u_all = append(u_all, F)
    u_all = append(u_vall, zF)
    u_all = append(u_all, F_0)

    xprime = col_cstr_model(t,X,u_all)
    return xprime
