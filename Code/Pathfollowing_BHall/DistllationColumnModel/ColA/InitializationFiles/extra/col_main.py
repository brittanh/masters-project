#!/opt/local/bin/python
# -*- encoding: ascii -*-
"""
    @purpose: Subroutine for simulation of distillation column using ode solver.
              Calls the distillation model col_model.py
    @author: Brittany Hall
    @date: 05.10.2017
    @version: 0.1
    @updates:
"""
from numpy import *
from scipy import integrate

def col_main(t,X):
    #Inputs and disturbances
    LT=2.70629                                                          #Reflux
    VB=3.20629                                                          #Boilup
    D=0.5                                                           #Distillate
    B=0.5                                                              #Bottoms
    F=1.0 + 0.00                                                      #Feedrate
    zF=0.5                                                    #Feed composition
    qF=1.0                                                #Feed liquid fraction

    #Storing all inputs and disturbances
    U[0] = LT
    U[1] = VB
    U[2] = D
    U[3] = B
    U[4] = F
    U[5] = zF
    U[6] = qF

    xprime = colamod(t,X,U)
    return xprime

