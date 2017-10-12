#!/opt/local/bin/python
# -*- encoding: ascii -*-
"""
    @purpose: Function used for simulation of distillation column with
    LV-configuration
    @author: Brittany Hall
    @date: 05.10.2017
    @version: 0.1
    @updates:
"""
from numpy import *
def col_LV(X):
    #Column Information
    """
        Inputs: reflux (LT) and boilup (VB)
        Disturbances: Feedrate (F) and feed composition (zF)
                      These are changed by changing the corresponding values
                      'col_model.m'.
        
        Outputs: x-liquid composition and liquid holdup for stages 1 through NT
    """
  #-------------------------------------------------------------------------#
    """
        Column specific properties
        Need to be changed if another column is used
    """
    global u         #Make perturbed inputs and disturbances available to model
    
    #P-controllers for control of reboiler and condenser hold up
    KcB = 10                                          #controller gain reboiler
    KcD = 10                                        #controller gain distillate
    MDs = 0.5                                        #Nominal holdup distillate
    MBs = 0.5                                            #Nominal holdup bottom
    Ds = 0.5                                           #Nominal distillate flow
    Bs = 0.5                                               #Nominal bottom flow
    MB = X(NT+1)                                        #Actual reboiler holdup
    MD = X(2*NT)                                       #Actual condenser holdup
    D = Ds + (MD-MDs)*KcD                           #Controlled distillate flow
    B = Bs + (MB-MBs)*KcB                              #Controlled bottoms flow

    #Storing all inputs and disturbances
    u_all = u[1:2]
    u_all = append(u_all,D)
    u_all = append(u_all,B)
    u_all = append(u_all,u[3:5])

    t= array([])
    xprime = col_model(t,X,u_all)



