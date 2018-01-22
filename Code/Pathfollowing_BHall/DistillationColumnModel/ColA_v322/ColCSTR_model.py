#!/opt/local/bin/python
# -*- encoding: ascii -*-
"""
    @purpose: CSTR model (stage NT+1) with a first order reaction (A-> B) plus
    nonlinear distillation column model with NT-1 theoretical stages including
    a reboiler (stage 1) plus a total condenser (stage NT).
    The model is based on column A in Skogestad and Postlethwaite (1996).
    @author: Brittany Hall
    @date: 31.10.2017
    @version: 0.2
    @updates: Fixed bug errors on index assignments
"""
from casadi import *
from numpy import array, Infinity

def ColCSTR_model(U,params):

    #Unpacking model parameters
    #===============Column Dependent Properties================#
    NC = params['dist']['NC']
    NF = params['dist']['NF']
    NT = params['dist']['NT']
    qF = params['dist']['qF']
    alpha = params['dist']['alpha']
    zF0 = params['dist']['zF']
    Muw = params['dist']['Muw']
    F_0 = U
    
    #Data for linearized Liquid flow dynamics
    #(does not apply to reboiler and condenser)
    taul = params['dist']['taul']
    F0 = params['dist']['F0']
    qF0 = params['dist']['qF0']
    L0 = params['dist']['L0']
    L0b = L0 + qF0*F0
    lam = params['dist']['lam']
    V0 = params['dist']['V0']
    V0t = V0 + (1-qF0)*F0
    
    #==========================================================#
    
    #States and Control Inputs
    x = SX.sym('x', NT+1, NC-1)                     #Composition
    M = SX.sym('M',NT+1, 1)                              #Holdup
    states = vertcat(x, M)
    L_T = SX.sym('L_T')                             #Liquid flow
    V_B = SX.sym('V_B')                              #Vapor flow
    F = SX.sym('F')                              #Feed to column
    D = SX.sym('D')                                  #Distillate
    B = SX.sym('B')                                      #Bottom
    inputs = vertcat(L_T,V_B)
    inputs = vertcat(inputs,F)
    inputs = vertcat(inputs,D)
    inputs = vertcat(inputs,B)

    t = SX.sym('t')                                         #Time
    y = SX.sym('y', NT-1, NC-1)                #Vapor composition
    Li = SX.sym('Li', NT, 1)               #Liquid flow on stages
    Vi = SX.sym('Vi', NT, 1)                #Vapor flow on stages

    dMdt = SX.sym('dMdt', NT+1, 1)            #Total Molar holdup
    dMxdt = SX.sym('dMxdt', NT+1, NC-1)    #Component wise holdup
    dxdt = SX.sym('dxdt', NT+1, NC-1)     #Rate of change of comp

    #Vapor flows (assumed constant, no dynamics)
    for i in range(1,NT):
        Vi[i-1] = V_B
        if i-1 >= NF:
            Vi[i-1] = Vi[i-1] + (1-qF)*F
    Vi[NT-1] = float('Inf')

    #Liquid flows (Wier formula)
    Li[0] = float('Inf')
    for i in range(1,NT):
        if i <= NF-1:
            Li[i] = L0b + (M[i]-Muw)/taul
        else:
            Li[i] = L0 + (M[i]-Muw)/taul

    #Top tray liquid
    Li[NT-1] = L_T

    #Vapor Liquid equilibrium
    for i in range(0,NT-1):
        for j in range(0,NC-1):
            y[i,j] = (x[i,j]*alpha)/(1+(alpha-1)*x[i,j])

    #Partial Reboiler
    dMdt[0] = Li[1] - Vi[0] - B
    for i in range(0,NC-1):
        dMxdt[0,i] = Li[1]*x[1,i] - Vi[0]*y[0,i] - B*x[0,j]

    #Stripping and Enrichment sections
    for i in range(1,NT-1):
        dMdt[i] = Li[i+1] - Li[i] + Vi[i-1]-Vi[i]
        for j in range(0,NC-1):
            dMxdt[i,j] = Li[i+1]*x[i+1,j] - Li[i]*x[i,j] + Vi[i-1]*y[i-1,j] -Vi[i]*y[i,j]
    #Correction for feed stage
    dMdt[NF-1] = dMdt[NF-1] + F
    for j in range(0, NC-1):
        dMxdt[NF-1,j] = dMxdt[NF-1, j] + F*x[NT]

    #Total Condenser
    dMdt[NT-1] = Vi[NT-2] - Li[NT-1] - D
    for j in range(0,NC-1):
        dMxdt[NT-1,j] = Vi[NT-2]*y[NT-2,j] - Li[NT-1]*x[NT-1,j]- D*x[NT-1,j]

    #CSTR Model
    k1 = params['cstr']['k1']                          
    dMdt[NT] = F_0 + D - F
    for j in range(0,NC-1):
        dMxdt[NT,j] = F_0*zF0[j] + D*x[NT-1,j] - F*x[NT,j]- k1*M[NT]*x[NT,j]

    for i in range(0, NT+1):
        for j in range(0, NC-1):
            dxdt[i,j] = (dMxdt[i,j]-x[i,j]*dMdt[i])/M[i]

    xdot = vertcat(dxdt,dMdt)

    return t, states, xdot, inputs 
