#!/opt/local/bin/python
# -*- encoding: ascii -*-
"""
    @purpose: CSTR model (stage NT+1) with a first order reaction (A-> B) plus
    nonlinear distillation column model with NT-1 theoretical stages including
    a reboiler (stage 1) plus a total condenser (stage NT).
    The model is based on column A in Skogestad and Postlethwaite (1996).
    @author: Brittany Hall
    @date: 06.10.2017
    @version: 0.1
    @updates:
"""
from casadi import *
from numpy import array, Infinity

def ColCSTR_model(U):

    #Model parameters
    #==================Column Dependent Properties============================#
    NC = 2                                                #Number of components
    NF = 21                                           #Stage number feed enters
    NT = 41                                         #Number of stages in column
    qF = 1.0                                              #Feed liquid fraction
    alpha = 1.5                                            #Relative volatility
    zF0 = 1.0                         #Feed composition to CSTR [mole fraction]
    Muw = 0.5
    F_0 = U                                  #Feed flow rate to CSTR [kmol/min]
    
    #Data for linearized Liquid flow dynamics
    #(does not apply to reboiler and condenser)
    taul = 0.063                       #Time constant for liquid dynamics [min]
    F0 = 1                                        #Nominal feed rate [kmol/min]
    qF0 = 1                                 #Nominal fraction of liquid in feed
    L0 = 2.70629                  #Nominal reflux flow (from steady-state data)
    L0b = L0 + qF0*F0                #Nominal liquid flow below feed [kmol/min]
    lam = 0                    #Effect of vapor flow on liquid flow (K2 effect)
    V0 = 3.20629                                               #Nominal boilup
    V0t = V0 + (1-qF0)*F0                                #Nominal vapor flows
    
    #Save to a parameter file so only have to define model parameters in one place
    
    #=========================================================================#
    
    #States and Control Inputs
    x = SX.sym('x', NT+1, NC-1)                                    #Composition
    M = SX.sym('M',NT+1,1)                                              #Holdup
    states = array([[x],[M]])
    L_T = SX.sym('L_T')                                            #Liquid flow
    V_B = SX.sym('V_B')                                             #Vapor flow
    F = SX.sym('F')                                             #Feed to column
    D = SX.sym('D')                                                 #Distillate
    B = SX.sym('B')                                                     #Bottom
    inputs = array([[L_T],[V_B],[F],[D],[B]])

    t = SX.sym('t')                                                       #Time
    y = SX.sym('y', NT-1, NC-1)                              #Vapor compisition
    Li = SX.sym('Li', NT, 1)                             #Liquid flow on stages
    Vi = SX.sym('Vi', NT, 1)                              #Vapor flow on stages

    dMdt = SX.sym('dMdt', NT+1, 1)                          #Total Molar holdup
    dMxdt = SX.sym('dMxdt', NT+1, NC-1)                  #Component wise holdup
    dxdt = SX.sym('dxdt', NT+1, NC-1)            #Rate of change of composition

    #Vapor flows (assumed constant, no dynamics)
    for i in range(0,NT-2):
        Vi[i] = V_B
        if i >= NF:
            Vi[i] = Vi[i] + (1-qF)*F
    Vi[NT] = Infinity

    #Liquid flows (Wier formula)
    Li[0] = Infinity

    for i in range(1,NT):
        if i <= NF:
            Li[i] = L0b + (M[i]-Muw)/taul
        else:
            Li[i] = L0 + (M[i]-Muw)/taul
    #Top tray liquid
    Li[NT] = L_T

    #Vapor Liquid equilibrium
    for i in range(0,NT-2):
        for j in range(0,NC-1):
            y[i,j] = (x[i,j]*alpha)/(1+(alpha-1)*x[i,j])

    #Partial Reboiler
    dMdt[0] = Li[1] -Vi[0]-B
    for i in range(0,NC-1):
        dMxdt[1,i] = Li[1]*x[1,i] -Vi[0]*y[1,i]-B*x[1,j]

    #Stripping and Enrichment sections
    for i in range(1,NT-2):
        dMdt[i] = Li[i+1] - Li[i] + Vi[i-1]-V[i]
        for j in range(0,NC-1):
            dMxdt[i,j] = Li[i+1]*x[i+1,j] - Li[i]*x[i,j] + Vi[i-1]*y[i-1,j] -Vi[i]*y[i,j]

    #Correction for feed stage
    dMdt[NF] = dMdt[NF] + F
    for j in range(0, NC-1):
        dMxdt[NT,j] = dMxdt[NF, j] + F*x[NT+1]

    #Total Condenser
    dMdt[NT] = Vi[NT-1] - Li[NT] - D
    for j in range(0,NC-1):
        dMxdt[NT,j] = Vi[NT-1]*y[NT-1,j] - Li[NT]*x[NT,j] - D*x[NT,j]

    #CSTR Model
    k1 = 34.1/60.0                                      #Reaction rate constant
    dMdt[NT+1] = F_0 + D - F
    for j in range(0,NC-1):
        dMxdt[NT+1,j] = F_0*z_F0[j,1] + D*x[NT,j] - F*x[NT+1,j] - k1*M[NT+1]*x[NT+1,j]

    for i in range(0, NT+1):
        for j in range(0, NC-1):
            dxdt[i,j] = (dMxdt[i,j]-x[i,j]*dMdt[i])/M[i]

    xdot = array([[dxdt],[dMdt]])

    return t, states, xdot, inputs 
