#!/opt/local/bin/python
# -*- encoding: ascii -*-
"""
    @purpose: Nonlinear distillation column model at steady state with NT-1
    theoretical stages including a reboiler (stage 1) plus a total condenser
    (stage NT).
    The model is based on column A in Skogestad and Postlethwaite (1996)
    @author: Brittany Hall
    @date: 04.10.2017
    @version: 0.1
    @updates:
"""
from numpy import *
def col_ss(t, X, U):
    #Column Information
    """
        Inputs: t - time [min]
                X - States, the first 41 states are compositions of light
                    component A with reboiler/bottom stage as X(0) and
                    condenser as X(40). X(41) is the holdup in the
                    reboiler/bottom stage and X(81) is the hold-up in condenser
                U[0] - reflux L
                U[1] - boilup V
                U[2] - top or distillate product flow D
                U[3] - bottom product flow B
                U[4] - feed rate F
                U[5] - feed composition zF
                U[6] - feed liquid fraction qF
        
        Outputs: xss- vector of steady state states
    """
    #NT = 82
    """
    Column specific properties
    Need to be changed if another column is used
    """
    NT = 41                                         #Number of stages in column
    NF = 21                       #Location of feed stage (counted from bottom)
    alpha = 1.5                                            #Relative volatility
    MO = zeros(NT)
    MO[0] = 0.5                                 #Nominal reboiler holdup [kmol]
    MO[1:NT-1] = 0.5    #Nominal stage (tray) holdup [kmol]
    MO[NT] = 0.5                               #Nominal condenser holdup [kmol]

    #Liquid flow dynamics (linear model)
    taul = 0.063                       #Time constant for liquid dynamics [min]
    FO = 1                                        #Nominal feed rate [kmol/min]
    qFO = 1                                 #Nominal fraction of liquid in feed
    L0 = 2.70629                  #Nominal reflux flow (from steady-state data)
    L0b = L0 + qF0*F0                #Nominal liquid flow below feed [kmol/min]
    lam = 0                    #Effect of vapor flow on liquid flow (K2 effect)
    #V0 = 3.20629
    #V0t = V0 + (1-qF0)*F0                                 #Nominal vapor flows\
                                                          (only needed if lam \
                                                           (lambda) is nonzero)
    #Dividing the states
    x = transpose(X[0:NT-1])             #Liquid composition from bottom to top
    M = transpose(X(NT+1:2*NT))              #Liquid hold up from bottom to top

    #Inputs and Disturbances
    LT = U[0]                                                  #Reflux flowrate
    VB = U[1]                                                  #Boilup flowrate
    D = U[2]                                               #Distillate flowrate
    B = U[3]                                                  #Bottoms flowrate
    F = U[4]                                                     #Feed flowrate
    zF = U[5]                                                 #Feed composition
    qF = U[6]                                             #Feed liquid fraction

    """
    Model Development
    """
    #Vapor-liquid equilibria
    y = zeros(len(x))
    y = [(alpha*x[i])/(1+(alpha-1)*x[i]) for i in range(0,NT-1)]

    #Vapor flows (assuming constant molar flow)
    V = zeros(len(x))
    V[0:NT-1] = VB*ones((1,NT-1))
    V[NF:NT-1] = V[NF:NT-1] + (1-qF)*F

    #Liquid Flows (linearized tray hydraulics)
    #only needed if lambda does NOT equal 0
#    L = zeros(len(x))
#    L[1:NF] = L0b +  (M[1:NF]-M0[1:NF])/taul + lam*(V[0:NF-1]-V0)
#    L[NF+1:NT-1] = L0 + (M[1:NF]-M0[1:NF])/taul + lam*(V[0:NF-1]-V0t)

    #Steady state
    V[0] = VB
    L[NT] = LT

    #Total Condenser
    V[NT-1] = L[NT] + D
    x[NT] = V[NT-1]*y[NT-1]/(LT+D)

    #Reboiler
    L[1] = V[0] + B
    x[0] = (L[1]*x[1]-V[0]*y[0])/B

    for i in range(NT-1,1,-1):
        if i == NF:
            L[i] = L[i+1] + V[i-1]-V[i] + F
            x[i] = (L[i+1]*x[i+1] + V[i-1]*y[i-1]-V[i]*y[i] + F*zF)/L[i]
        else:
            L[i] = L[i+1] + V[i-1] - V[i]
            x[i] = (L[i+1]*x[i+1] + V[i-1]*y[i-1]-V[i]*y[i])/L[i]

    xss = append(transpose(x), transpose(M))
    return xss
