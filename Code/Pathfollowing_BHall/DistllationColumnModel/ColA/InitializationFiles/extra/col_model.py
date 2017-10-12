#!/opt/local/bin/python
# -*- encoding: ascii -*-
"""
    @purpose: Nonlinear distillation column model with NT-1 theoretical
    stages including a reboiler (stage 1) plus a total condenser (stage NT).
    The model is based on column A in Skogestad and Postlethwaite (1996).
    @author: Brittany Hall
    @date: 05.10.2017
    @version: 0.1
    @updates:
"""
from numpy import *
def col_model(t, X, U):
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
        
        Outputs: xprime - vector of time derivative all states
    """
  #-------------------------------------------------------------------------#
    """
        Column specific properties
        Need to be changed if another column is used
    """
    NT = 41                                         #Number of stages in column
    NF = 21                                             #Location of feed stage
    alpha = 1.5                                            #Relative volatility
    MO = zeros(NT+1)
    MO[0] = 0.5                                 #Nominal reboiler holdup [kmol]
    MO[1:NT-1] = 0.5                        #Nominal stage (tray) holdup [kmol]
    MO[NT] = 0.5                               #Nominal condenser holdup [kmol]
    MO[NT+1] =0.5

    #Data for linearized Liquid flow dynamics\
    (does not apply to reboiler and condenser)
    taul = 0.063                       #Time constant for liquid dynamics [min]
    FO = 1                                        #Nominal feed rate [kmol/min]
    qFO = 1                                 #Nominal fraction of liquid in feed
    L0 = 2.70629                  #Nominal reflux flow (from steady-state data)
    L0b = L0 + qF0*F0                #Nominal liquid flow below feed [kmol/min]
    lam = 0                    #Effect of vapor flow on liquid flow (K2 effect)
    V0 = 3.20629                                               #Nominal boilup
    #V0t = V0 + (1-qF0)*F0                                #Nominal vapor flows\
                                                          (only needed if lam \
                                                           (lambda) is nonzero)
    #-------------------------------------------------------------------------#
    
    #Dividing the states
    x = transpose(X[0:NT])               #Liquid composition from bottom to top
    M = transpose(X(NT+1:2*NT))              #Liquid hold up from bottom to top

    #Inputs and Disturbances
    LT = U[0]                                                  #Reflux flowrate
    VB = U[1]                                                  #Boilup flowrate
    D = U[2]                                               #Distillate flowrate
    B = U[3]                                                  #Bottoms flowrate
    F = U[4]                                        #Distillation feed flowrate
    zF_0 = U[5]                                          #CSTR Feed composition
    qF = U[6]                                                   #CSTR flowrate
    """
    Model Development
    """
    #Vapor-liquid equilibria
    y = zeros(len(x))
    y[0:NT-1] = alpha*x[0:NT-1]/(1+(alpha-1)*x[0:NT-1])
    
    #Vapor flows (assuming constant molar flow)
    V = zeros(len(x))
    V[0:NT-1] = VB*ones((1,NT-1))
    V[NF:NT-1] = V[NF:NT-1] + (1-qF)*F

    #Liquid flows (assuming linearized tray hydraulics)
    L = zeros(len(x))
    L[1:NF] = L0b +  (M[1:NF]-M0[1:NF])/taul
    L[NF+1:NT-1] = L0 + (M[1:NF]-M0[1:NF])/taul
    L[NT] = LT

    """
    Time Derivatives of material balances for total holdup and component holdup
    """
    #Column
    dMdt[1:NT-1] = L[i+1] - L[i] + V[i-1] - V[i]
    dMxdt[1:NT-1] = L[i+1]*x[i+1]-L[i]*x[i]+V[i-1]*y[i-1]-V[i]*y[i]

    #Correction for feed at feed stage (feed is assumed to be mixed)
    dMdt[NF] = dMdt[NF] + F
    dMxdt[NF] = dMxdt[NF] + zF*F

    #Reboiler (assumed to be an equilibrium stage)
    dMdt[0] = L[1]-V[0]-B
    dMxdt[0] = L[1]*x[1]-V[0]*y[0] -B*x[0]

    #Total condensor (not an equilibrium stage)
    dMdt[NT] = V[NT-1] - LT - D
    dMxdt[NT] = V[NT-1]*y[NT-1]-LT*x[NT]-D*x[NT]

    #CSTR Model (inputs F_0 z_F0)
    k1 = 34.1/60.0                                           #Reaction rate [?]
    dMdt[NT+1] = F_0 + D - F
    dMxdt[NT+1] = F_0*z_F0 + D*x[NT] - F*x[NT+1] -k1*M[NT+1]*x[NT+1]
    
    #Calculating the derivative of the mole fractions
    dxdt[0:NT] = (dMxdt[0:NT]-x[0:NT]*dMdt[0:NT])/M[0:NT]

    xprime = append(transpose(dxdt),transpose(dMdt))
    return xprime

