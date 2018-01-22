#!/opt/local/bin/python
# -*- encoding: ascii -*-
"""
    @purpose: CSTR model (stage NT+1) with a first order
    reaction (A-> B) plus nonlinear distillation column model
    with NT-1 theoretical stages including a reboiler (stage 1)
    plus a total condenser (stage NT).
    The model is based on column A in Skogestad and Postlethwaite (1996).
    @author: Brittany Hall
    @date: 05.10.2017
    @version: 0.1
    @updates:
"""
from casadi import SX
from numpy import *
from params import *

def col_cstr_model(t, X, U):
    #Column Information
    """
        Inputs: t - time [min]
                X - States, the first 41 states are compositions
                    of light component A with reboiler/bottom
                    stage as X(0) and condenser as X(40). X(41) is
                    the holdup in the reboiler/bottom stage and X(81)
                    is the hold-up in condenser
                U[0] - reflux L
                U[1] - boilup V
                U[2] - top or distillate product flow D
                U[3] - bottom product flow B
                U[4] - feed rate F
                U[5] - feed composition zF
                U[6] - feed rate F0
        
        Outputs: xprime - vector of time derivative all states
    """
    #Unpacking model parameters
    #=============Column Dependent Properties================#
    NC = params['dist']['NC']
    NF = params['dist']['NF']
    NT = params['dist']['NT']
    qF = params['dist']['qF']
    alpha = params['dist']['alpha']
    zF0 = params['dist']['zF']
    M0 = params['dist']['MO']
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
    
    #=========================================================#
    
    #Dividing the states
    #Liquid composition of column plus composition in tank
    x = X[0:NT+1]

    #Liquid hold up from btm to top of col plus hold up in tank
    M = X[NT+1:]

    #Inputs and Disturbances
    LT = U[0]                                  #Reflux flowrate
    VB = U[1]                                  #Boilup flowrate
    D = U[2]                               #Distillate flowrate
    B = U[3]                                  #Bottoms flowrate
    F = U[4]                         #Distillation feed flowrate
    zF_0 = U[5]                           #CSTR Feed composition
    qF = params['dist']['qF']              #Feed liquid fraction
    F_0 = U[6]                                    #CSTR flowrate

    """
    Model Development
    """
    #Vapor-liquid equilibria
    y = []
    for i in range(0,NT-1):
        y.append(alpha*x[i]/(1+(alpha-1)*x[i]))
    
    #Vapor flows (assuming constant molar flow)
    V = []
    for i in range(0,NT-1):
        V.append(VB)
    for i in range(NF,NT-1):
        V[i] = V[i] + (1-qF)*F

    #Liquid flows (assuming linearized tray hydraulics)
    L = []
    L.append(0)
    for i in range(1,NF):
        L.append(L0b +  (M[i]-M0[i])/taul)

    for i in range(NF,NT-1):
        L.append(L0 + (M[i]-M0[i])/taul)

    L.append(LT)

    """
    Time Derivatives of material balances for total
    holdup and component holdup
    """
    #Column
    dMdt = []
    dMdt.append(0)
    dMxdt = []
    dMxdt.append(0)
    for i in range(1,NT-1):
        dMdt.append (L[i+1] - L[i] + V[i-1] - V[i])
        dMxdt.append(L[i+1]*x[i+1]-L[i]*x[i]+V[i-1]*y[i-1]
                     -V[i]*y[i])

    #Correction for feed at feed stage
    dMdt[NF-1] = dMdt[NF-1] + F
    dMxdt[NF-1] = dMxdt[NF-1] + x[NT]*F

    #Reboiler (assumed to be an equilibrium stage)
    dMdt[0] = L[1]-V[0]-B
    dMxdt[0] = L[1]*x[1]-V[0]*y[0] -B*x[0]

    #Total condensor (not an equilibrium stage)
    dMdt.append(V[NT-2] - LT - D)
    dMxdt.append(V[NT-2]*y[NT-2]-LT*x[NT-1]-D*x[NT-1])

    #CSTR Model (inputs F_0 z_F0)
    k1 = params['cstr']['k1']
    dMdt.append(F_0 + D - F)
    dMxdt.append(F_0*zF0[0,0] + D*x[NT-1] - F*x[NT]
                 -k1*M[NT]*x[NT])

    #Calculating the derivative of the mole fractions
    dxdt = []
    for i in range(0,NT+1):
        dxdt.append((dMxdt[i]-x[i]*dMdt[i])/M[i])

    xprime = append(dxdt,dMdt)
    return xprime


