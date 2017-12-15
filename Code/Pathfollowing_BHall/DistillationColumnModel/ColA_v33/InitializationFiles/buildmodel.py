#!/opt/local/bin/python
# -*- encoding: ascii -*-
"""
    @purpose: Creates objective function and constraints for Distillation column
    A model and CSTR
    @author: Brittany Hall
    @date: 11.10.2017
    @version: 0.1
    @updates:
"""
from casadi import *
from numpy import divide, multiply, zeros, array

def buildmodel(u, params):
    #Unpacking model parameters
    NT = params['dist']['NT']                                 #number of stages
    NF = params['dist']['NF']                          #stage where feed enters
    alpha = params['dist']['alpha'] #relative volatility
    Muw = params['dist']['Muw'] #nominal liquid hold ups
    taul = params['dist']['taul'] #time constant for liquid dynamics
    F = params['dist']['F'] #nominal distillation feed flowrate
    qF = params['dist']['qF'] #nominal distillation feed liquid fraction
    L0 = params['dist']['L0'] #nominal reflux flow
    L0b = params['dist']['L0b'] #nominal liquid flow below feed
    F_0 = params['dist']['F_0'] #nominal CSTR feed flowrate
    zF = params['dist']['zF'] #nominal feed composition

    #Inputs and disturbances
    LT = u[2*NT+2]                                                      #Reflux
    VB = u[2*NT+3]                                                      #Boilup
    F = u[2*NT+4]                                                     #Feedrate
    D = u[2*NT+5]                                                   #Distillate
    B = u[2*NT+6]                                                      #Bottoms

    """
    The Model
    """
    #Objective function
    pf = params['price']['pf']
    pV = params['price']['pV']
    pB = params['price']['pB']
    pD = params['price']['pD']
    J = pf*F_0 + pV*VB - pB*B - pD*D

    #Vapor and Liquid flowrates, composition, and holdups
    y = SX.zeros(NT-1)
    V = SX.zeros(NT-1)
    L = SX.zeros(NT)
    dMdt = SX.zeros(NT+1)
    dMxdt = SX.zeros(NT+1)
    for i in range(0,NT-1):
        y[i] = SX.sym('y_'+str(i+1),1)
        V[i] = SX.sym('V_'+str(i+1),1)
        L[i] = SX.sym('L_'+str(i+1),1)
        dMdt[i] = SX.sym('dMdt_'+str(i+1),1)
        dMxdt[i] = SX.sym('dMxdt_'+str(i+1),1)
    L[NT-1] = SX.sym('L_'+str(NT),1)
    dMdt[NT-1] = SX.sym('dMdt_'+str(NT),1)
    dMxdt[NT-1] = SX.sym('dMxdt_'+str(NT),1)
    dMdt[NT] = SX.sym('dMdt_'+str(NT+1),1)
    dMxdt[NT] = SX.sym('dMxdt_'+str(NT+1),1)

    #Vapor-liquid equilibria
    for i in range(0,NT-1): #don't calculate value for last stage NT
        y[i] = alpha*u[i]/(1+(alpha-1)*u[i])

    #Vapor flows (constant molar flows assumed)
    for i in range(0,NT-1):#don't calculate value for last stage NT
        if i >= NF-1:
            V[i] = VB + (1-qF)*F
        else:
            V[i] = VB

    #Liquid flows
    L[NT-1] = LT #last stage liquid
    for i in range(0,NT-1):#don't calculate value for last stage NT
        if i <= NF-1:
            L[i] = L0b + divide((u[NT+1+i] - Muw),taul)
        else:
            L[i] = L0 + divide((u[NT+1+i] - Muw),taul)

    #Time derivatives for material balances for total holdup and component
    for i in range(1,NT-1):
        dMdt[i] = L[i+1] - L[i] + V[i-1] - V[i]
        dMxdt[i] = multiply(L[i+1], u[i+1,0]) - multiply(L[i], u[i,0]) + multiply(V[i-1], y[i-1]) - multiply(V[i],y[i])

    #Correction for feed stage
    dMdt[NF-1] = dMdt[NF-1] + F
    dMxdt[NF-1] = dMxdt[NF-1] + F*u[NT]

    #Reboiler (assumed to be an equilibrium stage)
    dMdt[0] = L[1] - V[0] - B
    dMxdt[0] = L[1]*u[1] - V[0]*y[0] - B*u[0]

    #Total condenser (not an equilbrium stage)
    dMdt[NT-1] = V[NT-2] - LT - D
    dMxdt[NT-1] = V[NT-2]*y[NT-2] - LT*u[NT-1] - D*u[NT-1]

    #Compute the derivative for the mole fractions d(Mx) = xdM + Mdx
    ceq = SX.zeros(2*NT+2)
    for i in range(0,2*NT+2):
        ceq[i] = SX.sym('ceq_'+str(i+1),1)

    #CSTR model
    k1 = params['cstr']['k1']
    dMdt[NT] = F_0 + D - F
    dMxdt[NT] = F_0*zF + D*u[NT-1] - F*u[NT] - k1*u[2*NT+1]*u[NT]

    for i in range(0,NT+1):
        ceq[i] = dMxdt[i]
    
    for i in range(0,NT+1):
        ceq[NT+1+i] = dMdt[i]

    #Constraint bounds
    lbx = params['bounds']['lbx']
    ubx = params['bounds']['ubx']
    lbg = params['bounds']['lbg']
    ubg = params['bounds']['ubg']

    return J, ceq, lbx, ubx, lbg, ubg
