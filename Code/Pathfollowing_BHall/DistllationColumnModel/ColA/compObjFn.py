#!/opt/local/bin/python
# -*- encoding: ascii -*-
"""
    @purpose: Computing objective function values
    @author: Brittany Hall
    @date: 11.10.2017
    @version: 0.1
    @updates:
"""
from numpy import size, transpose, multiply
import scipy.io as spio
from itPredHorizon import *
from params import *

def compObjFn(uOpt,xActual):
    
    #prices
    pf = params['price']['pf']
    pV = params['price']['pV']
    pB = params['price']['pB']
    pD = params['price']['pD']
    
    #Setpoints
    F_0 = params['dist']['F_0']
    
    #Steady-state values
    data = spio.loadmat('CstrDistXinit.mat', squeeze_me=True)
    Xinit = data['Xinit']
    
    xs = Xinit[0:84]
    us = Xinit[84:]
    nx = size(xs, axis = 0)
    nu = size(us, axis = 0)

    #Loading in objective function weights
    data = spio.loadmat('Q.mat', squeeze_me = True)
    Qmax = data['Q']
    c1 = -0.05 #noise
    lss = -0.256905910000000 + c1 #steady state objective function value
    
    #Defining objective function
    Jecon = pf*F_0 + pV*uOpt[1] - pB*uOpt[4] - pD*uOpt[3]
    Jcontrol = mtimes(transpose(multiply(Qmax[nx:nx+nu], uOpt - us)), (uOpt - us))
    Jstate = mtimes(transpose(multiply(Qmax[0:nx],(xActual -xs))),
                    (xActual - xs))

    J = Jecon + Jcontrol + Jstate -lss
    
    print('----------------------------------------------------------------\n')
    print("Jecon: %f,\n Jcontrol: %f, \n Jstate: %f, \n" %(Jecon, Jcontrol, Jstate))
    Jobj = {}
    Jobj['reg'] = J
    Jobj['econ'] = Jecon
    
    return Jobj
