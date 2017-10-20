#!/opt/local/bin/python
# -*- encoding: ascii -*-
"""
    @purpose:
    @author: Brittany Hall
    @date: 06.10.2017
    @version: 0.1
    @updates:
"""
from numpy import reshape, tile
import scipy.io as spio
#user made functions
from optProblem import *
from system import *
from pfNMPC import *
from iNMPC import *
from params import *

#MPC iterations
MPCit = 150
#Prediction Horizon
N = 30
#Sampling time
T = 1                                                                    #[min]

#Loading in initial data (different initial conditions)
data = spio.loadmat('Xinit29.mat', squeeze_me = True)
Xinit = data['Xinit29']
u0 = Xinit[84:89]                                               #Initial inputs
u0 = u0.reshape(len(u0),1)
u0 = tile(u0,N)
tmeasure = 0.0                                                      #Start time
xmeasure = Xinit[0:84]                                          #Initial states
Uf = 0.3                                               #Feed rate to CSTR (F_0)
params['dist']['F_0'] = Uf

#Applying ideal NMPC
_, xmeasureAll, uAll, obj, optRes, params, runtime = iNMPC(optProblem, system, MPCit, N, T, tmeasure, xmeasure, u0, params)

print "iNMPC finished \n"

#Applying path-following NMPC
#_, xmeasureAll_pf, uAll_pf, obj_pf, optRest_pf, params_pf, runtime_pf = pfNMPC(optProblem, system, MPCit, N, T, tmeasure, xmeasure, u0, varargin)


