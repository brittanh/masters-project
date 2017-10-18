#!/opt/local/bin/python
# -*- encoding: ascii -*-
"""
    @purpose:
    @author: Brittany Hall
    @date: 06.10.2017
    @version: 0.1
    @updates:
"""
from numpy import *
import scipy.io as spio
#user made functions
from optProblem import *
from system import *
from pfNMPC import *
from iNMPC import *
from params import *

#Global variables
global N, params

#MPC iterations
MPCit = 150
#Prediction Horizeon
N = 30
#Sampling time
T = 1                                                                    #[min]

#Loading in initial data (different initial conditions)
data = spio.loadmat('Xinit29.mat', squeeze_me = True)
Xinit = data['Xinit29']
u0 = Xinit[85:89]                                               #initial inputs
tmeasure = 0.0                                                      #start time
xmeasure = Xinit[0:84]                                          #initial states

#Applying ideal NMPC
_, xmeasureAll, uAll, obj, optRes, params, runtime = iNMPC(optProblem, system, MPCit, N, T, tmeasure, xmeasure, u0, params)

#Applying path-following NMPC
#_, xmeasureAll_pf, uAll_pf, obj_pf, optRest_pf, params_pf, runtime_pf = pfNMPC(optProblem, system, MPCit, N, T, tmeasure, xmeasure, u0, varargin)


