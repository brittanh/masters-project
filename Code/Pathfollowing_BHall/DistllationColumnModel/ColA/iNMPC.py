#!/opt/local/bin/python
# -*- encoding: ascii -*-
"""
    @purpose: Ideal Nonlinear Model Predictive Control (iNMPC)
    @author: Brittany Hall
    @date: 11.10.2017
    @version: 0.1
    @updates:
"""
from numpy import size, zeros, append
import scipy.io as spio
from compObjFn import *
from solveOpt import *
from plotStates import *

def iNMPC(optProblem, system, MPCit, N, T, tmeasure, xmeasure, u0, params):
    
    #Unpacking required parameters
    NT = params['dist']['NT']
    
    #Constructing empty arrays for later use
    Tall = []
    Xall = zeros((MPCit, size(xmeasure, axis = 0)))
    Uall = zeros((MPCit, size(u0, axis = 0)))
    xmeasureAll = []
    uAll = []
    runtime = []
    u_nlp_opt = []
    x_nlp_opt = []

    #NMPC iteration
    iter = 1

    #Load in noise data
    data = spio.loadmat('noise1pct.mat', squeeze_me = True)
    noise = data['noise']

    while (iter <= MPCit):
        print "-----------------------------------------------------------\n"
        print "MPC iteration: %d \n" %(iter)

        #Obtaining new initial value
        def measureInitVal(tmeasure, xmeasure):
            t0 = tmeasure
            x0 = xmeasure
            return t0, x0
        t0,x0 = measureInitVal(tmeasure, xmeasure)

        #Measurement noise
        n_M = noise[:,iter-1]                                     #Holdup noise
        n_X = zeros((NT+1,1))                              #Concentration noise
        measure_noise = append(n_X, n_M)
        x0_measure = x0 + measure_noise     #Adding measurement noise to states

        #Solving NLP
        primalNLP, _, lb, ub, _, params, elapsedtime = solveOpt(optProblem, x0,
                                                     u0, N, x0_measure, params)

        #Re-arrange NLP solutions
        #(turning vectors into matrices to make easier to plot)
        u_nlp_opt, x_nlp_opt = plotStates(primalNLP, lb, ub, N, params)

        #Save open loop solution for error computation
        z1 = x_nlp_opt[0:nx,4]
        #Record information
        Tall = append(Tall, t0)
        Xall[iter+1,:] = transpose(x0)
        Uall[iter+1,:] = u0[:,0]

        #Applying control to process with optimized control
        def dynamic(system, T, t0, x0, u0):
            x = system(t0, x0, u0, T)
            x_intermediate = append(x0, x)
            t_intermediate = append(t0, t0+T)
            return x, t_intermediate, x_intermediate
        
        def applyControl(system, T, t0, x0, u0):
            xapplied = dynamic(system, T, t0, x0, u0[:,0])
            tapplied = t0 + T
            return xapplied, tapplied
        
        x0 = xmeasure
        tmeasure, xmeasure = applyControl(system, T, t0, x0, u_nlp_opt)
        print xmeasure
        raw_input()

        #Using actual state
        ObjVal = []
        ObjVal[iter] = compObjFn(u_nlp_opt[:,0], xmeasure)

        #Storing Output Variables
        xmeasureAll = append(xmeasureAll, xmeasure)
        uAll = append(uAll, u_nlp_opt[:,0])
        runtime = append(runtime, elapsedtime)

        def shiftHorizon(uopt):
            u0 = [u[:,1:size(u,axis = 1)], u[:,size(u,axis = 1)] ]
            return u0

        u0 = shiftHorizon(u_nlp_opt)

        iter += 1
            
    xmeasureAll = xmeasureAll.setshape(84, mpciterations)
    
    return Tall, xmeasureAll, uAll, ObjVal, primalNLP, params, runtime
