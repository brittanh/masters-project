#!/opt/local/bin/python
# -*- encoding: ascii -*-
"""
    @purpose: Path- following based Nonlinear Model Predictive Control (pfNMPC)
    @author: Brittany Hall
    @date: 07.10.2017
    @version: 0.1
    @updates:
"""
from solveOpt import solveOpt
from scipy.io import savemat, loadmat
from plotStates import plotStates
from ColCSTR_pf import ColCSTR_pf
from predictor_corrector import predictor_corrector 
from numpy import size, zeros, append, array, hstack, vstack
from compObjFn import *

def pfNMPC(optProblem, system, MPCit, N, T, tmeasure, xmeasure, u0, params):
    
    NT = params['dist']['NT']
    #Dimension of state and input
    nx = size(xmeasure) #Elements in state
    nu = size(u0, axis = 0) #Size of inputs
    #Constructing empty arrays for later use
    Tall = []
    Xall = zeros((MPCit, xmeasure.shape[0]))
    Uall = zeros((MPCit, u0.shape[0]))
    ObjVal = {}
    ObjVal['econ'] = []
    ObjVal['reg'] = []
    xmeasureAll = []
    uAll = []
    runtimepf = []
    u_pf_opt = []
    x_pf_opt = []

    #starting NMPC iteration
    iter = 1
    z1 = xmeasure
    #loading in noise data
    data = loadmat('noise1pct.mat', squeeze_me = True)
    noise = data['noise']
    runtime_pf = 0
    
    while (iter <= MPCit):
        print("------------------------------------\n")
        print("MPC iteration: %d\n" %iter)

        #Obtaining new initial value
        def measureInitVal(tmeasure, xmeasure):
            t0 = tmeasure
            x0 = xmeasure
            return t0, x0
        t0,x0 = measureInitVal(tmeasure, xmeasure)

        #adding measurement noise
        n_M = noise[:,iter-1]                 #Holdup noise
        n_X = zeros((NT+1,1))           #Concentration noise
        measure_noise = append(n_X, n_M)
        x0_measure = x0 + measure_noise  #Add measmt noise to states

        #advanced-step NMPC
        primalNLP, dualNLP, lb, ub, objVal, params,_ = solveOpt(optProblem, x0, u0, N, z1, params)
        
        #re-arrange NLP solutions
        _, x_nlp_opt = plotStates(primalNLP, lb, ub, N, params)

        p_init = primalNLP[0:nx]
        p_final = x0_measure
        xstart = primalNLP
        ystart = dualNLP
        
        delta_t = 0.5                                               #Step size
        lb_init = lb
        ub_init = ub

        #NLP sensitivity (predictor-corrector)
        primalPF, _, elapsedqp = predictor_corrector(lambda p:ColCSTR_pf(p),
            p_init, p_final, xstart, ystart, delta_t, lb_init, ub_init, 0, N)
            
        #Formatting variables for plotting
        u_pf_opt, x_pf_opt = plotStates(primalPF, lb, ub, N, params)
        z1 = x_pf_opt[0:nx,5]
    
        #Store output variables
        Tall = append(Tall, t0)
        Xall[iter-1,:] = transpose(x0)
        Uall[iter-1,:] = u0[:,0]
        
        #Apply control to process with optimized control from pf-algorithm
        x0 = xmeasure #from online step
        
        def dynamic(system, T, t0, x0, u0):
            x = system(t0, x0, u0, T)
            x_intermediate = vstack((x0, x))
            t_intermediate = hstack((t0, t0+T))
            return x, t_intermediate, x_intermediate
        
        def applyControl(system, T, t0, x0, u0):
            xapplied, _, _ = dynamic(system, T, t0, x0, u0[:,0])
            tapplied = t0 + T
            return tapplied, xapplied
        
        tmeasure, xmeasure = applyControl(system, T, t0, x0, u_pf_opt)
     
        #Using actual states to compute the objective function values
        Jobj = compObjFn(u_pf_opt[:,0],xmeasure)
        
        #Storing Output Variables
        ObjVal['econ'].append(float(Jobj['econ'][0]))
        ObjVal['reg'].append(float(Jobj['reg'][0]))
        
        #Collect variables
        xmeasureAll = append(xmeasureAll, xmeasure)
        uAll = append(uAll, u_pf_opt[:,0])
        runtime_pf = append(runtime_pf, elapsedqp)
        
        #Prepare restart
        def shiftHorizon(u):
            u0 = hstack((u[:,1:u.shape[1]], u[:,u.shape[1]-1]))
            return u0
        
        u0 = shiftHorizon(u_pf_opt)

        iter += 1

    xmeasureAll = reshape(xmeasureAll,(xmeasureAll.shape[0],1))
    xmeasureAll = reshape(xmeasureAll, (2*NT+2, MPCit))
    xmeasureAll = array(xmeasureAll)

    ObjReg = array(ObjVal['reg'])
    ObjEcon = array(ObjVal['econ'])

    pathfollowing = {
                'pf':{
                    'xmeasureAll': xmeasureAll,
                    'uAll': uAll,
                    'ObjReg': ObjReg,
                    'ObjEcon': ObjEcon,
                    'T': T,
                    'mpciterations': MPCit
                    }
                }

    savemat('pfNMPC.mat',pathfollowing)            #saving pfNMPC results

    return Tall, xmeasureAll, uAll, ObjVal, primalPF, params, runtime_pf
