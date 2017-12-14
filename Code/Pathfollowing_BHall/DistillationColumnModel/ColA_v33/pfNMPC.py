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
import scipy.io as spio
from plotStates import plotStates
from ColCSTR_pf import ColCSTR_pf
from predictor_corrector import predictor_corrector 
from numpy import size, zeros, append, array

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
    data = spio.loadmat('noise1pct.mat', squeeze_me = True)
    noise = data['noise']
  
    while (iter <= MPCit):
        print("--------------------------------------------------\n")
        print('MPC iteration: %d\n' %iter)

        #Obtaining new initial value
        def measureInitVal(tmeasure, xmeasure):
            t0 = tmeasure
            x0 = xmeasure
            return t0, x0
        t0,x0 = measureInitVal(tmeasure, xmeasure)

        #adding measurement noise
        n_M = noise[:,iter-1]                          #Holdup noise
        n_X = zeros((NT+1,1))                   #Concentration noise
        measure_noise = append(n_X, n_M)
        x0_measure = x0 + measure_noise  #Add measmt noise to states

        #advanced-step NMPC
        primalNLP,dualNLP,lb,ub,objVal,params,_=solveOpt(optProblem,
                                                x0,u0,N,z1,params)
        
        #re-arrange NLP solutions
        _, x_nlp_opt = plotStates(primalNLP, lb, ub, N, params)

        p_init = primalNLP[0:nx]
        p_final = x0_measure
        xstart = primalNLP
        ystart = dualNLP
        
        delta_t = 0.5                                  #Step size
        lb_init = lb
        ub_init = ub

        #NLP sensitivity (predictor-corrector)
        primalPF, _,elapsedqp=predictor_corrector(lambda p:ColCSTR_pf(p),
            p_init,p_final,xstart,ystart,delta_t,lb_init,ub_init,0,N)

        runtime_pf = append(runtime_pf, elapsedqp)
            
    return Tall,xmeasureAll,uAll,ObjVal,primalPF,params,runtime_pf
