#!/opt/local/bin/python
# -*- encoding: ascii -*-
"""
    @purpose: solving optimal control problem
    @author: Brittany Hall
    @date: 07.10.2017
    @version: 0.1
    @updates:
"""
from casadi import *
from numpy import transpose, shape, zeros
from optProblem import *
import time
from nlp_solve import *

def solveOpt(optProblem, system, N, t0, x0, u0, T, iter, u_pf_opt, x_pf_opt, z1):
    
    global ns
    x0_measure = z1
    x = zeros((N+1,84))
    x[0,:] = transpose(x0)
    for k in range(0,N):
        x[k+1,:] = transpose(x0)

    J, g, w0, w, lbg, ubg, lbw, ubw, params = optProblem(x, u0, N, x0_measure)

    #Solving the NLP
    prob = {'f': J, 'x': vertcat(w),'g': vertcat(g)}
    options = {}
    tic = time.clock()
    startnlp = tic
    sol = nlp_solve(prob, options)
    toc = time.clock()
    elapsednlp = toc - tic
    print('IPOPT solver run time = %f\n', elapsednlp)

    u = sol['x']
    lam = {}
    lam['lam_g'] = sol['lam_g']
    lam['lam_x'] = sol['lam_x']
    objVal = sol['f']

    return u, lam, lbw, ubw, objVal, params, elapsednlp
