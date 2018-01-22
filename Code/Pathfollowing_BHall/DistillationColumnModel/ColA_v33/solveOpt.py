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
from numpy import transpose, shape, zeros, savetxt
import numpy
numpy.set_printoptions(threshold=numpy.nan)
from optProblem import *
import time
from nlp_solve import *
from collections import *

def solveOpt(optProblem, x0, u0, N, z1, params):
    
    x0_measure = z1
    x = zeros((N+1,84))
    x[0,:] = transpose(x0)
    for k in range(0,N):
        x[k+1,:] = transpose(x0)

    J, g, w0, w, lbg, ubg, lbw, ubw, params = optProblem(x, u0, x0_measure, N, params)

    #Solving the NLP
    NLP = {'x': w, 'f': J, 'g': g}
    options = {}
    tic = time.clock()
    startnlp = tic
    sol = nlp_solve(NLP, options, w0, lbw, ubw, lbg, ubg)
    toc = time.clock()
    elapsednlp = toc - tic
    print "IPOPT solver run time = %f\n" %elapsednlp

    u = sol['x']
    lam = {}
    lam['lam_g'] = sol['lam_g']
    lam['lam_x'] = sol['lam_x']
    objVal = sol['f']

    return u, lam, lbw, ubw, objVal, params, elapsednlp, w
