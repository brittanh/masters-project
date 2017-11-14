#!/opt/local/bin/python
# -*- encoding: ascii -*-
"""
    @purpose: Path- following algorithm
    @author: Brittany Hall
    @date: 09.11.2017
    @version: 0.1
    @updates:
"""
from sympy import Symbol
from numpy import zeros
from qp_solve import qp_solve

def pathfollowing(problem, p_init, p_final, x_init, y_init, delta_t, lb_init,
                  ub_init, verbose_level, N):
    
    pp = Symbol('pp')
    theprob = lambda pp: problem(pp)
    prob = theprob(p_init)
    t = 0
    alpha_2 = 0.5

    #initializing some values
    iter = 0
    elapsedqp = 0
    numX = x_init.shape[0]
    x0 = zeros((numX,1))

    if verbose_level:
        print('Solving problem %s \n' %prob['name'])
        print ('Iteration delta_t t Success\n' %(delta_t,t))

    p_0 = p_init
    #Implementing path following algorithm

    while t<= 1:
        #Calculating the step
        tk = t + delta_t
        p_t = (1-tk)*p_0 + tk*p_final
        step = p_t - p_init
        #Update bound constraint
        if lb_init.size != 0:
            lb = lb_init - x_init
            ub = ub_init - x_init
        else:
            lb = []
            ub = []

        #Solving the QP problem
        qp_exit, y, x_qpopt, lam_qpopt, mu_qpopt, qp_run = qp_solve(prob, p_init, x_init, y_init, step, lb_init, ub_init, verbose_level, N)

        elapsedqp += qp_run
