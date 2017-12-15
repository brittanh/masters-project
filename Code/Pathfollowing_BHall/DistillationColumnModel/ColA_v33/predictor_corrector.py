#!/opt/local/bin/python
# -*- encoding: ascii -*-
"""
    @purpose: Predictor corrector
    @author: Brittany Hall
    @date: 08.10.2017
    @version: 0.1
    @updates:
"""
from casadi import *
from qp_solve import *
from numpy import zeros, shape

def predictor_corrector(problem, p_init, p_final, x_init, y_init, delta_t, lb_init, ub_init, verbose_level, N):
    
    p = p_init
    pp = SX.sym('pp')
    theprob = lambda p: problem(pp)
    prob = theprob(p)
    t = 0
    alpha_1 = 0.5
    iter = 0 #iteration number
    elapsedqp = 0
    numX = shape(x_init)[0]
    x0 = zeros(numX)

    if verbose_level:
        print('Solving problem %s \n', prob['name'])
        print('Iteration delta_t   t   Success\n')

    p_0 = p_init
    while t<=1:
        #Calculating the step
        tk = t + delta_t
        p_t = (1-tk)*p_0 +tk*p_final
        step = p_t + p_init
        
        #Updating bound constraints
        if lb_init.any():
            lb = lb_init-x_init
            ub = ub_init-x_init
        elif not lp_init:
            lb = array([])
            ub = array([])

        #Solve QP problem
        y, qp_val, qp_exit, lam_qpopt, mu_qpopt, qptime = qp_solve(prob, p, x_init, y_init, step, lb, ub, N, x0, lb_init, ub_init)
        elapsedqp += qptime
        raw_input()
        if qp_exit == 'infeasible':#QP infeasible
            delta_t = alpha_1*t                                   #shorten step
            t = t - delta_t
            #Print out iteration number and failure
            iter = iter + 1
            success = 0
            if verbose_level:
                print '%f    %f  %f   %d' %(iter, delta_t, t, success)
        else:#QP feasible
            #Update states, multipliers, parameter and time step
            x_init = x_init + y
            y_init['lam_x'] = y_init['lam_x'] + lam_qpopt['lam_x']
            t = t + delta_t
            p_init = p_t
            #Print out iteration number and success
            iter = iter + 1
            success = 1
            if verbose_level:
                print '%f    %f  %f   %d' %(iter, delta_t, t, success)

        if (1-t) <= 1e-5:
            break
                
    return x_init, y_init, elapsedqp
