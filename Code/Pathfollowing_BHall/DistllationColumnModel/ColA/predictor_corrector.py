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
from numpy import zeros

def predictor_corrector(problem, p_init, p_final, x_init, y_init, delta_t, lb_init, ub_init, verbose_level, N):
    
    p = p_init
    pp = SX.sym('pp')
    prob = problem(pp)
    t = 0
    alpha_2 = 0.5
    iter = 0 #number of iterations
    elapsedqp = 0
    numX = shape(x_init,0)
    x0 = zeros((numX,1))

    if verbose_level:
        print('Solving problem %s \n', prob['name'])
        print('Iteration delta_t   t   Success\n')

    p_0 = p_init
    while t<=1:
        #Calculating the step
        tk = t + delta_t
        print p_0
        print p_0.shape
        raw_input()
        p_t = (1-tk)*p_0 +tk*p_final
        step = p_t + p_init
        print step.shape
        raw_input()
        
        #Updating bound constraints
        if lb_init:
            lb = lb_init-x_init
            ub = ub_init-x_init
        elif not lp_init:
            lb = array([])
            ub = array([])
        
        #Solve QP problem
        qp_exit, y, qp_val, lam_qpopt, mu_qpopt, qprun = qp_solve(prob, p, x_init, y_init, step, lb, ub, N, x0, lb_init, ub_init)
        elapsedqp += qprun

        if qp_exit != 'optimal':
            #QP infeasible
            delta_t = alpha_2*t                                   #shorten step
            
            #Print out iteration number and failure
            iter = iter + 1
            success = 0
            if verbose_level:
                print '%f    %f  %f   %d' %(iter, delta_t, t, success)
        else:
            #QP is feasible
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

        if (1-t) <= 1e-5
            break
                
    return x_init, y_init, elapsedqp
