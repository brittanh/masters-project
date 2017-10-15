'''
Created on 9. nov. 2015

@author: suwartad
implementation of jpredictor.m in python
'''
# -*- coding: utf-8 -*-
import numpy as np
from qp_solve import *
from lp_solve import *

def jpredictor(prob, derivatives, p_init, p_final, x_init, y_init, delta_t, lb, ub, verbose_level):
    """
    please refer to file jpredictor_tr.m for more detailed explaination
    inputs:
    - problem: definition which includes gradient information (Jacobian, Hessian, constraint, etc)
    - p_init: initial value of parameter
    - p_final: final value of parameter
    - x_init: initial primal variable
    - y_init: initial dual variable
    - delta_t: initial increment of parameter change
    - lb and ub: lower and upper bounds for primal variable
    - verbose_level: display with print or not
    
    outputs:
    - x_init: primal variable
    - y_init: dual variable
    - info: information for plotting purposes
    """
    
    # initial info value
    in_info  = 0
    t_list   = np.array([])
    x_list_0 = np.array([])
    x_list_1 = np.array([])
    y_list   = np.array([])
    
    '''
    # error here fix it later...
    t_list[in_info] = 0
    x_list[in_info] = x_init
    y_list[in_info] = y_init
    '''
    
    # initial algorithm parameters
    t = 0.0
    alpha_1  = 0.25
    alpha_2  = 0.5
    iter     = 0
    gamma    = 1.5
    qpv_init = 0
    delta_t  = 0.1
    ok_code  = 0
    t_list   = np.append(t_list,t)
    x_list_0 = np.append(x_list_0,x_init[0])
    x_list_1 = np.append(x_list_1,x_init[1])
        
    while (t < 1):
        
        print "--------------------------------------------------------------------------------------"
        print "iteration number: %d \n" %(iter) 
        
        # calculate step
        step = delta_t * (p_final - p_init)
        # solve QP problem
        qp_exit, qp_val, y, dual, oqp, k_zero_tilde, k_plus_tilde, g = qp_solve(prob, derivatives, p_init, x_init, y_init, step, lb, ub)
        print qp_val
        raw_input()
        
        if (qp_exit == 'optimal'):
            # call LP solver 
            delta_l, lp_exit = lp_solve(oqp, y_init, step, y, k_zero_tilde, k_plus_tilde, g) #? 
            
            if (lp_exit == 'optimal'):
                # update states, multipliers, parameter, and time step
                x_init  = np.array(x_init + y)
                y_init  = delta_l.reshape(-1,1) #transpose
                t       = t + delta_t
                p_init  = p_init + step + t*(p_final - p_init)
                in_info = in_info + 1
                
                if t <= 1:
                    t_list = np.append(t_list,t)
                    x_list_0 = np.append(x_list_0,x_init[0])
                    x_list_1 = np.append(x_list_1,x_init[1])

		print "parameter values" 
		print p_init
		print "\n primal variables"
		print x_init
		print "\n dual variables"
		print y_init

            else:
                # shorten step
                delta_t = alpha_1 * t
                t       = t - alpha_1 * delta_t
            
        else:
            # shorten step
            delta_t = alpha_2*t
            t       = t - alpha_2 * delta_t
            
        iter += 1
        if t==1:
            pass
            
    #return x_init, y_init, info
    return x_init, y_init, t_list, x_list_0, x_list_1
