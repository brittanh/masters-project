#!/opt/local/bin/python
# -*- encoding: ascii -*-
"""
@purpose: Path-following algorithm (algorithm 2 from Suwartadi et al 2016)
@author: Brittany Hall
@date: 20.09.2017
@version: 0.1
@updates:
"""
from numpy import array, append, zeros
from nlp_solve import *
from qp_solve import *

def pathfollowing(p_init, p_final, x_init, x_opt, y_init, lam_opt, mu_opt, case):
    """
    Applying a path following algorithm to an NLP
    """

    #defining empty arrays
    t = 0.0
    t_list = array([])
    x_list_0 = array([])
    x_list_1 = array([])
    y_list = array([])
    lam_list = array([])
    mu_list = array([])
    iter = 1
    
    #appending initial values
    t_list = append(t_list, t)
    x_list_0 = append(x_list_0,x_init[0])
    x_list_1 = append(x_list_1,x_init[1])
    lam_list = append(lam_list,lam_opt)
    mu_list = append(mu_list,mu_opt)
    
    #initial algorithm parameters
    
    delta_t = 0.1                                                 #t step size
    N = int(1/delta_t)                                    #number of iterations
    alpha1 = 0.25
    p = zeros((N+1,2))
    p[0,:] = (1-t)*p_init + t*p_final
    for k in range(1,N+1):
        print "-----------------------------------------------------"
        print "Iteration number: %d \n" %(iter)

        #calculate step for p
        p[k,:] = (1-t)*p_init + t*p_final
        step = p[k,:] - p[k-1,:]
        print 't:', t
        print 'Step:', step
        print 'p:', p[k,:]
        if case == 'pure-predictor':
            param = p[k,:]
        elif case == 'predictor-corrector':
            param = p[k,:] + step
        #Solve QP problem
        qp_exit, optimal, x_qpopt, lam_qpopt, mu_qpopt = qp_solve(prob, obj, param, x_opt, y_init, lam_opt, mu_opt, case)
        print 'QP x:', x_qpopt
        #raw_input()
        
        #redefining variables
        del_x= x_qpopt
        del_lam = lam_qpopt
        del_mu = mu_qpopt

        if (qp_exit == 'optimal'):
            x_opt = x_opt + del_x
            if case == 'pure-predictor':
                lam_opt = lam_opt + del_lam * step
                mu_opt = mu_opt + del_mu * step
                lam_list = append(lam_list, lam_opt)
                mu_list = append(mu_list, mu_opt)
            elif case == 'predictor-corrector':
                lam_opt = del_lam
                mu_opt = del_mu
                lam_list = append(lam_list, lam_opt)
                mu_list = append(mu_list, mu_opt)
            t = t + delta_t
            t_list = append(t_list, t)
            x_list_0 = append(x_list_0, x_opt[0])
            x_list_1 = append(x_list_1, x_opt[1])
        else:
            delta_t = alpha1*delta_t
            t = t-alpha1*delta_t
        iter += 1
        
    return x_opt, y_init, t_list, x_list_0, x_list_1, lam_list, mu_list, p
