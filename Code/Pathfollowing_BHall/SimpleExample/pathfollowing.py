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

def pathfollowing(p_init, p_final, x_init, y_init, chi, lam_opt, mu_opt, case):
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
    iter = 0
    
    #appending initial values
    t_list = append(t_list, t)
    x_list_0 = append(x_list_0,x_init[0])
    x_list_1 = append(x_list_1,x_init[1])
    
    #initial algorithm parameters
    
    delta_t = 0.1                                                  #t step size
    N = int(1/delta_t)                                    #number of iterations
    alpha1 = 0.25
    p = zeros((N,2))
    p[0,:] = (1-t)*p_init + t*p_final
    
    for k in range(1,N):
        print "-----------------------------------------------------"
        print "Iteration number: %d \n" %(iter)

        #calculate step for p
        p[k] = (1-t)*p_init + t*p_final
        step = p[k] - p[k-1]
        
        #Solve QP problem
        qp_exit, optimal, x_qpopt, lam_qpopt, mu_qpopt = qp_solve(prob, obj, p_init, x_init, y_init, lam_opt, mu_opt)

        #redefining variables
        del_chi = x_qpopt
        del_lam = lam_qpopt
        del_mu = mu_qpopt

        if (qp_exit == 'optimal'):
            chi = chi + del_chi
            if case == 'pure-predictor':
                lam = lam_opt + del_lam * step
                mu = mu_opt + del_mu * step
                lam_list = append(lam_list, lam)
                mu_list[k] = append(mu_list, mu)
            elif case == 'predictor-corrector':
#                    lam = lam_init - del_lam
#                    mu = mu_init - del_mu
                lam = del_lam
                mu = del_mu
                lam_list = append(lam_list, lam)
                mu_list = append(mu_list, mu)
            t = t + delta_t
            t_list = append(t_list, t)
            x_list_0 = append(x_list_0, chi[0])
            x_list_1 = append(x_list_1, chi[1])
            #Updating values to send to the QP
            p_init = p[k]
            x_init = chi
            lam_opt = lam
            mu_opt = mu
        else:
            delta_t = alpha1*delta_t
            t = t-alpha1*delta_t
        k += 1
        
    return x_init, y_init, t_list, x_list_0, x_list_1, lam_list, mu_list, p
