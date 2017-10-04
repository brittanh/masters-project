#!/opt/local/bin/python
# -*- encoding: ascii -*-
"""
@purpose: Path-following algorithm
@author: Brittany Hall
@date: 20.09.2017
@version: 0.1
@updates:
"""
from numpy import *
from nlp_solve import *
from qp_solve import *

def pathfollowing(t, delta_t, alpha1, p_init, p_final, x_init, y_init, chi, lam_opt, mu_opt, case):
    """
    Applying a path following algorithm to an NLP
    """

    #defining empty arrays
    t_list = array([])
    x_list_0 = array([])
    x_list_1 = array([])
    y_list = array([])
    lam_list = array([])
    iter = 0
    #appending initial values
    t_list = append(t_list, t)
    x_list_0 = append(x_list_0,x_init[0])
    x_list_1 = append(x_list_1,x_init[1])
    
    while (t<1):
        print "-----------------------------------------------------"
        print "Iteration number: %d \n" %(iter)

        #calculate step for p
        del_p = delta_t*(p_final-p_init)

        #Solve QP problem
        optimal,x_qpopt,lam_qpopt,mu_qpopt = qp_solve(prob, obj, p_init, x_init, y_init, lam_opt, mu_opt)

        #redefining variables
        del_chi = x_qpopt
        del_lam = lam_qpopt
        del_mu = mu_qpopt
        #determining if a solution was found (IMPROVE)
        if del_chi.size1() == len(x_init):
            qp_exit = 'optimal'
        else:
            qp_exit = ''
        if (qp_exit == 'optimal'):
                chi = chi + del_chi
                if case == 'pure-predictor':
                    #calculate gradient with respect to p
                    lam_list[p_init + del_p] = lam_opt + del_lam* del_p
                    mu_list[p_init + del_p] = mu_opt + del_mu*del_p
                elif case == 'predictor-corrector':
                    lam = lam_init - del_lam
                    mu = mu_init - del_mu
                    lam_list[t] = lam
                    mu_list[t] = mu
                t = t + delta_t
                if t<= 1:
                    t_list = append(t_list, t)
                    x_list_0 = append(x_list_0, x_init[0])
                    x_list_1 = append(x_list_1, x_init[1])
        else:
            delta_t = alpha1*delta_t
            t = t-alpha1*delta_t
        iter += 1

    return x_init, y_init, t_list, x_list_0, x_list_1
