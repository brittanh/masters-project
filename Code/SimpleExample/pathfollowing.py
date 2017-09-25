##!/opt/local/bin/python
## -*- encoding: ascii -*-
#"""
#@purpose: Path-following algorithm
#@author: Brittany Hall
#@date: 20.09.2017
#@version: 1
#@updates:
#"""
#from numpy import *
#from nlp_solve import *
#
#def pathfollowing(t, delta_t, alpha1, p_init, p_final, x_init, y_int, chi, lam_opt, mu_ops):
#    """
#    Applying a path following algorithm to an NLP
#    """
#
#    #defining empty arrays
#    t_list = array([])
#    x_list_0 = array([])
#    x_list_1 = array([])
#    y_list = array([])
#    lam_list = array([])
#    iter = 0
#    #appending initial values
#    t_list = append(t_list, t)
#    x_list_0 = append(x_list_0,x_init[0])
#    x_list_1 = append(x_list_1,x_init[1])
#
#    raw_input()
#    
#    while (t<1):
#        print "-----------------------------------------------------"
#        print "Iteration number: %d \n" %(iter)
#
#        #calculate step for p
#        del_p = delta_t*(p_final-p_init)
#
#        #Solve QP problem
#        del_chi, lam , mu = qp_solve(prob, obj, p_init, x_init, y_init, del_p)
#        
#        if (qp_exit == 'optimal'):
#                chi = chi + del_chi
#                if case == 'pure-predictor':
#                    #calculate gradient with respect to p
#                    lam_list(p_init + del_p) = lam_opt + grad_p* del_p
#                    mu_list(p_init + del_p) = mu_opt + grad_p
#                elif case == 'predictor-corrector':
#                    del_lam = lam_init - lam
#                    del_mu = mu_init - mu
#                    lam_list(t) = del_lam
#                    mu_list(t) = del_mu
#                t = t + delta_t
#                k = k+1
#                if t<= 1:
#                    t_list = append(t_list, t)
#                    x_list_0 = append(x_list_0, x_init[0])
#                    x_list_1 = append(x_list_1, x_init[1])
#        else:
#            delta_t = alpha1*delta_t
#            t = t-alpha1*delta_t
#        iter += 1
#
#    return x_init, y_init, t_list, x_list_0, x_list_1
