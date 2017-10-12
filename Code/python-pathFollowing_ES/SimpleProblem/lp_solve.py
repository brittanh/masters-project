'''
Created on Nov 20, 2015

@author: detu
'''
# -*- coding: utf-8 -*-
import numpy as np
#import casadi as cas
from openopt import LP
from qp_solve import indices

def lp_solve(oqp, y_init, step, y, k_zero_tilde, k_plus_tilde, g):
    '''
    LP solver for pathfollowing algorithm to get dual variables
    '''
    ny      = y_init.shape[0]
    delta_l = np.zeros([ny,1])
    nx      = y.shape[0]
    #lb_delta_lamda = np.array([])
    #lb_delta_tau   = np.array([])
    
    '''
    CONSTRUCT CONSTRAINT
    '''
    if (len(k_zero_tilde) >0):
        #nk_zt      = size(k_zero_tilde,1)
        nk_zt      = len(k_zero_tilde)
        #k_zero_hat = zeros(nk_zt,1)
        k_zero_hat = np.zeros([nk_zt,1])
        
        #for i=1:nk_zt
        for i in np.arange(0,nk_zt):
            index = k_zero_tilde[i]
            #k_zero_hat[i] = oqp.cin[index] + oqp.J[index,:]*y + oqp.cp[index,:]*step
            #k_zero_hat[i] = oqp.cin[index] + np.dot(np.matrix(oqp.J[index,:]),y) + np.dot(np.matrix(oqp.cp[index,:]),step)
            k_zero_hat[i] = oqp['cin'][index] + np.dot(np.matrix(oqp['J'][index,:]),y) + np.dot(np.matrix(oqp['cp'][index,:]),step)
        
        #index_zh   = find(k_zero_hat < 1e-5)
        index_zh = indices(k_zero_hat, lambda x: x > 1e-5)
        nk_zh    = len(index_zh)
        
        #remove_index = zeros(nk_zh,1)
        remove_index = np.zeros([nk_zh,1])
        if (nk_zt == nk_zh):
            for i in np.arange(0,nk_zt):
                remove_index[i]          = k_zero_tilde[i]
                delta_l[k_zero_tilde[i]] = 0
    
    # equality constraint
    # construct matrix Aeq = [a11 a12; a21 a22]
    if (oqp['neq']):
        #a11 = oqp.Jeq'
        a11 = (oqp['Jeq']).T
        if (oqp['niq']):
            a12 = np.zeros([a11.shape])
        else:
            a12 = np.array([])
    else:
        a11            = np.array([])
        a12            = np.array([])
        #delta_lamda    = np.array([])
        #lb_delta_lamda = np.array([])
        
    if (oqp['niq']):
        #a22 = oqp.J'
        a22 = (oqp['J']).T
        #beq = -oqp.Lxp*step - oqp.H*y - g
        beq = np.dot(-np.matrix(oqp['Lxp']),step) - np.dot(np.matrix(oqp['H']),y) - g
        
        # add additional constraint delta_eta = 0
        if (ny > nx):
            if (nk_zh > 1):
                #adc = zeros(nk_zh,size(a22,2))
                adc = np.zeros(nk_zh,a22.shape[1])
                for i in np.arange(0,nk_zh):
                    adc[i,k_zero_tilde[i]] = 1
                    beq = np.vstack(np.array([beq]),0)
                #a22 = [a22;adc]
                a22 = np.vstack((a22,adc))
            else:
                #adc = zeros(1,size(a22,1))
                adc = np.zeros([1,a22.shape[0]])
                adc[k_zero_tilde] = 1
                #a22 = [a22;adc]
                a22 = np.vstack((a22,adc))
                #beq = [beq; 0]
                beq = np.vstack(np.array([beq]),0)
        
        # matrix arrangement
        if (oqp['neq']):
            #a21 = zeros(size(a22))
            a21 = np.zeros([a22.shape])
        else:
            a21 = np.array([])
    
    # setup Aeq and beq
    #Aeq = [a11 a12; a21 a22]
    #A11  = np.hstack((a11,a12))
    #A12  = np.hstack((a21,a22)) # error here ! dimension must be the same... continue in Bergen
    #Aeq  = np.vstack((A11,A12))
    # temporary solution !
    if a11.shape[0] > 0:
        Aeq = a11
    if a12.shape[0] > 0:
        Aeq = a12
    if a21.shape[0] > 0:
        Aeq = a21
    if a22.shape[0] > 0:
        Aeq = a22
    
    # bound constraint
    #ub = np.inf*np.ones([ny, 1])
    #lb = -np.inf*np.ones([ny, 1])
    ub = 1e6*np.ones([1,ny])
    lb = -1e6*np.ones([1,ny])
    
    # Inequality constraint
    A = np.array([])
    b = np.array([])
    
    '''
    CONSTRUCT OBJECTIVE FUNCTION
    consist of inequality and equality parts
    '''
    f = np.array([0.0])   # watch out might be error due to dimension !
    
    if (oqp['niq']):
        #in_p = step'*oqp.cp'
        #in_p = np.dot((step.T),(np.matrix(oqp['cp'])).T)
        in_p = np.dot((step.T),(oqp['cp']).T)
        f    = f + in_p
        #f    = np.array(f);
    
    if (oqp['neq']):
        #eq_p = step'*oqp.dpe'
        eq_p = np.array(np.dot((step.T),(np.matrix(oqp['dpe'])).T))
        f    = f + eq_p
        #f    = np.array(f)
    '''
    Solve LP problem; Maximization problem !
    '''
    #[delta_lp, ~, lp_exit] = linprog(-f,A,b,Aeq,beq,lb,ub,[],optionslin)
    beq = np.array(beq.T)  # everything must be array ! 
    #p = LP(f, A=A, Aeq=Aeq, b=b, beq=beq, lb=lb, ub=ub)
    f = f.reshape(ny,)
    p = LP(f, A=A, Aeq=Aeq, b=b, beq=beq, lb=lb, ub=ub)
    #r = p.maximize('pclp')  # check the solver ! otherwise use p.minimize('glpk') for CVXOPT
    r = p.maximize('glpk')
    
    lp_exit = r.msg
    delta_l = r.xf
    #if len(delta_l >0):
    #    lp_exit = 'optimal'
    
    '''
    DO IT LATER !
    # arrange delta_l based on index
    if (len(k_zero_tilde) > 0):
        if (lp_exit == 'optimal'):
            ndel = 0
            for i in np.arange(0,ny):
            #for i=1:ny:
                if (any(i == remove_index)):
                    ndel = ndel + 1
                    continue
                else:
                    ni = i - ndel
                    delta_l[i] = delta_lp[ni]
    else:
        delta_l = delta_lp
    '''
    return delta_l, lp_exit 
