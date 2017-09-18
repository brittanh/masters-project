'''
Created on 17. nov. 2015

@author: suwartad
'''
# -*- coding: utf-8 -*-
import numpy as np
#import casadi as cas
from openopt import QP

def indices(a, func):
    return [i for (i, val) in enumerate(a) if func(val)]
        

def qp_solve(prob, derivatives, p_init, x_init, y_init, step, lb, ub):
    '''
    QP solver for path-following algorithm.
    input: problem name (prob)
           parameter (p_init), 
           primal variable (x_init),
           dual variable (y_init),
           step,
           lower bound (lb),
           upper bound (up).
    output: solution primal variable (y),
            objective function value (qp_val),
            return status of QP solver (qp_exit),
            derivates of the prob (oqp),
            active-set index (k_zero_tilde),
            inactive-set index (k_plus_tilde),
            g (gradient of objective function
    '''
    # evaluate derivates of prob 
    fval,g,H,Lxp,cin,J,cp,Jeq,dpe = derivatives(x_init, y_init, p_init)
    #f = cas.mul(Lxp,step) + g
    f = np.dot(Lxp,step) + g
    
    n, neq, niq, name = prob()
    if (niq > 0):
        print "cin"
        print cin 
        k_plus_tilde = indices(y_init, lambda x: x > 1e-5)
        k_zero_tilde = indices(y_init, lambda x: x <= 1e-5)
        
        nk_pt = len(k_plus_tilde)
        nk_zt = len(k_zero_tilde)
        
        if (nk_zt > 0):
            if (nk_zt==niq):
                A = J
                #b = -cas.mul(np.matrix(cp),step)-cin
                b = -np.dot(np.matrix(cp),step)-cin
            else:
                A = np.zeros([nk_zt,n])
                b = np.zeros([nk_zt,1])
                for i in np.arange(0,nk_zt):
                    index  = k_zero_tilde[i]
                    A[i,:] = J[index,:]
                    #b[i,:] = -(cas.mul(np.matrix(cp[index,:]),step))-cin[index]
                    b[i,:] = -(np.dot(np.matrix(cp[index,:]),step))-cin[index]
        else:
            A = np.array([])
            b = np.array([])
    
    if (neq > 0):
        Aeq = Jeq
        #beq = -(cas.mul((dpe.T)*step))
        beq = -(np.dot((dpe.T)*step))
    else:
        if (nk_pt > 0):
            Aeq = np.zeros([nk_pt,n])
            beq = np.zeros([nk_pt,1])
            for i in np.arange(0,nk_pt):
                index    = k_plus_tilde[i]
                Aeq[i,:] = J[index,:]
                #beq[i,:] = -(cas.mul(np.matrix(cp[index,:]),step))-cin[index]
                beq[i,:] = -(np.dot(np.matrix(cp[index,:]),step))-cin[index]
        else:
            Aeq = np.array([])
            beq = np.array([])
            
    # storing variables for LP solver. for now skip this process...
    if lb:
        p = QP(H,f,A = A,b = b ,Aeq = Aeq,beq = beq, lb = lb, ub = ub)
    else:
        p = QP(H,f,A = A,b = b ,Aeq = Aeq,beq = beq)
 
    r = p._solve('cvxopt_qp', iprint = 0)
    
    #qp_val = r.ff
    primal = np.array([r.xf])
    primal = primal.T
    
    #qp_exit = r.msg
    dual  = np.array([r.duals])
    dual  = dual.T        

    oqp = {'neq':neq,'niq':niq,'cin':cin,'J':J,'cp':cp,'Jeq':Jeq,'Lxp':Lxp,'H':H}
    #return primal, val, success, k_zero_tilde, k_plus_tilde, g, dual
    return r.msg, r.ff, primal, dual, oqp, k_zero_tilde, k_plus_tilde, g
    
        
    
        