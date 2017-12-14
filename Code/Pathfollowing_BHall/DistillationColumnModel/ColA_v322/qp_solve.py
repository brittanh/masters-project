#!/opt/local/bin/python
# -*- encoding: ascii -*-
"""
    @purpose: Solving a QP
    @author: Brittany Hall
    @date: 08.10.2017
    @version: 0.1
    @updates:
"""
from numpy import where, multiply
from casadi import *
from params import params
from objective import objective
import time

def qp_solve(prob, p_init, x_init, y_init, step, lb, ub, N, x0, lb_init, ub_init):
    """
    QP solver for path-following algorithm
    inputs: prob - problem description
            p - parameters
            x_init - initial primal variable
            y_init - initial dual variable
            step - step to be taken (in p)
            lb_init - lower bounds
            ub_init - upper bounds
            verbose_level - amount of output text
            N - iteration number
    outputs: y - solution primal variable
            qp_val - objective function value
            qp_exit - return status of QP solver
                
    """
    
    #Importing problem to be solved

    neq = prob['neq']                           #Number of equality constraints
    niq = prob['niq']                         #Number of inequality constraints
    name = prob['name']                                        #Name of problem
    _, g, H, Lxp, cst, _, _, Jeq, dpe, _ = objective(x_init,y_init,p_init,N,params) #objective function

    #Setting up QP
    f = mtimes(Lxp,step) + g
    
    #Constraints
    ceq = cst
    Aeq = Jeq
    beq = mtimes(dpe,step) + ceq

    #Check Lagrange multipliers from bound constraints
    lamC = fabs(y_init['lam_x'])
    BAC = where(lamC >= 1e-3) #setting limits to determine if constraint is active

    #Finding active constraints
    numBAC = len(BAC)
    for i in range(0,numBAC):
        #Placing strongly active constraint on boundary
        indB = BAC[i]
        #Keeping upper bound on boundary
        ub[indB] = 0
        lb[indB] = 0

    #Solving the QP
    qp = {}
    qp['h'] = H.sparsity()
    qp['a'] = Aeq.sparsity()
    optimize = conic('optimize','qpoases',qp,{'sparse':True})
    #optimize = conic('optimize','gurobi',qp)
    startqp = time.time()
    optimal = optimize(h=H, g=f, a=Aeq, lba=beq, uba=beq, lbx=lb, ubx=ub)
    elapsedqp = time.time()-startqp
    x_qpopt = optimal['x']
    lam_qpopt = optimal['lam']
    mu_qpopt = optimal['mu']
    
    if x_qpopt.shape == x_init.shape:
        qp_exit = 'optimal'
    else:
        qp_exit = ''

    return qp_exit, y, x_qpopt, lam_qpopt, mu_qpopt, qprun
