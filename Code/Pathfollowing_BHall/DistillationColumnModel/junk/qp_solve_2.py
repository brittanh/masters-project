#!/opt/local/bin/python
# -*- encoding: ascii -*-
"""
    @purpose: Solving a QP using CVXOPT
    @author: Brittany Hall
    @date: 08.10.2017
    @version: 0.1
    @updates:
"""

from casadi import mtimes, fabs
from params import params
from objective import objective
import time
import sys
from cvxopt import matrix as cvxmat, sparse, spmatrix
from cvxopt.solvers import qp, options
from numpy import where, multiply, array, zeros, ones

def qp_solve_2(prob, p_init, x_init, y_init, step, lb, ub, N, x0, lb_init, ub_init):
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
    """
        minimize:
        (1/2)*x'*H*x + f'*x
        subject to:
        Aeq*x = beq
        lb <= x <= ub
    """
    #converting to cvxopt-style matrices
    def _convert(H, f, Aeq, beq, lb, ub, x0):
        """ Convert everything to cvxopt-style matrices """
        P = cvxmat(H)
        print H[0]
        print P[0]
        raw_input()
        q = cvxmat(f)
        A = cvxmat(Aeq)
        b = cvxmat(beq)
        n = lb.size
        G = sparse([-speye(n),speye(n)])
        h = cvxmat(np.vstack[-lb,ub])
        x0 = cvxmat(x0)
        return P, q, G, h, A, b, x0

    def speye(n):
        """Create a sparse identity matrix"""
        r = range(n)
        return spmatrix(1.0, r, r)

    P, q, G, h, A, b, x0 = _convert(H, f, Aeq, beq, lb, ub, x0)
    print P[0]
    raw_input()
    startqp = time.time()
    results = qp(P, q, G, h, A, b, x0)
    print results
    raw_input()
    
    elapsedqp = time.time()-startqp
    x_qpopt = optimal['x']
    if x_qpopt.shape == x_init.shape:
        qp_exit = 'optimal'
    else:
        qp_exit = ''


