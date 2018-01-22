#!/opt/local/bin/python
# -*- encoding: ascii -*-
"""
    @purpose: Solving a QP
    @author: Brittany Hall
    @date: 08.10.2017
    @version: 0.1
    @updates:
"""
from casadi import *
from numpy import where, multiply, shape, all, isnan, array, vstack, hstack, ones, reshape, zeros
from params import params
from objective import objective
import time
import osqp
import scipy.sparse as sparse

def qp_solve(prob, p_init, x_init, y_init, step, lb, ub, N, x0, lb_init, ub_init):
    """
    QP solver for path-following algorithm
    inputs: prob - problem description
            p_init - initial parameters
            x_init - initial primal variable
            y_init - initial dual variable
            step - step to be taken (in p)
            lb -
            ub -
            N - iteration number
            x0 - initial guess for primal variable
            lb_init - lower bounds
            ub_init - upper bounds
    outputs: y - solution primal variable
            qp_val - objective function value
            qp_exit - return status of QP solver
                
    """
    
    ##=================Importing problem to be solved========================##

    neq = prob['neq']                           #Number of equality constraints
    niq = prob['niq']                         #Number of inequality constraints
    name = prob['name']                                        #Name of problem
    
    _, g, H, Lxp, cst, _, _, Jeq, dpe, _ = objective(x_init, y_init, p_init, N, params)

    #Setting up QP
    f = mtimes(Lxp,step) + g

    #Constraints
    ceq = cst
    Aeq = Jeq
    beq = mtimes(dpe,step) + ceq

    #Check Lagrange multipliers from bound constraints
    lamC = fabs(y_init['lam_x'])
    
    #Setting limits to determine if constraint is active
    BAC = where(lamC >= 1e-3)
    BAC = BAC[0]

    #Finding active constraints
    numBAC = len(BAC)
    for i in range(0,numBAC):
        #Placing strongly active constraint on boundary
        indB = BAC[i]
        #Keeping upper bound on boundary
        ub[indB] = 0
        lb[indB] = 0

    ##===============Solving the QP using OSQP===============================##
    #Setting up proper format for solver
    nx = params['prob']['nx']
    nu = params['prob']['nu']
    Ax = sparse.eye(shape(Aeq)[1])
    leq = reshape(beq, shape(beq)[0])
    lb = reshape(lb, shape(lb)[0])
    ueq = reshape(beq, shape(beq)[0])
    ub = reshape(ub, shape(ub)[0])
    l = hstack([leq, lb])
    u = hstack([ueq, ub])
    A = sparse.vstack([Aeq, Ax])
    A = sparse.csc_matrix(A)
    P = sparse.csc_matrix(H)
    q = array(f)

    #Starting solver
    prob = osqp.OSQP()
    #ftol = 1e-7
    #xtol = 1e-7
    prob.setup(P, q, A, l, u, eps_rel = 1e-7)
    results = prob.solve()
    if results.info.status != 'solved':
        print "OSQP did not solve the problem!\n"
        qp_exit = 'infeasible'
    else:
        print "OSQP solver runtime: %f\n" %results.info.solve_time
        qp_exit = 'feasible'
    #Results
    y = results.x
    lam_qpopt = results.y[0:shape(Aeq)[0]]
    mu_qpopt = results.y[shape(Aeq)[0]:]
    qp_val = results.info.obj_val
    elapsedqp = results.info.solve_time


    ##===============Solving QP using QPOASES/GUROBI========================##
#    qp = {}
#    qp['h'] = H.sparsity()
#    qp['a'] = Aeq.sparsity()
#    #optimize = conic('optimize','qpoases',qp,{'sparse':True})
#    startqp = time.time()
#    optimal = optimize(h=H, g=f, a=Aeq, lba=beq, uba=beq, lbx=lb, ubx=ub, x0=x0)
#    elapsedqp = time.time()-startqp
#    x_qpopt = optimal['x']                                     #primal solution
#    y = x_qpopt
#    qp_val = optimal['cost']                                      #optimal cost
#    lam_qpopt = optimal['lam_a']                   #dual solution-linear bounds
#    mu_qpopt = optimal['lam_x']                    #dual solution-simple bounds
#
#    if isnan(array(x_qpopt[0])):
#        qp_exit = 'infeasible'
#    else:
#        qp_exit = 'optimal'
#    print qp_val
#    print y
#    raw_input()
    ##====================Solving QP using IPOPT============================##
#    w = MX()
#    nx = params['prob']['nx']
#    nu = params['prob']['nu']          #Number of states
#    nk = params['prob']['nk']
#    X0 = MX.sym('X0', nx)
#    w =  vertcat(w,X0)
#    _, C, D, d = collocationSetup()
#    for iter in range(0,N):
#        for k in range(0,nk):
#            Xkj = {}
#            for j in range(0,d):
#                Xkj[str(j)] = MX.sym('X_' + str((iter)*nk + k)+'_'+str(j+1), nx)
#            #New NLP variable for control
#            Uk = MX.sym('U_'+str((iter)*nk+k),nu)
#            w = vertcat(w,Uk)
#            for j in range(0,d):
#                w = vertcat(w, Xkj[str(j)])
#                Xk = MX.sym('X_'+ str((iter)*nk + k), nx)
#            w = vertcat(w, Xk)
#
#    F = 0.5*mtimes(w.T,mtimes(H,w)) + mtimes(f.T,w)
#    startqp = time.time()
#    qp = {'x':w , 'f':F , 'g': mtimes(Aeq,w)}
#    optimize = nlpsol('optimize', 'ipopt', qp)
#    sol = optimize(lbx=lb, ubx=ub, lbg = beq, ubg=beq)
#    elapsedqp = time.time()-startqp
#    y = sol['x']
#    qp_val = sol['f']
#    lam_x = sol['lam_x']
#    lam_g = sol['lam_g']


    return y, qp_val, qp_exit, lam_qpopt, mu_qpopt, elapsedqp
