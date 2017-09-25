#!/opt/local/bin/python
# -*- encoding: ascii -*-
"""
    @purpose: Solving a QP
    @author: Brittany Hall
    @date: 20.09.2017
    @version: 0.1
    @updates:
"""
from numpy import *
from casadi import *
from problem import prob, obj

p_init = array([1,-4])                                 #initial parameter value
x_init = array([0.5,0.6])                              #initial primal variable
y_init = array([1.2])                                    #initial dual baraible

#defining indices
def indices(a, func):
    return [i for (i, val) in enumerate(a) if func(val)]

#QP solver
def qp_solve(prob, obj, p_init, x_init, y_init, lam_opt, mu_opt):
    """
    QP solver for path-following algorithm
    inputs: prob - problem description
            obj - problem equations
            p_init - initial parameter
            x_init - initial primal variable
            y_init - initial dual variable
            del_p - step
    outputs: y - solution primal variable
            qp_val - objective function value
            qp_exit - return status of QP solver
            deriv - derivatives of the problem
            k_zero_tilde - active set index
            k_plus_tilde - inactive set index
            grad - gradient of objective function
    """
    #Importing problem to be solved
    n, np, neq, niq, name = prob()
    x, p, f, f_fun, con, conf, ubx, lbx, ubg, lbg = obj(x_init,y_init,p_init, neq, niq, n, np)

    #Calculating Lagrangian
    #Deteriming constraint types
    eq_con_ind = array([])
    iq_con_ind = array([])
    eq_con = array([])
    iq_con = array([])
    for i in range(0,len(lbg[0])):
        if lbg[0,i] == 0:
            eq_con = vertcat(eq_con,con[i])
            eq_con_ind = append(eq_con_ind,i)
        elif lbg[0,i] > 0:
            iq_con = vertcat(iq_con,con[i])
            iq_con_ind = append(iq_con_ind,i)

    #Evaluating constraints at current iteration point
    val = conf(x_init,p_init)

    #Determining which inequality constraints are active
    k_plus_tilde = array([])                                 #active constraint
    k_zero_tilde = array([])                               #inactive constraint
    for i in range(0, len(iq_con_ind)):
        if val[i] >= ubg[0,i]:                               #active constraint
            k_plus_tilde = append(k_plus_tilde,i)
        elif val[i] < ubg[0,i]:                            #inactive constraint
            k_zero_tilde = append(k_zero_tilde,i)
    nk_pt = len(k_plus_tilde)
    nk_zt = len(k_zero_tilde)

    #Calculating Lagrangian
    lag = f_fun(x_init,p_init) - lam_opt.T*eq_con-mu_opt.T*iq_con

    #Calculating derivatives
    dx = jacobian(f, x)                       #derivative of objective function
    ddx = hessian(lag,x)                   #second derivative of the Lagrangian
    deq = jacobian(eq_con,x)                #derivative of equality constraints
    diq = jacobian(iq_con,x)              #derivative of inequality constraints

    if (niq>0):
        #Determining active constraints

#        nk_pt = len(k_plus_tilde)
#        nk_zt = len(k_zero_tilde)
#        
#        if (nk_zt >0):
#            if (nk_zt == niq):
#                A = array([deq,diq,diq])
#    
#                lba = array([-con_eq,-con_iq, 0])
#                uba = array([-con_eq, -con_iq, -con_iq])
#            else:
#                A = zeros([nk_zt, n])
#                lba = array()
#        else:
#            A = array([])
#            lba = array([])
#            uba = array([])
#    if (neq>0):
#        Aeq
#    #Constructing matrices
#    H = 2*DMatrix(ddx)                                                #H matrix
#    g = DMatrix(g)                                                    #g matrix
#    lbx = DMatrix(lbx)
#    ubx = DMatrix(ubx)
#    lba = DMatrix(lbg)
#    ubg = DMatrix(ubg)
#    A = DMatrix(A)
#    A = reshape((A.size1(),H.size1()))            #making sure dimensions match
#
#    #QP structure
#    qp = qpStruct(h=H.sparsity(),a=A.sparsity())
#    solver = QpSolver("qpoases",qp)
#        
#    #Initialize the solver
#    solver.init()
#        
#    #Pass problem data
#    solver.setInput(H,"h")
#    solver.setInput(g,"g")
#    solver.setInput(A,"a")
#    solver.setInput(lb,"lb")
#    solver.setInput(ub,"ub")
#        
#    #Solve the QP
#    solver.evaluate()
#        
#    return array(solver.getOutput("x"))

qp_solve(prob, obj, p_init, x_init, y_init)
