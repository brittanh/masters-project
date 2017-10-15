#!/opt/local/bin/python
# -*- encoding: ascii -*-
"""
    @purpose: Solving a QP
    @author: Brittany Hall
    @date: 20.09.2017
    @version: 0.1
    @updates:
"""
from numpy import array, append, zeros
from casadi import vertcat, gradient, jacobian, hessian, Function, conic, SX, mtimes
from problem import prob, obj

#p_init = array([1,-4])                                 #initial parameter value
#x_init = array([0.5,0.6])                              #initial primal variable
#y_init = array([1.2])                                    #initial dual variable
#lam_opt = array([0])
#mu_opt = array([0])

#QP solver
def qp_solve(prob, obj, p_init, x_init, y_init, lam_opt, mu_opt):
    """
    QP solver for path-following algorithm
    inputs: prob - problem description
            obj - problem equations
            p_init - initial parameter
            x_init - initial primal variable
            y_init - initial dual variable
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
    x, p, f, f_fun, con, conf, ubx, lbx, ubg, lbg = obj(x_init, y_init, p_init, neq, niq, n, np)

    #Deteriming constraint types
    eq_con_ind = array([])
    iq_con_ind = array([])
    eq_con = array([])
    iq_con = array([])
    tol = 1e-6
    for i in range(0,len(lbg[0])):
        if lbg[0,i] == 0 + tol:
            eq_con = vertcat(eq_con,con[i])
            eq_con_ind = append(eq_con_ind,i)
        elif lbg[0,i] < 0:
            iq_con = vertcat(iq_con,con[i])
            iq_con_ind = append(iq_con_ind,i)
    #Evaluating constraints at current iteration point
    val = conf(x_init,p_init)

    #Determining which inequality constraints are active
    k_plus_tilde = array([])                                 #active constraint
    k_zero_tilde = array([])                               #inactive constraint
    for i in range(0, len(iq_con_ind)):
        if val[i] >= ubg[0,i]:
            k_zero_tilde = append(k_zero_tilde,i)
        elif val[i] < ubg[0,i]:
            k_plus_tilde = append(k_plus_tilde,i)
    nk_pt = len(k_plus_tilde)                     #number of active constraints
    nk_zt = len(k_zero_tilde)                   #number of inactive constraints

    #Calculating Lagrangian
    lam = SX.sym('lam',neq)         #Lagrangian multiplier equality constraints
    mu = SX.sym('mu',niq)         #Lagrangian multiplier inequality constraints
    lag_f = f + mtimes(lam.T,eq_con)+ mtimes(mu.T,iq_con)  #Lagrangian equation

    #Calculating derivatives
    g = gradient(f, x)                        #Derivative of objective function
    g_fun = Function('g_fun',[x,p], [gradient(f, x)])
    H = 2*jacobian(gradient(lag_f,x),x)    #Second derivative of the Lagrangian
    H_fun = Function('H_fun',[x,p,lam,mu],[jacobian(jacobian(lag_f,x),x)])

    if len(eq_con_ind)>0:
        deq = jacobian(eq_con,x)            #Derivative of equality constraints
    else:
        deq = array([])
    if len(iq_con_ind)>0:
        diq = jacobian(iq_con,x)          #Derivative of inequality constraints
    else:
        diq = array([])

    #Creating constraint matrices
    nc = niq + neq                                 #Total number of constraints
    if (niq>0) and (neq>0):                #Equality and inequality constraints
        if (nk_zt >0):                              #Inactive constraints exist
            A = SX.zeros((nc,n))
            A[0,:] = deq                                              #A matrix
            lba = SX.zeros((nc,1))
            lba[0,:] = -eq_con                                #lower bound of A
            uba = SX.zeros((nc,1))
            uba[0,:]= -eq_con                                 #upper bound of A
            for j in range(0,nk_pt): #adding active constraints
                A[neq+j+1,int(k_plus_tilde[j])] = diq[int(k_plus_tilde[j])]
                lba[neq+j+1,int(k_plus_tilde[j])] = -iq_con[int(k_plus_tilde[j])]
                #uba[neq+1,:] = zeros((nk_pt,1))
            for i in range(0,nk_zt): #adding inactive constraints
                A[neq+nk_pt+i+1,int(k_zero_tilde[i])] = diq[int(k_zero_tilde[i])]
                lba[neq+nk_pt+i+1,int(k_zero_tilde[i])] = -iq_con[int(k_zero_tilde[i])]
                uba[neq+nk_pt+i+1,int(k_zero_tilde[i])] = -iq_con[int(k_zero_tilde[i])]
        else:                                          #Active constraints only
            A = vertcat(deq,diq)
            lba = vertcat(-eq_con,-iq_con)
            uba = vertcat(-eq_con,-iq_con)
    elif (niq>0) and (neq==0):                           #Inquality constraints
        if (nk_zt >0):                              #Inactive constraints exist
            A = SX.zeros((nc,n))
            lba = SX.zeros((nc,1))
            uba =  SX.zeros((nc,1))
            for j in range(0,nk_pt): #adding active constraints
                A[j,int(k_plus_tilde[j])] = diq[int(k_plus_tilde[j])]
                lba[j] = -iq_con[int(k_plus_tilde[j])]
                uba[j] = -iq_con[int(k_plus_tilde[j])]
            for i in range(0,nk_zt): #adding inactive constraints
                A[nk_pt+i,int(k_zero_tilde[i])] = diq[int(k_zero_tilde[i])]
                lba[nk_pt+i] = -iq_con[int(k_zero_tilde[i])]
                #uba[nk_pt+i] = -iq_con[int(k_zero_tilde[i])]
        else:
            A = vertcat(deq,diq)
            lba = -iq_con
            uba = -iq_con
    elif (niq==0) and (neq>0):                            #Equality constriants
        A = deq
        lba = -eq_con
        uba = -eq_con
    A_fun = Function('A_fun',[x,p],[A])
    lba_fun = Function('lba_fun',[x,p],[lba])
    uba_fun = Function('uba_fun',[x,p],[uba])

    #Checking that matrices are correct sizes and types
    if (H.size1() != n) or (H.size2() != n) or (H.is_dense()=='False'):
        #H matrix should be a sparse (nxn) and symmetrical
        print('WARNING: H matrix is not the correct dimensions or matrix type')
    if (g.size1() != n) or (g.size2() != 1) or g.is_dense()=='True':
        #g matrix should be a dense (nx1)
        print('WARNING: g matrix is not the correct dimensions or matrix type')
    if (A.size1() !=(neq+niq)) or (A.size2() != n) or (A.is_dense()=='False'):
        #A should be a sparse (nc x n)
        print('WARNING: A matrix is not the correct dimensions or matrix type')
    if lba.size1() !=(neq+niq) or (lba.size2() !=1) or lba.is_dense()=='False':
        print('WARNING: lba matrix is not the correct dimensions or matrix type')
    if uba.size1() !=(neq+niq) or (uba.size2() !=1) or uba.is_dense()=='False':
        print('WARNING: uba matrix is not the correct dimensions or matrix type')

    #Evaluating QP matrices at optimal points
    H_opt = H_fun(x_init,p_init,lam_opt,mu_opt)
    g_opt = g_fun(x_init, p_init)
    A_opt = A_fun(x_init,p_init)
    lba_opt = lba_fun(x_init,p_init)
    uba_opt = uba_fun(x_init,p_init)

    #Defining QP structure
    qp = {}
    qp['h'] = H_opt.sparsity()
    qp['a'] = A_opt.sparsity()
    optimize = conic('optimize','qpoases',qp)
    optimal = optimize(h=H_opt, g=g_opt, a=A_opt, lba=lba_opt, uba=uba_opt, x0=x_init)

    x_qpopt = optimal['x']
    if x_qpopt.shape == x_init.shape:
        qp_exit = 'optimal'
    else:
        qp_exit = ''
    lag_qpopt = optimal['lam_a']
    print lag_qpopt
    #determing Lagrangian constants
    lam_qpopt = zeros((nk_pt,1)) #Lagrange multiplier of active constraints
    mu_qpopt = zeros((nk_zt,1)) #Lagrange multiplier of inactive constraints
    if nk_pt > 0:
        for j in range(0,len(k_plus_tilde)):
            lam_qpopt[j] = lag_qpopt[int(k_plus_tilde[j])]
    print nk_zt

    if nk_zt > 0:
        for k in range(0,len(k_zero_tilde)):
            print lag_qpopt[int(k_zero_tilde[k])]
            mu_qpopt[k] = lag_qpopt[int(k_zero_tilde[k])]
    raw_input()
    return qp_exit, optimal, x_qpopt, lam_qpopt, mu_qpopt
