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
def qp_solve(prob, obj, p_init, x_init, y_init, lam_opt, mu_opt,case):
    """
    QP solver for path-following algorithm
    inputs: prob - problem description
            obj - problem equations
            p_init - initial parameter
            x_init - initial primal variable
            y_init - initial dual variable
            lam_opt - Lagrange multipliers of equality and active constraints
            mu_opt - Lagrange multipliers of inequality constraints
    outputs: y - solution primal variable
            qp_val - objective function value
            qp_exit - return status of QP solver
            deriv - derivatives of the problem
            k_zero_tilde - active set index
            k_plus_tilde - inactive set index
            grad - gradient of objective function
    """
    print 'Current point x:', x_init
    #Importing problem to be solved
    nx, np, neq, niq, name = prob()
    x, p, f, f_fun, con, conf, ubx, lbx, ubg, lbg = obj(x_init, y_init, p_init, neq, niq, nx, np)

    #Deteriming constraint types
    eq_con_ind = array([]) #indices of equality constraints
    iq_con_ind = array([]) #indices of inequality constraints
    eq_con = array([]) #equality constraints
    iq_con = array([]) #inequality constraints
    
    for i in range(0,len(lbg[0])):
        if lbg[0,i] == 0:
            eq_con = vertcat(eq_con,con[i])
            eq_con_ind = append(eq_con_ind,i)
        elif lbg[0,i] < 0:
            iq_con = vertcat(iq_con,con[i])
            iq_con_ind = append(iq_con_ind,i)
#    print 'Equality Constraint:', eq_con
#    print 'Inequality Constraint:', iq_con

#    if case == 'pure-predictor':
    
#       return qp_exit, optimal, x_qpopt, lam_qpopt, mu_qpopt
    if case == 'predictor-corrector':
        #Evaluating constraints at current iteration point
        con_vals = conf(x_init,p_init)
        #Determining which inequality constraints are active
        k_plus_tilde = array([])                             #active constraint
        k_zero_tilde = array([])                           #inactive constraint
        tol = 10e-5 #tolerance
        for i in range(0, len(iq_con_ind)):
            if ubg[0,i] - tol <= con_vals[i] and  con_vals[i] <= ubg[0,i]+tol:
                k_plus_tilde = append(k_plus_tilde,i)
            else:
                k_zero_tilde = append(k_zero_tilde,i)
#        print 'Active constraints:', k_plus_tilde
#        print 'Inactive constraints:', k_zero_tilde
#        print 'Constraint values:', con_vals

        nk_pt = len(k_plus_tilde)                 #number of active constraints
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
            #this part needs to be tested
            if (nk_zt >0):                          #Inactive constraints exist
                A = SX.zeros((nc,nx))
                print deq
                A[0,:] = deq                                          #A matrix
                lba = -1e16*SX.zeros((nc,1))
                lba[0,:] = -eq_con                            #lower bound of A
                uba = 1e16*SX.zeros((nc,1))
                uba[0,:]= -eq_con                             #upper bound of A
                for j in range(0,nk_pt): #adding active constraints
                    A[neq+j+1,:] = diq[int(k_plus_tilde[j]),:]
                    lba[neq+j+1] = -iq_con[int(k_plus_tilde[j])]
                    uba[neq+j+1] = -iq_con[int(k_plus_tilde[i])]
                for i in range(0,nk_zt): #adding inactive constraints
                    A[neq+nk_pt+i+1,:] = diq[int(k_zero_tilde[i]),:]
                    uba[neq+nk_pt+i+1] = -iq_con[int(k_zero_tilde[i])]
                    #inactive constraints don't have lower bounds
            else:                                      #Active constraints only
                A = vertcat(deq,diq)
                lba = vertcat(-eq_con,-iq_con)
                uba = vertcat(-eq_con,-iq_con)
        elif (niq>0) and (neq==0):                           #Inquality constraints
            if (nk_zt >0):                              #Inactive constraints exist
                A = SX.zeros((nc,nx))
                lba = -1e16*SX.ones((nc,1))
                uba =  1e16*SX.ones((nc,1))
                for j in range(0,nk_pt): #adding active constraints
                    A[j,:] = diq[int(k_plus_tilde[j]),:]
                    lba[j] = -iq_con[int(k_plus_tilde[j])]
                    uba[j] = -iq_con[int(k_plus_tilde[j])]
                for i in range(0,nk_zt): #adding inactive constraints
                    A[nk_pt+i,:] = diq[int(k_zero_tilde[i]),:]
                    uba[nk_pt+i] = -iq_con[int(k_zero_tilde[i])]
                    #inactive constraints don't have lower bounds
            else:
                raw_input()
                A = vertcat(deq,diq)
                lba = -iq_con
                uba = -iq_con
        elif (niq==0) and (neq>0):                        #Equality constriants
            A = deq
            lba = -eq_con
            uba = -eq_con
        A_fun = Function('A_fun',[x,p],[A])
        lba_fun = Function('lba_fun',[x,p],[lba])
        uba_fun = Function('uba_fun',[x,p],[uba])
        #Checking that matrices are correct sizes and types
        if (H.size1() != nx) or (H.size2() != nx) or (H.is_dense()=='False'):
            #H matrix should be a sparse (nxn) and symmetrical
            print('WARNING: H matrix is not the correct dimensions or matrix type')
        if (g.size1() != nx) or (g.size2() != 1) or g.is_dense()=='True':
            #g matrix should be a dense (nx1)
            print('WARNING: g matrix is not the correct dimensions or matrix type')
        if (A.size1() !=(neq+niq)) or (A.size2() != nx) or (A.is_dense()=='False'):
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

#        print 'Lower bounds', lba_opt
#        print 'Upper bounds', uba_opt
#        print 'Bound matrix', A_opt

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

        #Determing Lagrangian multipliers (lambda and mu)
        lam_qpopt = zeros((nk_pt,1))     #Lagrange multiplier of active constraints
        mu_qpopt = zeros((nk_zt,1))    #Lagrange multiplier of inactive constraints
        if nk_pt > 0:
            for j in range(0,len(k_plus_tilde)):
                lam_qpopt[j] = lag_qpopt[int(k_plus_tilde[j])]
        if nk_zt > 0:
            for k in range(0,len(k_zero_tilde)):
                print lag_qpopt[int(k_zero_tilde[k])]
        return qp_exit, optimal, x_qpopt, lam_qpopt, mu_qpopt
