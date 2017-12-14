#!/opt/local/bin/python
# -*- encoding: ascii -*-
"""
    @purpose: Solving optimal control problem
    @author: Brittany Hall
    @date: 10.11.2017
    @version: 0.1
    @updates:
    """
from casadi import *
from numpy import ones, zeros, multiply, append
import scipy.io as spio

def itPredHorizon_pf(Xk, V, cons, obj, params, iter, ssoftc):

    #Extracting parameters
    NT = params['dist']['NT']
    sf = params['model']['sf']
    xdot_val_rf_ss = params['model']['xdot_val_rf_ss']
    u_opt = params['model']['u_opt']

    pf = params['price']['pf']
    pV = params['price']['pV']
    pB = params['price']['pB']
    pD = params['price']['pD']
    F_0 = params['dist']['F_0']

    C = params['colloc']['C']
    D = params['colloc']['D']
    h = params['colloc']['h']

    delta_t = params['weight']['delta_time']
    Qmax = params['Qmax']

    nx = params['prob']['nx']
    nu = params['prob']['nu']
    nk = params['prob']['nk']
    d = params['prob']['d']
    ns = params['prob']['ns']

    count = 0
    for k in range(0,nk):
        #New NLP variable for control
        Uk = MX.sym('U_'+str((iter)*nk+k), nu)
        V = vertcat(V,Uk)
        Jcontrol = mtimes(transpose(multiply(Qmax[nx:nx+nu], Uk - u_opt)), (Uk - u_opt))

        #State at collocation points
        SumX1 = 0
        Xkj = {}
        for j in range(0,d):
            Xkj[str(j)] = MX.sym('X_' + str((iter)*nk + k) +'_'+str(j+1), nx)
            V = vertcat(V,Xkj[str(j)])
            count += 1

        #Loop over collocation points
        Xk_end = D[0] * Xk
        for j in range(0,d):
            xp = C[0,j+1] * Xk
            for r in range(0,d):
                xp = xp + C[r+1,j+1] * Xkj[str(r)]
            #Append collocation equations
            fj = sf(Xkj[str(j)],Uk)
            cons = vertcat(cons, h*fj - xp)

            #Add contribution to the end state
            Xk_end = Xk_end + D[j+1]*Xkj[str(j)]

        #New NLP variable for state at end of interval
        Xk = MX.sym('X_'+ str((iter)*nk + k), nx)
        V = vertcat(V, Xk)
        #Add equality constraint
        cons = vertcat(cons, Xk_end-Xk)

        Jecon = (pf*F_0 + pV*Uk[1] - pB*Uk[4] - pD*Uk[3])*delta_t
        Jstate =  mtimes(transpose(multiply(Qmax[0:nx],(Xk -xdot_val_rf_ss))),
                 (Xk - xdot_val_rf_ss))*delta_t

        #Compute rotate cost function
        fm = sf(Xk, Uk)
        alpha = 1
        beta = 1
        gamma = 1
        
        obj = obj + alpha*Jcontrol + gamma*Jstate + beta*Jecon

    return obj, cons, V, Xk, params, ssoftc
