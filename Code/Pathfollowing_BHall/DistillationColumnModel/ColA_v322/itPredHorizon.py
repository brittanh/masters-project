#!/opt/local/bin/python
# -*- encoding: ascii -*-
"""
    @purpose: solving optimal control problem
    @author: Brittany Hall
    @date: 07.10.2017
    @version: 0.1
    @updates:
"""
from casadi import *
from numpy import ones, zeros, multiply, append
import scipy.io as spio

def itPredHorizon(Xk, w, w0, lbw, ubw, lbg, ubg, g, J, params, iter, count, ssoftc, d):
    
    #extracting parameter variables
    nx = params['prob']['nx']           #Number of states
    nu = params['prob']['nu']           #Number of inputs
    nk = params['prob']['nk']
    tf = params['prob']['tf']
    h =  params['prob']['h']
    ns = params['prob']['ns']
    
    x_min = params['bounds']['x_min']
    x_max = params['bounds']['x_max']
    u_min = params['bounds']['u_min']
    u_max = params['bounds']['u_max']
    
    NT = params['model']['NT']
    f = params['model']['f']
    xdot_val_rf_ss = params['model']['xdot_val_rf_ss']
    x = params['model']['x']
    u = params['model']['u']
    u_opt = params['model']['u_opt']
    
    pf = params['price']['pf']
    pV = params['price']['pV']
    pB = params['price']['pB']
    pD = params['price']['pD']
    
    F_0 = params['dist']['F_0']
    
    MDs = params['gain']['MDs']
    MBs = params['gain']['MBs']
    Ds = params['gain']['Ds']
    Bs = params['gain']['Bs']
    
    C = params['colloc']['C']
    D = params['colloc']['D']
    h = params['colloc']['h']
    
    delta_t = params['weight']['delta_t']
    alpha = params['weight']['alpha']
    beta = params['weight']['beta']
    gamma = params['weight']['gamma']
    Qmax = params['Qmax']
    
    for k in range(0,nk):
        #New NLP variable for control
        Uk = MX.sym('U_'+str((iter)*nk+k),nu)
        w = vertcat(w,Uk)
        lbw = append(lbw,u_min)
        ubw = append(ubw,u_max)
        
        indexU = iter*nk + k
        w0 = append(w0,u[:,indexU])
        Jcontrol = mtimes(transpose(multiply(Qmax[nx:nx+nu],
                                    Uk - u_opt)), (Uk - u_opt))
        
        #State at collocation points
        SumX1 = 0
        Xkj = {}
        for j in range(0,d):
            Xkj[str(j)] = MX.sym('X_' + str((iter)*nk + k)
                                 +'_'+str(j+1), nx)
            w = vertcat(w, Xkj[str(j)])
            lbw = append(lbw, x_min)
            ubw = append(ubw, x_max)
            w0 = append(w0, x[iter+1,:])
            count += 1
    
        #Loop over collocation points
        Xk_end = D[0] * Xk
        for j in range(0,d):
            xp = C[0,j+1] * Xk
            for r in range(0,d):
                xp = xp + C[r+1,j+1] * Xkj[str(r)]
            
            #Append collocation equations
            fj = f(Xkj[str(j)],Uk)
            g = vertcat(g, h*fj-xp)
            lbg = append(lbg, zeros((nx,1)))
            ubg = append(ubg, zeros((nx,1)))
            #Add contribution to the end state
            Xk_end = Xk_end + D[j+1]*Xkj[str(j)]
        
        #New NLP variable for state at end of interval
        Xk = MX.sym('X_'+ str((iter)*nk + k), nx)
        w = vertcat(w, Xk)
        lbw = append(lbw, x_min)
        x_maxEnd = ones((2*NT+2,1))
        x_maxEnd[0,0] = 0.1
        x_maxEnd[2*NT+1,0] = 0.7
        ubw = append(ubw, x_maxEnd)
        w0 = append(w0, x[iter+1,:])
        w0 = w0.reshape(len(w0),1)
        count += 1
        
        #Add equality constraint
        g = vertcat(g, Xk_end-Xk)
        lbg = append(lbg, zeros((nx,1)))
        ubg = append(ubg, zeros((nx,1)))
        
        Jecon = (pf*F_0 + pV*Uk[1] - pB*Uk[4]
                 - pD*Uk[3]) * delta_t
        Jstate = mtimes(transpose(multiply(Qmax[0:nx],
                (Xk -xdot_val_rf_ss))),(Xk - xdot_val_rf_ss))*delta_t

        J = J + alpha*Jcontrol + gamma*Jstate + beta*Jecon
    
    return J, g, w0, w, lbg, ubg, lbw, ubw, Xk, params, count, ssoftc
