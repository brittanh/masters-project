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
from numpy import append, ones, zeros, multiply
import scipy.io as spio

def itPredHorizon(Xk, w, w0, lbw, ubw, lbg, ubg, g, J, params, i, count, ssoftc):
    global N, nx, nu, nk, d, ns
    #extracting parameter variables
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
    F_0 = params['price']['F_0']
    
    KcB = params['gain']['KcB']
    KcD = params['gain']['KcD']
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

    for k in range(0,nk-1):
        #New NLP variable for control
        Uk = Mx.sym(['U_'+str((i-1)*nk+k)],nu)
        w = append(w,Uk)
        lbw = append(lbw, u_min)
        ubw = append(ubw, u_max)
        indexU = (i-1)*nk + (k+1)
        w0 = append(w0, u[:,indexU])
        
        Jcontrol = transpose(multiply(Qmax[nx+1:nx+nu, 1],(Uk-u_opt)))*(Uk - u_opt)
        
        #State at collocation points
        Xkj = arrays([])
        SumX1 = 0
        for j in range(0,d):
            Xkj[j] = MX.sym(['X_'+str((i-1)*nk+k)+'_'+str(j)],nx)
            w = append(w,Xkj[j])
            lbw = append(lbw, x_min)
            ubw = append(ubw, x_max)
            count += 1
        #Loop over collocation points
        XK_end = D[0] * Xk
        for j in range(0,d):
            xp = C[0,j+1] * Xk
            for r in range(0,d):
                xp = xp + C[r+1,j+1]*Xkj[r]
            #Append collocation equations
            fj = f[Xkj[j],Uk]
            g = append(g, h*fj-xp)
            lbg = append(lbg, zeros((nx,1)))
            ubg = append(ubg, zeros((nx,1)))
            #Add contribution to the end state
            Xk_end = Xk_end + D[j+1]*Xkj[j]
        
            #New NLP variable for state at end of interval
            Xk = Mx.sym(['X_'+str((i-1)*nk+k)],nx)
            w = append(w, Xk)
            lbw = append(lbw, x_min)
            x_maxEnd = ones((84,1))
            x_maxEnd[0,0] = 0.1
            x_maxEnd[84,0] = 0.7
            ubw = append(ubw, x_maxEnd)
            w0 = append(w0, transpose(x[i+1,:]))
            count += 1
            
            #Add equality constraint
            g = append(g, Xk_end-Xk)
            lbg = append(lbg, zeros((nx,1)))
            ubg = append(ubg, zeros((nx,1)))
            
            Jecon = (pf*F_0 + pV*Uk[1] - pB*Uk[4] - pD*Uk[3])*delta_t
            Jstate = transpose(multiply(Qmax[0:nx,0],(Xk -xdot_val_rf_ss)))*(Xk - xdot_val_rf_ss)*delta_t
        #Compute rotated cost function
        fm = f[Xk, Uk]
           
        #Load Lagrange multipliers from steady-state optimization
        data = spio.loadmat('LambdaCstrDist.mat', squeeze_me = True)
        lam = data['lambda']
        Jmodel = transpose(lam)*fm
           
        J = J + alpha*Jcontrol + gamma*Jstate + beta*Jecon
    return J, g, w0, w, lbg, ubg, lbw, ubw, Xk, params, count, ssoftc

