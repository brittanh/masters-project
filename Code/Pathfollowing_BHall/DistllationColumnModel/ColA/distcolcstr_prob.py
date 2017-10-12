#!/opt/local/bin/python
# -*- encoding: ascii -*-
"""
    @purpose: Formatting CSTR and distillation column problem
    @author: Brittany Hall
    @date: 08.10.2017
    @version: 0.1
    @updates:
"""

from casadi import *
from numpy import zeros, append, tranpose
from objective import *

def itPredHor(Xk, V, cons, obj, params, iter, ssoftc):
    NT = params['model']['NT']
    sf = params['model']['sf']
    xdot_val_rf_ss = params['model']['xdot_val_rf_ss']
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
    Qmax = params['Qmax']

    global nx nu nk d ns N
    count = 1
    for k in range(0,nk-1):
        Uk = MX.sym('U_' str((iter-1)*nk+k), nu)
        V = append(V, Uk)
        Jcontrol = transpose((Qmax[nx+1:nx+nu,0].* (Uk-u_opt)))*(Uk - u_opt)
        #state at collocation points
        Xkj = array([])
        SumX1 = 0
        for j in range(0,d):
            Xkj[j] = MX.sym('X_'+str((iter-1)*nk+k)+'_'+str(j),nx)
            V = append(V[:],Xkj[j])
        #Loop over collocation points
        Xk_end = D[0]*Xk
        for j in range(0,d):
            xp = C[0,j+1]*Xk
            for r in range[0,d]:
                xp = xp + C[r+1, j+1]*Xkj[r]
            #append collocation equations
            fj = sf[Xkj[j],Uk]
            cons = append(cons[:], h*fj - xp)
            #Add contribution to the end state
            Xk_end = Xk_end + D[j+1]*Xkj[j]

        Xk = MX.sym('X_'+str((iter-1)*nk+k, nx)
        V = append(V[:],Xk)

        #Add equality constraint
        cons = append(cons, Xk_end-Xk)
                    
        Jecon = (pf*F_0 + pV*Uk[1] - pB*Uk[4] - pD*Uk[3])*delta_t
                    
        Jstate = transpose(Qmax[0:nx,0].*(Xk - xdot_val_rf_ss))*(Xk - xdot_val_rf_ss)*delta_t
        
        #Compute rotated cost function
        fm = sf(Xk,Uk)
        #Load Lagrange multipliers from steady-state optimization
        data = spio.loadmat('LambdaCstrDist.mat', squeeze_me = True)
        lam = data['lambda']
        Jmodel = transpose(lam)*fm

        alpha = params['weight']['alpha']
        beta = params['weight']['beta']
        gamma = params['weight']['gamma']
                    
        opt = opt + alpha*Jcontrol + gamma*Jstate + beta*Jecon
        return obj, cons, V, Xk, params, ssoftc
                    
                    
def distcolcstr_prob(p):
    prob = {'neq':{0},'niq':{0},'cin':{0}, 'ceq':{0}, 'dp_in':{0},'dp_eq':{0}, 'hess':{0}, 'lxp':{0},'x': 0, 'name': 0}
    prob['neq'] = 2000                          #Number of equality constraints
    prob['niq'] = 0                           #Number of inequality constraints
    prob['name'] = 'CSTR + Distillation Column Model'
    prob['x'] = zeros((2,1))
    prob['obj'] = objective(x,y,p,N)



