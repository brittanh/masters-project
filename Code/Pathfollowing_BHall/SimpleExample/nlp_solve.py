#!/opt/local/bin/python
# -*- encoding: ascii -*-
"""
    @purpose: NLP solver
    @author: Brittany Hall
    @date: 18.09.2017
    @version: 0.1
    @updates:
"""
from casadi import nlpsol

def nlp_solve(prob, obj, p_init, x_init, y_init):
    """
    NLP solver for initial conditions to path-following algorithm
    """
    nx, np, neq, niq, name = prob()
    if niq >0:
        x, p, f, f_fun, con, conf, ubx, lbx, ubg, lbg = obj(x_init, y_init, p_init, neq, niq, nx, np)
        
        #Formulating NLP to solve
        #All constraints must be formatted as inequality constraints for this solver
        nlp = {'x':x, 'p':p, 'f':f, 'g':con}
        solver = nlpsol('solver', 'ipopt', nlp)
        sol = solver(x0 = x_init, p = p_init,
                    lbg = lbg, ubg = ubg, ubx = ubx, lbx = lbx)
        x_opt = sol['x']                                         #Solving for x
        lagmul = sol['lam_g']
        #Determining active constraints
        #(necessary to determine which multipliers are a lambda and which are a mu)
        con_vals = conf(x_opt,p_init)
        tol = 1e-6
        for k in range(0,len(con_vals)):
            if con_vals[k] >= 0 + tol or con_vals[k] >= 0 - tol: #active constraint
                lam_opt = lagmul[k]
            else: #inactive constraint
                mu_opt = lagmul[k]
        #print('x_opt:',x_opt,'lambda:',lam_opt,'mu:',mu_opt)
        return x_opt, lam_opt, mu_opt, con 
