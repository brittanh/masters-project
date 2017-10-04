#!/opt/local/bin/python
# -*- encoding: ascii -*-
"""
    @purpose: NLP solver
    @author: Brittany Hall
    @date: 18.09.2017
    @version: 0.1
    @updates:
"""
from casadi import *

def nlp_solve(prob, obj, p_init, x_init, y_init):
    """
    NLP solver for initial conditions to path-following algorithm
    """
    n, np, neq, niq, name = prob()
    if niq >0:
        x, p, f, f_fun, con, conf, ubx, lbx, ubg, lbg = obj(x_init, y_init, p_init, neq, niq, n, np)
        
        #Formulating NLP to solve

        nlp ={'x':x, 'p':p, 'f':f, 'g':con}
        solver = nlpsol('solver', 'ipopt', nlp)
        sol = solver(x0 = x_init, p = p_init,
                    lbg = lbg, ubg = ubg, ubx = ubx, lbx = lbx)
        x_opt = sol['x']                                         #solving for x
        lagm = sol['lam_x']
        lam_opt = lagm[0]
        mu_opt = lagm[1]
        #print('x_opt:',x_opt,'lambda:',lam_opt,'mu:',mu_opt)
        return(x_opt, lam_opt, mu_opt, con)
