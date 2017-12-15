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

def nlp_solve(prob, options, w0, lbw, ubw, lbg, ubg):
    """
    NLP solver for initial conditions to path-following algorithm
    """
    #Formulating NLP to solve
    solver = nlpsol('solver', 'ipopt', prob, options)
    sol = solver(x0 = w0, lbx = lbw, ubx = ubw, lbg=lbg, ubg=ubg)
    return sol
