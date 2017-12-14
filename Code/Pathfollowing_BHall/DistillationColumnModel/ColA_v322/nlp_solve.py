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

def nlp_solve(problem, options, x0, lbx, ubx, lbg, ubg):
    """
    NLP solver for initial conditions to path-following algorithm
    """
    #Formulating NLP to solve
    solver = nlpsol('solver', 'ipopt', problem, options)
    sol = solver(x0=x0, lbx=lbx, ubx=ubx, lbg=lbg, ubg=ubg)
    return sol
