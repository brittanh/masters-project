#!/opt/local/bin/python
# -*- encoding: ascii -*-
"""
    @purpose: Distillation column and CSTR model to be used in pathfollowing
              method
    @author: Brittany Hall
    @date: 09.11.2017
    @version: 0.1
    @updates:
"""
from numpy import zeros
from objective import *

def ColCSTR_pf(p):
    prob = {'neq': 0, 'niq': 0, 'cin': 0, 'ceq': 0, 'dp_in':0, 'dp_eq':0,
            'hess':0, 'lxp':0, 'x':0, 'name':0}

    prob['neq'] = 2000                          #Number of equality constraints
    prob['niq'] = 0                           #Number of inequality constraints
    prob['name'] = 'Distillation Column A + CSTR Model'
    prob['x'] = zeros((2,1))
    prob['obj'] = lambda x,y,p,N: objective(x,y,p,N)

    return prob
