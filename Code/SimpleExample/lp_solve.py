#!/opt/local/bin/python
# -*- encoding: ascii -*-
"""
    @purpose: LP solver to obtain dual variables for path following algorithm
    @author: Brittany Hall
    @date: 18.09.2017
    @version: 1
    @updates:
"""
import sys
sys.path.append('/Library/gurobi751/mac64/lib/python2.7/gurobipy/')
from gurobipy import *

#Creating model
m = Model("lp_solve")

#Creating variables
x = m.addVariable(vtype = GRB.BINARY, name="x")
