"""
================================================================================================
                       Test for iNMPC module using nl interface
                           Case: simple penicillin reactor         
                       Federico Lozano Santamaria - November 2015                           
================================================================================================
"""

## Imports
from __future__ import division
from pyomo.environ import *
from pyomo.dae import *
import math
 
## NEW IMPORTS FOR (N)MPC
from pyomo.iNMPC_nl_interface import *

 
## Model
m = ConcreteModel()
 
## Sets
m.to = 0
m.tf = 0.02
m.t = ContinuousSet(bounds=(m.to, m.tf)) # This set constains zero and the sampling time

m.y1 = Var(m.t, within=NonNegativeReals,bounds=(0,None),initialize=0.8)
m.y2 = Var(m.t, within=NonNegativeReals,bounds=(0,None),initialize=0.0)
m.theta = ManipulatedVar(m.t,within=NonNegativeReals,bounds=(20,30),initialize=27)
 
# Initial conditions
m.y1_init = InitialCondition(m.y1, time_set=m.t, within=NonNegativeReals, initialize=0.02) 
m.y2_init = InitialCondition(m.y2, time_set=m.t, within=NonNegativeReals, initialize=0.0)
 
## Derivative variables
m.dy1dt = DerivativeVar(m.y1,wrt=m.t)
m.dy2dt = DerivativeVar(m.y2,wrt=m.t)
 
## Constraints
def _init_y1(m):
  return m.y1[m.to] == m.y1_init
m.EqInity1 = Constraint(rule=_init_y1)
def _init_y2(m):
  return m.y2[m.to] ==  m.y2_init
m.EqInity2 = Constraint(rule=_init_y2)
 
def _b1(m,i):
    return 13.1*( (1-(0.005*((m.theta[i]-30)**2))) / (1-(0.005*((25-30)**2))) )
m.b1 = Expression(m.t, rule=_b1)
 
def _b2(m,i):
    return 0.94*( (1-(0.005*((m.theta[i]-30)**2))) / (1-(0.005*((25-30)**2))) )
m.b2 = Expression(m.t, rule=_b2)
 
def _Dy1(m,i):  
  if i == m.to:
    return Constraint.Skip
  else:
    return m.dy1dt[i] == (m.b1[i]*m.y1[i]) - ((m.b1[i]/m.b2[i])*(m.y1[i]**2))    
m.EqDy1 = Constraint(m.t,rule=_Dy1)
 
def _b3(m,i):
    return 1.71*( (1-(0.005*((m.theta[i]-20)**2))) / (1-(0.005*((25-20)**2))) )
m.b3 = Expression(m.t, rule=_b3)
 
def _Dy2(m,i):  
  if i == m.to:
    return Constraint.Skip
  else:
    return m.dy2dt[i] == m.b3[i]*m.y1[i]
m.EqDy2 = Constraint(m.t,rule=_Dy2)
 
## Objective
def _obj(m):      
  return -m.y2[m.t[-1]]  
m.OBJ = Objective(rule=_obj,sense=minimize)
 
## NMPC
CONTROLLER = IdealNMPC(model=m, model_time=m.t, prediction_horizon=50,
    control_horizon=25, simulation_periods=50, default_discretization=True)
CONTROLLER.solve('ipopt', options=[('linear_solver', 'ma57')],
    file_name='Penicillin_iNMPC_nl_test', dynamic_plot=True)
