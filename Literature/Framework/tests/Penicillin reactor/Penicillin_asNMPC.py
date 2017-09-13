"""
================================================================================================
                                 Test for asNMPC
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
from pyomo.asNMPC import *

 
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
CONTROLLER = AdvancedStepNMPC(model=m, model_time=m.t, prediction_horizon=50,
    control_horizon=25, simulation_periods=50, default_discretization=True)

def plant(ManVar_dict, IC_dict, current_time):
    """
    A function of the real plant. This method must return a dictionary between
    a string corresponding to the intial conditions name and their value
    """
    #Model
    p = ConcreteModel()
    p.t = ContinuousSet(bounds=(0, 0.02))
    p.y1 = Var(p.t)
    p.y2 = Var(p.t)
    p.dy1dt = DerivativeVar(p.y1)
    p.dy2dt = DerivativeVar(p.y2)     
    theta = ManVar_dict['theta'][0].value
    y1_init = IC_dict['y1_init'].value
    y2_init = IC_dict['y2_init'].value
    def _b1(p,i):
        return 15.1*( (1-(0.005*((theta-30)**2))) / (1-(0.005*((25-30)**2))) )
    p.b1 = Expression(p.t, rule=_b1)	 
    def _b2(p,i):
        return 0.74*( (1-(0.005*((theta-30)**2))) / (1-(0.005*((25-30)**2))) )
    p.b2 = Expression(p.t, rule=_b2)	 
    def _Dy1(p,i):  
        if i == 0:
	    return p.y1[i] == y1_init
	else:
	    return p.dy1dt[i] == (p.b1[i]*p.y1[i]) - ((p.b1[i]/p.b2[i])*(p.y1[i]**2))    
    p.EqDy1 = Constraint(p.t,rule=_Dy1)	 
    def _b3(p,i):
        return 1.91*( (1-(0.005*((theta-20)**2))) / (1-(0.005*((25-20)**2))) )
    p.b3 = Expression(p.t, rule=_b3)	 
    def _Dy2(p,i):  
        if i == 0:
            return p.y2[i] == y2_init
        else:
	    return p.dy2dt[i] == p.b3[i]*p.y1[i]
    p.EqDy2 = Constraint(p.t,rule=_Dy2)    
    #def _obj(p):      
    #    return -p.y2[p.t[-1]]  
    #p.OBJ = Objective(rule=_obj,sense=minimize)  
    #Discretize
    from pyomo.dae.plugins.colloc import Collocation_Discretization_Transformation
    discretize = Collocation_Discretization_Transformation()
    disc = discretize.apply(p,nfe=1,ncp=3,scheme='LAGRANGE-RADAU')
    #Solution
    from pyomo.opt import SolverFactory
    opt = SolverFactory('ipopt')
    results = opt.solve(p, tee=False)    
    #Dictionary file - THIS DICTIONARY MUST CONTAINT ALL THE INITIAL CONDITIONS
    #IT MUST INCLUDE ALL THE INDEXES
    return_dict = {}
    return_dict['y1_init'] = p.y1[0.02].value
    return_dict['y2_init'] = p.y2[0.02].value        
    return return_dict

CONTROLLER.solve('ipopt_sens', options=[('linear_solver', 'ma57'), 
		('mu_init','1E-4')], file_name='Penicillin', real_plant = plant, 
		dynamic_plot=True)
