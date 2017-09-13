"""
================================================================================================
                                 Test for iNMPC module
                       Case: car position control - Grune's book eg 2.3         
                        Federico Lozano Santamaria - November 2015                           
================================================================================================
"""

## Imports
from __future__ import division
from pyomo.environ import *
from pyomo.dae import *
from pyomo.iNMPC import *

 
# Model
m = ConcreteModel()
 
# Sets
m.t = ContinuousSet(bounds=(0, 1.0)) 

# Variables
m.x = Var(m.t, bounds=(None,None), initialize=0.0)
m.y = Var(m.t, bounds=(None,None), initialize=0.0)
m.ux = ManipulatedVar(m.t, bounds=(None,None), initialize=0.02)
m.uy = ManipulatedVar(m.t, bounds=(None,None), initialize=0.0)
 
# Derivative variables
m.dxdt = DerivativeVar(m.x, wrt=m.t)
m.dydt = DerivativeVar(m.y, wrt=m.t)	

# Initial conditions
m.x_init = InitialCondition(m.x, time_set=m.t, initialize=0.0) 
m.y_init = InitialCondition(m.y, time_set=m.t, initialize=0.0)

# Set-point change
def sp_change(c_time):
    if c_time <= 150:
        return 1
    else:
    	return 0
m.y_sp = Disturbance(m.t, rule=sp_change, forecast='CONSTANT')    	

## Constraints

#Rail
def _c1(m, l):    
    return (m.x[l]**2) + 4*((m.y[l]-0.5)**2) == 1    
m.C1 = Constraint(m.t, rule=_c1)    

#Max velocity
def _vn(m, l):
    return (m.ux[l]**2) + (m.uy[l]**2) <= (0.02)**2
m.vn = Constraint(m.t, rule=_vn)    


#Differential equations
def _dx(m, l):
    if l == 0:
        return m.x[0] == m.x_init
    else:
    	return m.dxdt[l] == m.ux[l]
m.dx = Constraint(m.t, rule=_dx)

def _dy(m, l):
    if l == 0:
        return m.y[0] ==  m.y_init
    else:
    	return m.dydt[l] == m.uy[l]
m.dy = Constraint(m.t, rule=_dy)
 
## Objective
def _obj(m):              
    return sum( (m.x[l]-0.0)**2 + 5*(m.y[l]-m.y_sp[l])**2 + (m.ux[l]**2) + (m.uy[l]**2) for l in m.t)
m.OBJ = Objective(rule=_obj, sense=minimize)
 
## NMPC
CONTROLLER = IdealNMPC(model=m, model_time=m.t, prediction_horizon=4,
    control_horizon=4, simulation_periods=300, default_discretization=False)
CONTROLLER.discretize_model(discretization_type='FiniteDifferences', 
				wrt=m.t, scheme='BACKWARD')
CONTROLLER.solve('ipopt', options=[('halt_on_ampl_error', 'yes'),
				   ('linear_solver','ma27')],
    file_name='car_iNMPC_test', dynamic_plot=True)
