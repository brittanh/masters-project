"""
================================================================================================
                                   Test for asNMPC
            Case: extractive distillation column of ethanol-water system
         Index 2 DAE model reduced to index 1 using differentiation and subtitution
                       Federico Lozano Santamaria - November 2015                           
================================================================================================
"""

from __future__ import division

## ================================================================================================
## Solution of steady state for initial conditions
## ================================================================================================
from SteadyState import ms
ss = ms.create_instance('SteadyState.dat')
from pyomo.opt import SolverFactory
opt = SolverFactory('ipopt')
opt.options['Linear_Solver'] = 'ma57'
ss.preprocess()
results = opt.solve(ss, tee=True)

## ================================================================================================
## Imports and model creation
## ================================================================================================
from pyomo.environ import *
from pyomo.dae import *
from pyomo.asNMPC import *

m = AbstractModel()

## ================================================================================================
## Sets
## ================================================================================================
m.comp = Set()
m.NT = Param(within=PositiveIntegers) #Total number of stages
m.Net = RangeSet(1,m.NT)

m.tf = 5.0
m.to = 0
m.t = ContinuousSet( bounds=(m.to, m.tf) )

## ================================================================================================
## Parameters
## ================================================================================================

# Operation
m.FA = Param(m.Net, within=NonNegativeReals, initialize=0, default=0, mutable=True) # Azeotropic mix feed flow rate [kmol/min]. Disturbance
m.TfeedA = Param(within=PositiveReals) # Azeotropic mixture feed temperature

m.FG = Param(m.Net, within=NonNegativeReals, initialize=0) # Glycerol feed flow rate [kmol/min]
m.zg = Param(m.comp, within=NonNegativeReals) # Glycerol feed composition
m.TfeedG = Param(within=PositiveReals) # Glycerol feed temperature

m.xD_Ethanol = Param(within=PositiveReals) # Ethanol composition in distillate

m.P = Param(m.Net, within=NonNegativeReals) # Preassure [bar]

m.ML_C = Param(within=PositiveReals) # Molar liquid hold-up in condenser [kmol]
m.ML_R = Param(within=PositiveReals) # Molar liquid hold-up in reboiler [kmol]

# Tray and column geometry
m.CD = Param(within=PositiveReals) # Column diameter [m]
m.Lw = Param(within=PositiveReals) # Weir length [m]
m.hw = Param(within=PositiveReals) # Weir height [m]
m.AC = Param(within=PositiveReals) # Column area [m2]
m.AT = Param(within=PositiveReals) # Active area [m2]
m.A0 = Param(within=PositiveReals) # Holey area [m2]
m.da = Param(within=PositiveReals) # Hole diameter [m]
m.ep = Param(within=PositiveReals) # Tray thikness [m]
m.pitch = Param(within=PositiveReals) # Distance between holes [m]
m.HS = Param(within=PositiveReals) # Distance between trays [m]

## ================================================================================================
## Variables
## ================================================================================================
# Original variables
def i_ML0(m,l,i):
  return ss.ML[i].value
m.ML = Var(m.t, m.Net, bounds=(0,5), initialize=i_ML0) # [kmol]
def i_L0(m,l,i):
  return ss.L[i].value/60
m.L = Var(m.t, m.Net, bounds=(0,5), initialize=i_L0) # [kmol/min]
def i_V0(m,l,i):
  return ss.V[i].value/60
m.V = Var(m.t, m.Net, bounds=(0,5), initialize=i_V0) # [kmol/min]
def i_x0(m,l,i,j):
  return ss.x[i,j].value
m.x = Var(m.t, m.Net, m.comp, bounds=(0,100), initialize=i_x0) # 0-100
def i_y0(m,l,i,j):
  return ss.y[i,j].value
m.y = Var(m.t, m.Net, m.comp, bounds=(0,100), initialize=i_y0) # 0-100
def i_Temp0(m,l,i):
  return ss.Temp[i].value
m.Temp = Var(m.t, m.Net, bounds=(273,500), initialize=i_Temp0) # [K]
m.D = Var(m.t, bounds=(0,5), initialize=ss.D.value/60) # [kmol/min]
m.Qc = Var(m.t, bounds=(0,200), initialize=ss.Qc.value/60) # [MJ/min]
m.Qr = ManipulatedVar(m.t, bounds=(0,200), initialize=ss.Qr.value/60) # [MJ/min]
m.LR = ManipulatedVar(m.t, bounds=(0,5), initialize=ss.L[1].value/60) # [kmol/min]

# Derivatives variables
m.dMLdt = DerivativeVar(m.ML, wrt=m.t) 
m.dxdt = DerivativeVar(m.x, wrt=m.t) 

# Initial conditions
def ic_x0(m,i,j):
  return ss.x[i,j].value
m.x_init = InitialCondition(m.x, time_set=m.t, within=NonNegativeReals, initialize=ic_x0) # Initial conditions for liquid mole fractions
def ic_ML0(m,i):
  return ss.ML[i].value
m.ML_init = InitialCondition(m.ML, time_set=m.t, within=NonNegativeReals, initialize=ic_ML0) # Initial condition for the molar liquid hold-up

## ================================================================================================
## Constants for physical properties
## ================================================================================================

# Constants for liquid density
m.MW = Param(m.comp)
m.C1r = Param(m.comp)
m.C2r = Param(m.comp)
m.C3r = Param(m.comp)
m.C4r = Param(m.comp)
m.C5r = Param(m.comp)

# Constants for saturation preassure
m.C1sp = Param(m.comp)
m.C2sp = Param(m.comp)
m.C3sp = Param(m.comp)
m.C4sp = Param(m.comp)
m.C5sp = Param(m.comp)
m.C6sp = Param(m.comp)
m.C7sp = Param(m.comp)

# Constants for NRTL
m.a = Param(m.comp,m.comp)
m.b = Param(m.comp,m.comp)
m.c = Param(m.comp,m.comp)

# Constants for vapor enthalpy
m.H0 = Param(m.comp)
m.C1c = Param(m.comp)
m.C2c = Param(m.comp)
m.C3c = Param(m.comp)
m.C4c = Param(m.comp)
m.C5c = Param(m.comp)
m.C6c = Param(m.comp)

# Constants for liquid enthalpy
m.Tcrit = Param(m.comp)
m.C1v = Param(m.comp)
m.C2v = Param(m.comp)
m.C3v = Param(m.comp)
m.C4v = Param(m.comp)
m.C5v = Param(m.comp)

## ================================================================================================
## Construct the model
## ================================================================================================
m = m.create_instance('All_DAE2r.dat') 

## ================================================================================================
## Disturbance
## ================================================================================================
m.za = {}

def _za_e(c_time):
  return 80 - 5.0*sin(c_time*3.14159/60)
m.za_e = Disturbance(m.t, rule=_za_e, forecast='CONSTANT')
def _za_w(c_time):
  return 20 + 5.0*sin(c_time*3.14159/60)
m.za_w = Disturbance(m.t, rule=_za_w, forecast='CONSTANT')  
def _za_g(c_time):
  return 0
m.za_g = Disturbance(m.t, rule=_za_g, forecast='CONSTANT')  

m.za['Ethanol'] = m.za_e
m.za['Water'] = m.za_w
m.za['Glycerol'] = m.za_g

## ================================================================================================
## Reflux flow rate
## ================================================================================================

def _EqLR(m,l):
  return m.L[l,1] == m.LR[l]
m.RefluxFlowRate = Constraint(m.t, rule=_EqLR)

## ================================================================================================
## Purity constraint
## ================================================================================================
#'''
m.slack_var_x = Var(m.t, bounds=(0,None), initialize=0)
m.penalty = Var(bounds=(0,None))
m.s_cons = Constraint(expr = m.slack_var_x[0] == 0)
def _purity(m,k,i,j):
  if k != m.to and i == 1 and j == 'Ethanol':
    return m.x[k,i,j] == m.xD_Ethanol + m.slack_var_x[k] - m.penalty
  else:
    return Constraint.Skip
m.Purity = Constraint(m.t, m.Net, m.comp, rule=_purity)
#'''
## ================================================================================================
# Liquid flow rate
## ================================================================================================

# Liquid density equation
def _EqRhoL(m,l,i):
  rho = {}
  for j in m.comp:
    if j == 'Water':
      Tr = 1 - ((m.Temp[l,i]-273.15)/373.946)
      rho[j] = m.C1r[j] + (m.C2r[j]*(Tr**0.35)) + (m.C3r[j]*(Tr**(2.0/3.0))) + (m.C4r[j]*Tr) + (m.C5r[j]*(Tr**(4.0/3.0))) 
    else:
      rho[j] = m.C1r[j] / ( m.C2r[j]**( 1+ ((1-(m.Temp[l,i]/m.C3r[j]))**m.C4r[j]) ) )
  rhoL_mol = sum( (m.x[l,i,j]/100)*rho[j] for j in m.comp)    
  return rhoL_mol
m.rhoL = Expression(m.t, m.Net, rule=_EqRhoL)

# Vapor density equation
def _EqRhoV(m,l,i):
  return m.P[i] / (0.08314*m.Temp[l,i]) # R = 0.08314 m3*bar/K*kmol  
m.rhoV = Expression(m.t, m.Net, rule=_EqRhoV)

# Flow factpr
def _EqFP(m,l,i):
  if i == 1:
    return 0
  elif i == m.NT:
    return 0
  else:
    return ((m.L[l,i]/m.rhoL[l,i])/(m.V[l,i]/m.rhoV[l,i]))*(((m.rhoL[l,i]*sum((m.x[l,i,j]/100)*m.MW[j] for j in m.comp))/(m.rhoV[l,i]*sum((m.y[l,i,j]/100)*m.MW[j] for j in m.comp)))**(0.5)) 
m.FP = Expression(m.t, m.Net, rule=_EqFP)

# Liquid flow rate
def _LF(m,l,i):
  if i == 1: # Condenser (steady-state)
    return m.V[l,i+1] == m.L[l,i] + m.D[l]
  elif i == m.NT: # Reboiler (steady-state)
    return m.L[l,i-1] == m.V[l,i] + m.L[l,i]
  else: # Column    
    return m.ML[l,i] == 0.6*m.rhoL[l,i]*m.AT*(m.hw**0.5)*((m.FP[l,i]*m.AT*m.pitch/m.Lw)**0.25)
m.LiquidFlowRate = Constraint(m.t, m.Net, rule=_LF)

## ================================================================================================
## MESH equations
## ================================================================================================

# -------------------------------------------------------------------------------------------------
# Thermodynamic equilibrium
# -------------------------------------------------------------------------------------------------

# Saturation preassure
def _EqPsat(m,l,i,j):
  return exp( m.C1sp[j] + (m.C2sp[j]/( m.Temp[l,i]+m.C3sp[j])) + (m.C4sp[j]* m.Temp[l,i]) + (m.C5sp[j]*log( m.Temp[l,i])) + (m.C6sp[j]*( m.Temp[l,i]**m.C7sp[j])) ) 
m.Psat = Expression(m.t, m.Net, m.comp, rule=_EqPsat)

# Activity coefficient
def _EqTau(m,l,k,i,j):
  return m.a[i,j] + (m.b[i,j]/m.Temp[l,k])
m.Tau = Expression(m.t, m.Net, m.comp, m.comp, rule=_EqTau )

def _EqG(m,l,k,i,j):
  return exp(-m.c[i,j]*m.Tau[l,k,i,j])
m.G = Expression(m.t, m.Net, m.comp, m.comp, rule=_EqG )

def _EqS(m,l,k,i):
  return sum( (m.x[l,k,j]/100)*m.G[l,k,j,i] for j in m.comp )
m.S = Expression(m.t, m.Net, m.comp, rule=_EqS)

def _EqC(m,l,k,i):
  return sum( (m.x[l,k,j]/100)*m.G[l,k,j,i]*m.Tau[l,k,j,i] for j in m.comp )
m.C = Expression(m.t, m.Net, m.comp, rule=_EqC)

def _EqEps(m,l,j,i,k):
  return (m.G[l,j,i,k]/m.S[l,j,k])*(m.Tau[l,j,i,k] - (m.C[l,j,k]/m.S[l,j,k]))
m.epsilon = Expression(m.t, m.Net, m.comp, m.comp, rule=_EqEps )

def _Eqgamma(m,l,j,i):
  return exp( (m.C[l,j,i]/m.S[l,j,i]) + sum( (m.x[l,j,k]/100)*m.epsilon[l,j,i,k] for k in m.comp ) )
m.gamma = Expression(m.t, m.Net, m.comp, rule=_Eqgamma )

# Equation
def _EK(m,l,i,j):  
  return m.y[l,i,j]*m.P[i] == m.Psat[l,i,j]*m.gamma[l,i,j]*m.x[l,i,j]
m.ThermodynamicEquilibrium = Constraint(m.t, m.Net, m.comp, rule=_EK)

# -------------------------------------------------------------------------------------------------
# Phase equilibrium error
# -------------------------------------------------------------------------------------------------
def _EE(m,l,i):
  return sum( m.y[l,i,j] - m.x[l,i,j] for j in m.comp ) == 0
m.EquilibriumError = Constraint(m.t, m.Net, rule=_EE)

# -------------------------------------------------------------------------------------------------
# Total mass balance
# -------------------------------------------------------------------------------------------------

# Initial conditions  
def _ML0(m,j):
  return m.ML[m.to,j] == m.ML_init[j]#m.ML_0[j]
m.IC_ML = Constraint(m.Net, rule=_ML0)

# Equations
def _TMB(m,l,j):
  if l == m.to:
    return Constraint.Skip
  else:
    if j == 1: # Condenser
      return m.dMLdt[l,j] == 0
    elif j == m.NT: # Reboiler
      return m.dMLdt[l,j] == 0
    else: # Column
      return m.dMLdt[l,j] == m.FA[j] + m.FG[j] + m.V[l,j+1] + m.L[l,j-1] - m.V[l,j] - m.L[l,j]
m.TotalMassBalance = Constraint(m.t, m.Net, rule=_TMB)

# No vapor flow rate in condenser
def _V0(m,l,j):
  if j == 1:
    return m.V[l,j] == 0
  else:
    return Constraint.Skip
m.VaporCondenser = Constraint(m.t, m.Net, rule=_V0)

# -------------------------------------------------------------------------------------------------
# Partial mass balance
# -------------------------------------------------------------------------------------------------

# Initial conditions  
def _x0(m,i,j):
  return m.x[m.to,i,j] == m.x_init[i,j]#m.x_0[i,j]
m.IC_x = Constraint(m.Net, m.comp, rule=_x0)

# Equation
def _PMB(m,l,i,j):
  #if l == m.to:
  #  return Constraint.Skip
  #else:
  if i == 1: # Condenser
    return m.dxdt[l,i,j] == ((m.y[l,i+1,j]*m.V[l,i+1]) - ((m.L[l,i] + m.D[l])*m.x[l,i,j]))/m.ML[l,i]
  elif i == m.NT: # Reboiler
    return m.dxdt[l,i,j] == ( (m.x[l,i-1,j]*m.L[l,i-1]) - (m.x[l,i,j]*m.L[l,i]) - (m.y[l,i,j]*m.V[l,i]) )/m.ML[l,i]
  else: # Column
    return m.dxdt[l,i,j] == ( (m.FA[i]*(m.za[j][l]-m.x[l,i,j])) + (m.FG[i]*(m.zg[j]-m.x[l,i,j])) + (m.L[l,i-1]*(m.x[l,i-1,j]-m.x[l,i,j])) + (m.V[l,i+1]*(m.y[l,i+1,j]-m.x[l,i,j])) - (m.V[l,i]*(m.y[l,i,j]-m.x[l,i,j])) )/m.ML[l,i]
m.PartialMassBalance = Constraint(m.t, m.Net, m.comp, rule=_PMB)

# -------------------------------------------------------------------------------------------------
# Energy balances
# -------------------------------------------------------------------------------------------------

# Vapor enthalpy 
Tref = 298.15

def _EqHVi(m,l,i,j):
  return (m.C1c[j]*(m.Temp[l,i]-Tref)) + ((m.C2c[j]/2)*((m.Temp[l,i]**2)-(Tref**2))) + ((m.C3c[j]/3)*((m.Temp[l,i]**3)-(Tref**3))) + ((m.C4c[j]/4)*((m.Temp[l,i]**4)-(Tref**4))) + ((m.C5c[j]/5)*((m.Temp[l,i]**5)-(Tref**5))) + ((m.C6c[j]/6)*((m.Temp[l,i]**6)-(Tref**6))) + m.H0[j]
m.HVi = Expression(m.t, m.Net, m.comp, rule=_EqHVi) # [MJ/kmol]

def _EqHV(m,l,i):
  return sum( m.HVi[l,i,j]*m.y[l,i,j]/100 for j in m.comp)
m.HV = Expression(m.t, m.Net, rule=_EqHV ) # [MJ/kmol]

# Liquid enthalpy
def _EqDHvapi(m,l,i,j):
  Tr = m.Temp[l,i]/m.Tcrit[j]
  return m.C1v[j]*( (1-Tr)**(m.C2v[j] + (m.C3v[j]*Tr) + (m.C4v[j]*(Tr**2)) + (m.C5v[j]*(Tr**3)) ) )
m.DHvapi = Expression(m.t, m.Net, m.comp, rule=_EqDHvapi ) # [MJ/kmol]

def _EqHLi(m,l,i,j):
  return m.HVi[l,i,j] - m.DHvapi[l,i,j]
m.HLi = Expression(m.t, m.Net, m.comp, rule=_EqHLi ) # [MJ/kmol]

def _EqHL(m,l,i):
  return sum( m.HLi[l,i,j]*m.x[l,i,j]/100 for j in m.comp)
m.HL = Expression(m.t, m.Net, rule=_EqHL ) # [MJ/kmol]

# Derivatives for index reduction
def _EqdPsatdT(m,l,j,i):
  return m.Psat[l,j,i]*( (-m.C2sp[i]/((m.Temp[l,j]+m.C3sp[i])**2)) + m.C4sp[i] + (m.C5sp[i]/m.Temp[l,j]) + (m.C7sp[i]*m.C6sp[i]*(m.Temp[l,j]**(m.C7sp[i]-1))) )
m.dPsatdT = Expression(m.t, m.Net, m.comp, rule=_EqdPsatdT)



def _EqdgammadX(m,l,n,i,j):
  return m.gamma[l,n,i]*( m.epsilon[l,n,i,j] + m.epsilon[l,n,j,i] - sum( ((m.x[l,n,k]/100)/m.S[l,n,k])*( (m.G[l,n,i,k]*m.epsilon[l,n,j,k]) + (m.G[l,n,j,k]*m.epsilon[l,n,i,k]) ) for k in m.comp ) )
m.dgammadX = Expression(m.t, m.Net, m.comp, m.comp, rule=_EqdgammadX)



def _EqdTaudT(m,l,k,i,j):
  return -m.b[i,j]/(m.Temp[l,k]**2)
m.dTaudT = Expression(m.t, m.Net, m.comp, m.comp, rule=_EqdTaudT )

def _EqdGdT(m,l,k,i,j):
  return -m.c[i,j]*m.dTaudT[l,k,i,j]*m.G[l,k,i,j]
m.dGdT = Expression(m.t, m.Net, m.comp, m.comp, rule=_EqdGdT )

def _EqdSdT(m,l,k,i):
  return sum( (m.x[l,k,j]/100)*m.dGdT[l,k,j,i] for j in m.comp )
m.dSdT = Expression(m.t, m.Net, m.comp, rule=_EqdSdT)

def _EqdCdT(m,l,k,i):
  return sum( (m.x[l,k,j]/100)*( (m.G[l,k,j,i]*m.dTaudT[l,k,j,i]) + (m.Tau[l,k,j,i]*m.dGdT[l,k,j,i]) ) for j in m.comp )
m.dCdT = Expression(m.t, m.Net, m.comp, rule=_EqdCdT)

def _EqdEpsdT(m,l,j,i,k):
  return (1/(m.S[l,j,k]**2))*( (m.S[l,j,k]*( (m.G[l,j,i,k]*( m.dTaudT[l,j,i,k] - ( ( (m.S[l,j,k]*m.dCdT[l,j,k]) - (m.C[l,j,k]*m.dSdT[l,j,k]) )/(m.S[l,j,k]**2) ) )) + (m.dGdT[l,j,i,k]*(m.Tau[l,j,i,k]-(m.C[l,j,k]/m.S[l,j,k]))) )) - (m.G[l,j,i,k]*(m.Tau[l,j,i,k]-(m.C[l,j,k]/m.S[l,j,k]))*m.dSdT[l,j,k]) )
m.depsilondT = Expression(m.t, m.Net, m.comp, m.comp, rule=_EqdEpsdT )
 
def _EqdgammadT(m,l,j,i):
  return m.gamma[l,j,i]*( ( ( (m.S[l,j,i]*m.dCdT[l,j,i]) - (m.C[l,j,i]*m.dSdT[l,j,i]) )/(m.S[l,j,i]**2) ) + sum( (m.x[l,j,k]/100)*m.depsilondT[l,j,i,k] for k in m.comp ) )
m.dgammadT = Expression(m.t, m.Net, m.comp, rule=_EqdgammadT )



def _EqdTdt(m,l,j):
  return ( (((-1/m.P[j])*(sum( (m.x[l,j,i]/100)*m.Psat[l,j,i]*(sum( m.dgammadX[l,j,i,jj]*(m.dxdt[l,j,jj]/100) for jj in m.comp)) for i in m.comp)))) - sum( (m.Psat[l,j,i]*m.gamma[l,j,i]/m.P[j])*(m.dxdt[l,j,i]/100) for i in m.comp) ) / ( (1/m.P[j])*( sum((m.x[l,j,i]/100)*m.Psat[l,j,i]*m.dgammadT[l,j,i] for i in m.comp) + sum((m.x[l,j,i]/100)*m.gamma[l,j,i]*m.dPsatdT[l,j,i] for i in m.comp) ) )
m.dTdt = Expression(m.t, m.Net, rule=_EqdTdt)



def _EqdHVidT(m,l,j,i):
  return m.C1c[i] + (m.C2c[i]*m.Temp[l,j]) + (m.C3c[i]*(m.Temp[l,j]**2)) + (m.C4c[i]*(m.Temp[l,j]**3)) + (m.C5c[i]*(m.Temp[l,j]**4)) + (m.C6c[i]*(m.Temp[l,j]**5))
m.dHVidT = Expression(m.t, m.Net, m.comp, rule = _EqdHVidT)

def _EqdHvapidT(m,l,j,i):
  Tr = m.Temp[l,j]/m.Tcrit[i]
  return (m.DHvapi[l,j,i]/m.Tcrit[i])*( ((m.C3v[i] + (2*m.C4v[i]*Tr) + (3*m.C5v[i]*(Tr**2)))*log(m.C1v[i]*(1-Tr))) - ( (m.C2v[i] + (m.C3v[i]*Tr) + (m.C4v[i]*(Tr**2)) + (m.C5v[i]*(Tr**3)) )/(1-Tr) ) )
m.dHvapidT = Expression(m.t, m.Net, m.comp, rule = _EqdHvapidT)

def _EqdHLidT(m,l,j,i):
  return m.dHVidT[l,j,i] - m.dHvapidT[l,j,i]
m.dHLidT = Expression(m.t, m.Net, m.comp, rule = _EqdHLidT)

# Equation
from Enthalpies import Eq_HL
def _EB(m,l,j):
  za_t = {}
  if j == 1: # Condenser
    return m.ML[l,j]*( sum( m.HLi[l,j,i]*(m.dxdt[l,j,i]/100) for i in m.comp) + (m.dTdt[l,j]*sum((m.x[l,j,i]/100)*m.dHLidT[l,j,i] for i in m.comp)) )  == m.V[l,j+1]*m.HV[l,j+1] - (m.L[l,j] + m.D[l])*m.HL[l,j] - m.Qc[l]
  elif j == m.NT: # Reboiler
    return m.ML[l,j]*( sum( m.HLi[l,j,i]*(m.dxdt[l,j,i]/100) for i in m.comp) + (m.dTdt[l,j]*sum((m.x[l,j,i]/100)*m.dHLidT[l,j,i] for i in m.comp)) )  == m.L[l,j-1]*m.HL[l,j-1] + m.Qr[l] - m.V[l,j]*m.HV[l,j] - m.L[l,j]*m.HL[l,j]
  else: # Column
    for k in m.comp:
      za_t[k] = m.za[k][l]
    HFA = Eq_HL(m.comp, za_t, m.TfeedA, m.H0, m.C1c, m.C2c, m.C3c, m.C4c, m.C5c, m.C6c, m.Tcrit, m.C1v, m.C2v, m.C3v, m.C4v, m.C5v)
    HFG = Eq_HL(m.comp, m.zg, m.TfeedG, m.H0, m.C1c, m.C2c, m.C3c, m.C4c, m.C5c, m.C6c, m.Tcrit, m.C1v, m.C2v, m.C3v, m.C4v, m.C5v)      
    return m.ML[l,j]*( sum( m.HLi[l,j,i]*(m.dxdt[l,j,i]/100) for i in m.comp) + (m.dTdt[l,j]*sum((m.x[l,j,i]/100)*m.dHLidT[l,j,i] for i in m.comp)) )  == (m.FA[j]*(HFA-m.HL[l,j])) + (m.FG[j]*(HFG-m.HL[l,j])) + (m.V[l,j+1]*(m.HV[l,j+1]-m.HL[l,j])) + (m.L[l,j-1]*(m.HL[l,j-1]-m.HL[l,j])) - (m.V[l,j]*(m.HV[l,j]-m.HL[l,j]))
m.EnergyBalance = Constraint(m.t, m.Net, rule=_EB)

## ================================================================================================
## Objective
## ================================================================================================
def _obj(m):
  CostQr = 1E-3 # Cost of the reboiler duty [$/MJ]
  CostD = 30 # Cost of the distillate [$/kmol]
  alphaX = 1E4 # Weight for the tracking of ethanol purity in distillate
  alphaQr = 500 # Weight for the tracking of the reboiler duty
  alphaRR = 100 # Weight for the tracking of the reflux ratio	
  t_tmp = sorted(m.D.keys())
  f = {}
  for i in t_tmp:
    f[i] = (CostQr*m.Qr[i]) - (CostD*m.D[i]) + (alphaX*(( m.x[i,1,'Ethanol'] - 99.5 )**2)) + (alphaQr*(( m.Qr[i] - 93.2505 )**2)) + (alphaRR*(( (m.L[i,1]/m.D[i]) - 0.57906  )**2))     
  return (1E6*m.penalty) + sum( 0.5*(f[i]+f[t_tmp[t_tmp.index(i)-1]])*(i - t_tmp[t_tmp.index(i)-1]) for i in t_tmp if i>m.to)  
m.ObjectiveFunction = Objective(rule=_obj, sense=minimize)

## ================================================================================================
## asNMPC solution
## ================================================================================================
CONTROLLER = AdvancedStepNMPC(model=m, model_time=m.t, prediction_horizon=10,
    control_horizon=10, simulation_periods=60, default_discretization=True)

import random
def plant(ManVar_dict, IC_dict, current_time):
    """
    A function of the real plant. This method must return a dictionary between
    a string corresponding to the intial conditions name and their value
    """    

    #Dictionary file   
    return_dict = {}    
    for stage in m.Net:
    	if stage == 1 or stage == m.NT:
    	    return_dict['ML_init['+str(stage)+']'] = m.ML_init[stage].value
    	else:
            return_dict['ML_init['+str(stage)+']'] = m.ML_init[stage].value*(1 + random.uniform(-0.05, 0.05))
        for component in m.comp:
            return_dict['x_init['+str(stage)+','+str(component)+']'] = m.x_init[stage, component].value   
    return return_dict

CONTROLLER.solve('ipopt', options=[('linear_solver', 'ma57'), 
					   ('mu_init', 1E-5), 
					   ('warm_start_init_point', 'yes'),
					   ('warm_start_bound_push', 1e-5),
					   ('warm_start_mult_bound_push', 1e-5),
					   ('halt_on_ampl_error', 'yes')],
    file_name='EDC_asNMPC_MPM_test', real_plant = plant, dynamic_plot=True)

## ================================================================================================
## END
