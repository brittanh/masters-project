"""
================================================================================================
                 Steady state model for the extractive distillation column
                         Federico Lozano Santamaria - June 2015                           
================================================================================================
"""

## ================================================================================================
## Imports and model creation
## ================================================================================================
from __future__ import division
from pyomo.environ import *
ms = AbstractModel()

## ================================================================================================
## Sets
## ================================================================================================
ms.comp = Set()
ms.NT = Param(within=PositiveIntegers) #Total number of stages
ms.Net = RangeSet(1,ms.NT)

## ================================================================================================
## Parameters
## ================================================================================================
ms.FA = Param(ms.Net, within=NonNegativeReals, initialize=0) # Azeotropic mix feed flow rate [kmol/h]
ms.za = Param(ms.comp, within=NonNegativeReals) # Azeotropic mixture composition
ms.TfeedA = Param(within=PositiveReals) # Azeotropic mixture feed temperature

ms.FG = Param(ms.Net, within=NonNegativeReals, initialize=0) # Glycerol feed flow rate [kmol/h]
ms.zg = Param(ms.comp, within=NonNegativeReals) # Glycerol feed composition
ms.TfeedG = Param(within=PositiveReals) # Glycerol feed temperature

ms.xD_Ethanol = Param(within=PositiveReals) # Ethanol composition in distillate

ms.P = Param(ms.Net, within=NonNegativeReals) # Preassure [bar]

ms.ML_C = Param(within=NonNegativeReals)
ms.ML_R = Param(within=NonNegativeReals)

# Tray and column geometry
ms.CD = Param(within=PositiveReals)
ms.Lw = Param(within=PositiveReals)
ms.hw = Param(within=PositiveReals)
ms.AC = Param(within=PositiveReals)
ms.AT = Param(within=PositiveReals)
ms.A0 = Param(within=PositiveReals)
ms.da = Param(within=PositiveReals)
ms.ep = Param(within=PositiveReals)
ms.pitch = Param(within=PositiveReals)
ms.HS = Param(within=PositiveReals)

## ================================================================================================
## Variables
## ================================================================================================
ms.ML = Var(ms.Net, bounds=(0,20), initialize=0.2) # [kmol]
ms.L = Var(ms.Net, bounds=(0,5000), initialize=200) # [kmol/h]
ms.V = Var(ms.Net, bounds=(0,5000), initialize=100) # [kmol/h]
ms.D = Var(bounds=(0,5000), initialize=80) # [kmol/h]
ms.x = Var(ms.Net, ms.comp, bounds=(0,100), initialize=30) # 0-100
ms.y = Var(ms.Net, ms.comp, bounds=(0,100), initialize=30) # 0-100
ms.Temp = Var(ms.Net, bounds=(273,500), initialize=330) # [K]
ms.Qc = Var(bounds=(0,1E8), initialize=5000) # [MJ/h]
ms.Qr = Var(bounds=(0,1E8), initialize=5000) # [MJ/h]

## ================================================================================================
## Purity constraint
## ================================================================================================
def _purity(ms,i,j):
  if i == 1 and j == 'Ethanol':
    return ms.x[i,j] >= ms.xD_Ethanol
  else:
    return Constraint.Skip
ms.Purity = Constraint(ms.Net, ms.comp, rule=_purity)

## ================================================================================================
## MESH equations
## ================================================================================================

# -------------------------------------------------------------------------------------------------
# Thermodynamic equilibrium
# -------------------------------------------------------------------------------------------------
# Constants for saturation preassure
ms.C1sp = Param(ms.comp)
ms.C2sp = Param(ms.comp)
ms.C3sp = Param(ms.comp)
ms.C4sp = Param(ms.comp)
ms.C5sp = Param(ms.comp)
ms.C6sp = Param(ms.comp)
ms.C7sp = Param(ms.comp)
# Constants for NRTL
ms.a = Param(ms.comp,ms.comp)
ms.b = Param(ms.comp,ms.comp)
ms.c = Param(ms.comp,ms.comp)

# Equation
# Saturation preassure
def _EqPsat(ms,i,j):
  return exp( ms.C1sp[j] + (ms.C2sp[j]/( ms.Temp[i]+ms.C3sp[j])) + (ms.C4sp[j]* ms.Temp[i]) + (ms.C5sp[j]*log( ms.Temp[i])) + (ms.C6sp[j]*( ms.Temp[i]**ms.C7sp[j])) ) 
ms.Psat = Expression(ms.Net, ms.comp, rule=_EqPsat)

# Activity coefficient
def _EqTau(ms,k,i,j):
  return ms.a[i,j] + (ms.b[i,j]/ms.Temp[k])
ms.Tau = Expression(ms.Net, ms.comp, ms.comp, rule=_EqTau )

def _EqG(ms,k,i,j):
  return exp(-ms.c[i,j]*ms.Tau[k,i,j])
ms.G = Expression(ms.Net, ms.comp, ms.comp, rule=_EqG )

def _EqS(ms,k,i):
  return sum( (ms.x[k,j]/100)*ms.G[k,j,i] for j in ms.comp )
ms.S = Expression(ms.Net, ms.comp, rule=_EqS)

def _EqC(ms,k,i):
  return sum( (ms.x[k,j]/100)*ms.G[k,j,i]*ms.Tau[k,j,i] for j in ms.comp )
ms.C = Expression(ms.Net, ms.comp, rule=_EqC)

def _EqEps(ms,j,i,k):
  return (ms.G[j,i,k]/ms.S[j,k])*(ms.Tau[j,i,k] - (ms.C[j,k]/ms.S[j,k]))
ms.epsilon = Expression(ms.Net, ms.comp, ms.comp, rule=_EqEps )

def _Eqgamma(ms,j,i):
  return exp( (ms.C[j,i]/ms.S[j,i]) + sum( (ms.x[j,k]/100)*ms.epsilon[j,i,k] for k in ms.comp ) )
ms.gamma = Expression(ms.Net, ms.comp, rule=_Eqgamma )

# Equation
def _EK(ms,i,j):  
  return ms.y[i,j]*ms.P[i] == ms.Psat[i,j]*ms.gamma[i,j]*ms.x[i,j]
ms.ThermodynamicEquilibrium = Constraint(ms.Net, ms.comp, rule=_EK)


# -------------------------------------------------------------------------------------------------
# Phase equilibrium error
# -------------------------------------------------------------------------------------------------
def _EE(ms,i):
  return sum( ms.y[i,j] - ms.x[i,j] for j in ms.comp ) == 0
ms.EquilibriumError = Constraint(ms.Net, rule=_EE)

# -------------------------------------------------------------------------------------------------
# Total mass balance
# -------------------------------------------------------------------------------------------------
def _TMB(ms,j):
  if j == 1: # Condenser
    return ms.V[j+1] == ms.L[j] + ms.D
  elif j == ms.NT: # Reboiler
    return ms.L[j-1] == ms.L[j] + ms.V[j]
  else: # Column
    return ms.FA[j] + ms.FG[j] + ms.V[j+1] + ms.L[j-1] == ms.V[j] + ms.L[j]
ms.TotalMassBalance = Constraint(ms.Net, rule=_TMB)

def _V0(ms,j):
  if j == 1:
    return ms.V[j] == 0
  else:
    return Constraint.Skip
ms.VaporCondenser = Constraint(ms.Net, rule=_V0)

# -------------------------------------------------------------------------------------------------
# Partial mass balance
# -------------------------------------------------------------------------------------------------
def _PMB(ms,i,j):
  
  if i == 1: # Condenser
    return ms.y[i+1,j] == ms.x[i,j]
  elif i == ms.NT: # Reboiler
    return ms.x[i-1,j]*ms.L[i-1] == ms.x[i,j]*ms.L[i] + ms.y[i,j]*ms.V[i]
  else: # Column
    return ms.za[j]*ms.FA[i] + ms.zg[j]*ms.FG[i] + ms.y[i+1,j]*ms.V[i+1] + ms.x[i-1,j]*ms.L[i-1] == ms.y[i,j]*ms.V[i] + ms.x[i,j]*ms.L[i]
ms.PartialMassBalance = Constraint (ms.Net, ms.comp, rule=_PMB)

# -------------------------------------------------------------------------------------------------
# Liquid hold-up
# -------------------------------------------------------------------------------------------------

# Constants for liquid density
ms.MW = Param(ms.comp)
ms.C1r = Param(ms.comp)
ms.C2r = Param(ms.comp)
ms.C3r = Param(ms.comp)
ms.C4r = Param(ms.comp)
ms.C5r = Param(ms.comp)

# Liquid density equation
def _EqRhoL(ms,i):
  rho = {}
  for j in ms.comp:
    if j == 'Water':
      Tr = 1 - ((ms.Temp[i]-273.15)/373.946)
      rho[j] = ms.C1r[j] + (ms.C2r[j]*(Tr**0.35)) + (ms.C3r[j]*(Tr**(2.0/3.0))) + (ms.C4r[j]*Tr) + (ms.C5r[j]*(Tr**(4.0/3.0))) 
    else:
      rho[j] = ms.C1r[j] / ( ms.C2r[j]**( 1+ ((1-(ms.Temp[i]/ms.C3r[j]))**ms.C4r[j]) ) )
  rhoL_mol = sum( (ms.x[i,j]/100)*rho[j] for j in ms.comp)    
  return rhoL_mol
ms.rhoL = Expression(ms.Net, rule=_EqRhoL)

# Vapor density equation
def _EqRhoV(ms,i):
  return ms.P[i] / (0.08314*ms.Temp[i]) # R = 0.08314 m3*bar/K*kmol  
ms.rhoV = Expression(ms.Net, rule=_EqRhoV)

# Flow factpr
def _EqFP(ms,i):
  if i == 1:
    return 0
  elif i == ms.NT:
    return 0
  else:
    return ((ms.L[i]/ms.rhoL[i])/(ms.V[i]/ms.rhoV[i]))*(((ms.rhoL[i]*sum((ms.x[i,j]/100)*ms.MW[j] for j in ms.comp))/(ms.rhoV[i]*sum((ms.y[i,j]/100)*ms.MW[j] for j in ms.comp)))**(0.5)) 
ms.FP = Expression(ms.Net, rule=_EqFP)

# Liquid flow rate
def _LF(ms,i):
  if i == 1: # Condenser (steady-state)
    return ms.ML[i] == ms.ML_C
  elif i == ms.NT: # Reboiler (steady-state)
    return ms.ML[i] == ms.ML_R
  else: # Column    
    return ms.ML[i] == 0.6*ms.rhoL[i]*ms.AT*(ms.hw**0.5)*((ms.FP[i]*ms.AT*ms.pitch/ms.Lw)**0.25)
ms.LiquidFlowRate = Constraint(ms.Net, rule=_LF)
  
# -------------------------------------------------------------------------------------------------
# Energy balances
# -------------------------------------------------------------------------------------------------

# Constants for vapor enthalpy
ms.H0 = Param(ms.comp)
ms.C1c = Param(ms.comp)
ms.C2c = Param(ms.comp)
ms.C3c = Param(ms.comp)
ms.C4c = Param(ms.comp)
ms.C5c = Param(ms.comp)
ms.C6c = Param(ms.comp)
# Constants for liquid enthalpy
ms.Tcrit = Param(ms.comp)
ms.C1v = Param(ms.comp)
ms.C2v = Param(ms.comp)
ms.C3v = Param(ms.comp)
ms.C4v = Param(ms.comp)
ms.C5v = Param(ms.comp)
# Equation

# Vapor enthalpy 
Tref = 298.15

def _EqHVi(ms,i,j):
  return (ms.C1c[j]*(ms.Temp[i]-Tref)) + ((ms.C2c[j]/2)*((ms.Temp[i]**2)-(Tref**2))) + ((ms.C3c[j]/3)*((ms.Temp[i]**3)-(Tref**3))) + ((ms.C4c[j]/4)*((ms.Temp[i]**4)-(Tref**4))) + ((ms.C5c[j]/5)*((ms.Temp[i]**5)-(Tref**5))) + ((ms.C6c[j]/6)*((ms.Temp[i]**6)-(Tref**6))) + ms.H0[j]
ms.HVi = Expression(ms.Net, ms.comp, rule=_EqHVi) # [MJ/kmol]

def _EqHV(ms,i):
  return sum(ms.HVi[i,j]*ms.y[i,j]/100 for j in ms.comp)
ms.HV = Expression(ms.Net, rule=_EqHV ) # [MJ/kmol]

# Liquid enthalpy
def _EqDHvapi(ms,i,j):
  Tr = ms.Temp[i]/ms.Tcrit[j]
  return ms.C1v[j]*( (1-Tr)**(ms.C2v[j] + (ms.C3v[j]*Tr) + (ms.C4v[j]*(Tr**2)) + (ms.C5v[j]*(Tr**3)) ) )
ms.DHvapi = Expression(ms.Net, ms.comp, rule=_EqDHvapi ) # [MJ/kmol]

def _EqHLi(ms,i,j):
  return ms.HVi[i,j] - ms.DHvapi[i,j]
ms.HLi = Expression(ms.Net, ms.comp, rule=_EqHLi ) # [MJ/kmol]

def _EqHL(ms,i):
  return sum( ms.HLi[i,j]*ms.x[i,j]/100 for j in ms.comp)
ms.HL = Expression(ms.Net, rule=_EqHL ) # [MJ/kmol]

# Equation
from Enthalpies import Eq_HL
def _EB(ms,j):
  za_t = {}
  if j == 1: # Condenser
    return ms.V[j+1]*ms.HV[j+1] - (ms.L[j] + ms.D)*ms.HL[j] - ms.Qc == 0
  elif j == ms.NT: # Reboiler
    return ms.L[j-1]*ms.HL[j-1] + ms.Qr - ms.V[j]*ms.HV[j] - ms.L[j]*ms.HL[j] == 0
  else: # Column
    for k in ms.comp:
      za_t[k] = ms.za[k]
    HFA = Eq_HL(ms.comp, za_t, ms.TfeedA, ms.H0, ms.C1c, ms.C2c, ms.C3c, ms.C4c, ms.C5c, ms.C6c, ms.Tcrit, ms.C1v, ms.C2v, ms.C3v, ms.C4v, ms.C5v)
    HFG = Eq_HL(ms.comp, ms.zg, ms.TfeedG, ms.H0, ms.C1c, ms.C2c, ms.C3c, ms.C4c, ms.C5c, ms.C6c, ms.Tcrit, ms.C1v, ms.C2v, ms.C3v, ms.C4v, ms.C5v)      
    return ms.FA[j]*HFA + ms.FG[j]*HFG + ms.V[j+1]*ms.HV[j+1] + ms.L[j-1]*ms.HL[j-1] - ms.V[j]*ms.HV[j] - ms.L[j]*ms.HL[j] == 0  
ms.EnergyBalance = Constraint(ms.Net, rule=_EB)


## ================================================================================================
## Objective
## ================================================================================================
def _obj(ms):
  return 1E-3*ms.Qr - 30*ms.D
ms.ObjectiveFunction = Objective(rule=_obj, sense=minimize)

## ================================================================================================
## END

