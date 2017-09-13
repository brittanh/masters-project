## Imports
from __future__ import division
from pyomo.environ import *

def Eq_HV(comp, y, T, H0, C1, C2, C3, C4, C5, C6):
  """
  Function that calculates the enthalpy of a vapor mixture
    *Inputs:
      comp: set of components defined as the keys of a dictionary
      y: molar fraction of the vapor mixture (0 - 100). Dictionary indexed by comp
      T: temperature [K]
      H0: enthalpy of formation [kJ/kmol]
      C1..C6: constant of the equation (integration of Cp polynomial)
    *Output:
      HV: enthalpy of vapor mixture [kJ/kmol]
  The equation is:
  Cp(comp) = C1 + C2*T + C3*T**2 + C4*T**3 + C5*T**4 + C6*T**5
  HV(comp) = int(Cp,T)
  HV = sum( y*HV )
  """
  
  Tref = 298.15

  HVi = {}
  for i in comp:
    HVi[i] = (C1[i]*(T-Tref)) + ((C2[i]/2)*((T**2)-(Tref**2))) + ((C3[i]/3)*((T**3)-(Tref**3))) \
             + ((C4[i]/4)*((T**4)-(Tref**4))) + ((C5[i]/5)*((T**5)-(Tref**5))) \
             + ((C6[i]/6)*((T**6)-(Tref**6))) + H0[i]
  
  HV = sum( HVi[i]*y[i]/100 for i in comp)
  
  return HV

def Eq_HL(comp, x, T, H0, C1c, C2c, C3c, C4c, C5c, C6c, Tcrit, C1v, C2v, C3v, C4v, C5v):
  """
  Function that calculates the enthalpy of a liquid mixture
    *Inputs:
      comp: set of components defined as the keys of a dictionary
      x: molar fraction of the liquid mixture (0 - 100). Dictionary indexed by comp
      T: temperature [K]
      H0: enthalpy of formation [kJ/kmol]
      C1c..C6c: constant of the equation (integration of Cp polynomial)
      Tcrit: critical temperature [K]
      C1v..C5v: constant of the equation (enthalpy of vaporization)
    *Output:
      HL: enthalpy of liquid mixture [kJ/kmol]
  The equation is:
  Cp(comp) = C1c + C2c*T + C3c*T**2 + C4c*T**3 + C5c*T**4 + C6c*T**5
  HV(comp) = int(Cp,T)
  DHvap(comp) = C1v*( (1-Tr)**(C2v + C3v*Tr + C4v*Tr**2 + C5v*Tr**3 ) )
  HL(comp) = HV - DHvap
  HL = sum( x*HL)
  """
  
  Tref = 298.15

  HVi = {}
  HLi = {}
  DHvap = {}
  Tr = {}
  for i in comp:
    Tr[i] = T/Tcrit[i]            
    HVi[i] = (C1c[i]*(T-Tref)) + ((C2c[i]/2)*((T**2)-(Tref**2))) + ((C3c[i]/3)*((T**3)-(Tref**3)))\
             + ((C4c[i]/4)*((T**4)-(Tref**4))) + ((C5c[i]/5)*((T**5)-(Tref**5))) \
             + ((C6c[i]/6)*((T**6)-(Tref**6))) + H0[i]
    DHvap[i] = C1v[i]*( (1-Tr[i])**(C2v[i] + (C3v[i]*Tr[i]) + (C4v[i]*(Tr[i]**2)) + (C5v[i]*(Tr[i]**3)) ) )
    HLi[i] = HVi[i] - DHvap[i]

  HL = sum( HLi[i]*x[i]/100 for i in comp)

  return HL

## Test
'''
m = AbstractModel()

m.comp = Set()

m.H0 = Param(m.comp)

m.C1c = Param(m.comp)
m.C2c = Param(m.comp)
m.C3c = Param(m.comp)
m.C4c = Param(m.comp)
m.C5c = Param(m.comp)
m.C6c = Param(m.comp)

m.Tcrit = Param(m.comp)

m.C1v = Param(m.comp)
m.C2v = Param(m.comp)
m.C3v = Param(m.comp)
m.C4v = Param(m.comp)
m.C5v = Param(m.comp)

m.y = Param(m.comp)
m.x = Param(m.comp)

m1 = m.create('Enthalpies.dat')

print Eq_HV(m1.comp, m1.y, 351.111, m1.H0, m1.C1c, m1.C2c, m1.C3c, m1.C4c, m1.C5c, m1.C6c)

print Eq_HL(m1.comp, m1.x, 351.111, m1.H0, m1.C1c, m1.C2c, m1.C3c, m1.C4c, m1.C5c, m1.C6c, m1.Tcrit, m1.C1v, m1.C2v, m1.C3v, m1.C4v, m1.C5v)
'''
