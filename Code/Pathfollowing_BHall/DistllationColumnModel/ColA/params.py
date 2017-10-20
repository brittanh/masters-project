#!/opt/local/bin/python
# -*- encoding: ascii -*-
"""
    @purpose: Distillation Column A and CSTR model parameters
    @author: Brittany Hall
    @date: 11.10.2017
    @version: 0.1
    @updates:
"""
from numpy import zeros, ones, concatenate, array
params = {}
#---------------------Distillation column parameters-------------------------#
NC = 2                                                    #Number of components
NT = 41                                                       #Number of stages
NF = 21                                                 #Location of feed stage
LT = 2.827                                                              #Reflux
VB = 3.454                                                              #Boilup
F = 1.0                                                               #Feedrate
zF = array([[1.0],[0.0]])                      #Feed composition (# components)
D = 0.5                                                    #Distillate flowrate
B = 0.5                                                       #Bottoms flowrate
qF = 1.0                                                  #Feed liquid fraction
F_0 = 0.3                                                       #CSTR Feed rate
F0 = F                                             #Nominal feed rate to column
qF0 = qF
alpha = 1.5                                                #Relative volatility
#Nominal liquid holdups
Muw = 0.5
MO = zeros(NT+1)
MO[0] = 0.5                                     #Nominal reboiler holdup [kmol]
MO[1:NT-2] = 0.5                            #Nominal stage (tray) holdup [kmol]
MO[NT-1] = 0.5                                 #Nominal condenser holdup [kmol]
MO[NT] =0.5                                        #Nominal CSTR hold up [kmol]
#Linearized flow dynamics (NA to reboiler and condenser)
taul = 0.063                           #Time constant for liquid dynamics [min]
L0 = 2.70629
L0b = L0 + qF*F0                     #Nominal liquid flow below feed [kmol/min]
lam = 0
V0 = 3.206
#-------------------------CSTR parameters------------------------------------#
#Reaction
k1 = 34.1/60.0
#----------------------Objective Function & Constraints-----------------------#
#Prices
pf = 1
pV = 0.02
pB = 2
pD = 0
#Gains
KcB = 10
KcD = 10
#Nominal values
MDs = 0.5
MBs = 0.5
#Nominal flows
Ds = 0.5
Bs = 0.5
#Constraint bounds
lbu = array([[0.1], [0.1], [0.1], [0.1], [0.1]])
ubu = array([[10],[4.008],[10],[1.0],[1.0]])
#State bounds
x_min = zeros((2*NT+2,1))
x_max = ones((2*NT+2,1))
xB_max = 0.1
x_max[0] = xB_max
x_min[2*NT+1] = 0.3
x_max[2*NT+1] = 0.7
lbx = concatenate((x_min, lbu))
ubx = concatenate((x_max, ubu))
lbg = zeros((2*NT+2,1))
ubg = zeros((2*NT+2,1))

#Collecting all parameters into a dictionary
params = {}
params['dist'] = {'NC':NC,'F_0': F_0, 'NT': NT, 'zF': zF, 'qF': qF, 'NF': NF,
            'VB': VB, 'LT': LT, 'F': F, 'alpha': alpha, 'B': B, 'D': D, 'zF': zF,
            'Muw': Muw, 'L0': L0, 'L0b': L0b, 'qF0': qF0, 'F0': F0, 'taul': taul,
                'V0':V0, 'lam':lam, 'MO': MO}
params['cstr'] = {'k1': k1}
params['price'] = {'pf': pf, 'pV': pV, 'pB': pV, 'pD': pD}
params['bounds'] = {'x_min':x_min, 'x_max':x_max, 'lbu': lbu, 'ubu': ubu, 'lbx': lbx, 'ubx': ubx, 'ubg': ubg, 'lbg': lbg}
params['gain']={'MDs':MDs,'MBs':MBs,'Ds':Ds,'Bs':Bs, 'KcD':KcD, 'KcB':KcB}


