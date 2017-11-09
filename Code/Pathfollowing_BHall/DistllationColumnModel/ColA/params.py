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
MO[1:NT-1] = 0.5                            #Nominal stage (tray) holdup [kmol]
MO[NT-1] = 0.5                                 #Nominal condenser holdup [kmol]
MO[NT] = 0.5                                        #Nominal CSTR hold up [kmol]
#Linearized flow dynamics (NA to reboiler and condenser)
taul = 0.063                           #Time constant for liquid dynamics [min]
L0 = 2.70629
L0b = L0 + qF*F0                     #Nominal liquid flow below feed [kmol/min]
lam = 0
V0 = 3.206
VB_max = 4.008
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
#Nominal holdup values
MDs = 0.5
MBs = 0.5
#Nominal flow rates
Ds = 0.5
Bs = 0.5
#Constraint bounds
u_min = array([[0.1], [0.1], [0.1], [0.1], [0.1]])
u_max = array([[10],[VB_max],[10],[1.0],[1.0]])
#State bounds
x_min = zeros((2*NT+2,1))
x_max = ones((2*NT+2,1))

lbx = concatenate((x_min, u_min))
ubx = concatenate((x_max, u_max))
lbg = zeros((2*NT+2,1))
ubg = zeros((2*NT+2,1))
#Problem Dimensions
nx = 2*NT+2                      #Number of states (CSTR + Distillation Column)
nu = 5                                      #Number of inputs (LT, VB, F, D, B)
nk = 1
tf = 1
h = tf/nk
ns = 0
#Collecting all parameters into a dictionary
params = {}
params['dist'] = {'NC':NC,'F_0': F_0, 'NT': NT, 'zF': zF, 'qF': qF, 'NF': NF,
            'VB': VB, 'LT': LT, 'F': F, 'alpha': alpha, 'B': B, 'D': D, 'zF': zF,
            'Muw': Muw, 'L0': L0, 'L0b': L0b, 'qF0': qF0, 'F0': F0, 'taul': taul,
                'V0':V0, 'lam':lam, 'MO': MO}
params['cstr'] = {'k1': k1}
params['price'] = {'pf': pf, 'pV': pV, 'pB': pB, 'pD': pD}
params['bounds'] = {'x_min':x_min, 'x_max':x_max, 'u_min': u_min, 'u_max': u_max, 'lbx': lbx, 'ubx': ubx, 'ubg': ubg, 'lbg': lbg}
params['gain'] = {'MDs':MDs,'MBs':MBs,'Ds':Ds,'Bs':Bs, 'KcD':KcD, 'KcB':KcB}
params['prob'] = {'nx':nx, 'nu':nu, 'nk':nk, 'tf': tf, 'h': h, 'ns':ns}

