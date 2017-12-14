#!/opt/local/bin/python
# -*- encoding: ascii -*-
"""
    @purpose: Generates noise for states
    @author: Brittany Hall
    @date: 23.10.2017
    @version: 0.1
    @updates:
"""
import scipy.io as spio
from numpy import zeros, array, append, random

mpc_iter = 500
noiselevel = 0.1 # 1 percent noise

#Load in steady state data
data = spio.loadmat('CstrDistXinit.mat', squeeze_me = True)
Xinit = data['Xinit']
xf = Xinit[0:2*NT+2]
xholdup = xf[NT+1:-1]

noise = array([])
for i in range(0,mpc_iter):
    noise = append(noise, noiselevel*xholdup*random.randn(NT+1,1))

print noise
raw_input()
spio.savemat('noise1pct',noise)
