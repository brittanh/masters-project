#!/opt/local/bin/python
# -*- encoding: ascii -*-
"""
    @purpose: Plots the results
    @author: Brittany Hall
    @date: 08.11.2017
    @version: 0.1
    @updates:
"""
import matplotlib.pyplot as plt
import scipy.io as spio
from numpy import reshape, append, hstack, linspace


def plotting(MPCit,T):
    #Loading in steady state data
    data = spio.loadmat('CstrDistXinit.mat', squeeze_me = True)
    Xinit = data['Xinit']
    xf = Xinit[0:84]
    u_opt = Xinit[84:]
    
    #Loading in iNMPC results
    #data_ideal = spio.loadmat('iNMPC.mat', squeeze_me = True)
    #uAll = data_ideal[iNMPC]['uAll']
    #xmeasureAll = data_ideal['iNMPC']['xmeasureAll']
    #ObjReg = data_ideal['iNMPC']['ObjReg']
    #ObjEcon = data_ideal['iNMPC']['ObjEcon']
    
    #Loading in pfNMPC results
    #data_pf = spio.loadmat('pfNMPC.mat', squeeze_me = True)
    #uAll_pf = data_pf[pfNMPC]['uAll']
    #xmeasureAll_pf = data_pf['pfNMPC']['xmeasureAll']
    #ObjReg_pf = data_pf['pfNMPC']['ObjReg']
    #ObjEcon_pf = data_pf['pfNMPC']['ObjEcon']
    

    nu = u0.shape[0]
    uAll = reshape(uAll, (nu, MPCit))
    uAll_pf = reshape(uAll_pf,(nu,MPCit))

    #Add initial control
    uAll = append(u0[:,0] uAll )
    uAll_pf = append(u0[:,0], uAll_pf)

    #Add initial states
    xmeasureAll = hstack(xmeasure,xmeasureAll)
    xmeasureAll_pf = hstack(xmeasure,xmeasureAll_pf)

    x = linspace(1,MPCit, MPCit/T)
    xi = append(0,x)

    #Plotting
    plt.plot(x,)



