#!/opt/local/bin/python
# -*- encoding: ascii -*-
"""
    @purpose: Used to reshape the data to make it easier for plotting
    @author: Brittany Hall
    @date: 08.10.2017
    @version: 0.1
    @updates:
"""
from numpy import array, zeros, reshape, delete, size
from casadi import *

def plotStates(data, lb, ub, N, params):
    #unpacking params
    nu = params['prob']['nu']
    nx = params['prob']['nx']
    ns = params['prob']['ns']
    nk = params['prob']['nk']
    d = params['prob']['d']
    
    #Optimized initial state
    x0_opt = data[0:nx]
    index = range(0,nx)
    data = delete(data,index)
    data = reshape(data, ((nu + (nx+ns)*d + (nx+ns)),N*nk))
    u_nlp_opt = data[0:nu,0:N*nk]
    data = data[nu:,:]
    
    lb0 = lb[0:nx+ns]
    lb = delete(lb,range(0,nx))
    #print where(lb!=0)[0]
    lb = reshape(lb,(nu+(nx+ns)*d+(nx+ns),N*nk))
    lbU = lb[0:nu, 0:N*nk]
    lb = lb[nu:,:]
    ub0 = ub[0:nx+ns]
    ub = ub[nx:]
    ub = reshape(ub, (nu+(nx+ns)*d+(nx+ns),N*nk))
    ubU = ub[0:nu,0:N*nk]
    ub = ub[nu:,:]

    #Preparing matrix for plotting
    nState = (nx+ns) + N*nk*(d+1)*(nx+ns)
    nPoint = nState/(nx+ns)
    plotState = zeros((nx+ns,nPoint))
    for i in range(0,nx):
        plotState[i,0] = x0_opt[i]
    plotLb = zeros((nx+ns, nPoint))
    plotLb[:,0] = lb0
    plotUb = zeros((nx+ns,nPoint))
    plotUb[:,0] = ub0

    #Extract states from each collocation point at each time horizon
    sInd = 1 #initial index row
    for i in range(0,N*nk-1):
        temp = data[:,i]
        numCol = size(temp,axis=0)
        numRow = numCol/(nx+ns)
        temp = reshape(temp,(nx+ns,numRow))
        plotState[:,sInd:(numRow+sInd)] = temp
        tempLb = lb[:,i]
        tempLb = reshape(tempLb, (nx+ns,numRow))
        plotLb[:,sInd:(numRow+sInd)] = tempLb
        tempUb = ub[:,i]
        tempUb = reshape(tempUb, (nx+ns,numRow))
        plotUb[:,sInd:(numRow+sInd)] = tempUb

        sInd += numRow

    return u_nlp_opt, plotState
