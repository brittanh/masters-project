#!/opt/local/bin/python
# -*- encoding: ascii -*-
"""
    @purpose: Used to reshape the data to make it easier for plotting
    @author: Brittany Hall
    @date: 08.10.2017
    @version: 0.1
    @updates:
"""
from numpy import array, zeros, reshape
from casadi import *

def plotStates(data, lb, ub, N, params):
    #unpacking params
    nu = params['prob']['nu']
    nx = params['prob']['nx']
    ns = params['prob']['ns']
    nk = params['prob']['nk']
    
    #Optimized initial state
    x0_opt = data[0:nx]
    data = data[nx-1:-1]
    print data.shape
    print data
    raw_input()
    data = reshape(data, ((nu + (nx+ns)*d + (nx+ns)),N*nk))
    u_nlp_opt = data[0:nu,0:N*nk]
    data = data[nu-1:-1]

    lb0 = lb[0:nx+ns]
    lb[0:nx] = array([])
    lb = reshape(lb,(nu+(nx+ns)*d+(nx+ns),N*nk))
    lbU = lb[0:nu, 0:N*nk]
    lb[0:nu,:] = array([])
    ub0 = ub[0:nx+ns]
    ub[0:nx] = array([])
    ub = reshape(ub, (nu+(nx+ns)*d+(nx+ns),N*nk))
    ubU = ub[0:nu,0:N*nk]
    ub[0:nu,:] = []

    #Preparing matrix for plotting
    nState = (nx+ns) + N*nk*(d+1)*(nx+ns)
    nPoint = nState/(nx+ns)
    plotState = zeros((nx+ns,nPoint))
    plotState[0:nx,0] = x0_opt

    plotLb = zeros((nx+ns, nPoint))
    plotLb[:,0] = lb0
    plotUb = zeros((nx+ns,nPoint))
    plotUb[:,0] = ub0

    #Extract states from each collocation point at each time horizon
    sInd = 2 #initial index row
    for i in range(0,N*nk):
        temp = data[:,i]
        numCol = shape(temp,0)
        numRow = numCol/(nx+ns)
        temp = reshape(temp,(nx+ns,numRow))
        plotState[:,sInd:(numRow+sInd-1)] = temp

        tempLb = lb[:,i]
        tempLb = reshape(tempLb, (nx+ns,numRow))
        plotLb[:,sInd:(numRow+sInd-1)] = tempLb
        tempUb = ub[:,i]
        tempUb = reshape(tempUb, (nx+ns,numRow))
        plotUb[:,sInd:(numRow+sInd-1)] = tempUb

        sInd += numRow

    return u_nlp_opt, plotState
