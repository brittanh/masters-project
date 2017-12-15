#!/opt/local/bin/python
# -*- encoding: ascii -*-
"""
    @purpose: Plots the results (iNMPC vs pfNMPC)
    @author: Brittany Hall
    @date: 08.11.2017
    @version: 0.1
    @updates:
"""
import matplotlib.pyplot as plt
import scipy.io as spio
from numpy import reshape, append, hstack, linspace, ones, transpose, vstack
from params import params

def plotting(u0, xmeasure, MPCit, T):

    NT = params['dist']['NT']
    NF = params['dist']['NF']
    
    #Loading in steady state data
    data = spio.loadmat('CstrDistXinit.mat', squeeze_me = True, struct_as_record=False)
    Xinit = data['Xinit']
    xf = Xinit[0:84]
    u_opt = Xinit[84:]
    
    #Loading in iNMPC results
    data_ideal = spio.loadmat('iNMPC.mat', squeeze_me = False)
    uAll = data_ideal['ideal']['uAll']
    uAll = uAll[0,0]
    xmeasureAll = data_ideal['ideal']['xmeasureAll']
    xmeasureAll = xmeasureAll[0,0]
    ObjReg = data_ideal['ideal']['ObjReg']
    ObjReg = transpose(ObjReg[0,0])
    ObjEcon = data_ideal['ideal']['ObjEcon']
    ObjEcon = transpose(ObjEcon[0,0])
    T = data_ideal['ideal']['T']
    mpcit = data_ideal['ideal']['mpciterations']
    
    #Loading in iNMPC MATLAB results
    data_iMat = spio.loadmat('iNmpcData.mat', squeeze_me = False)
    xmeasureAll_mat = data_iMat['xmeasureAll']
    uAll_mat = data_iMat['uAll']
    ObjReg_mat = transpose(data_iMat['ObjReg'])
    ObjEcon_mat = transpose(data_iMat['ObjEcon'])

    #Loading in pfNMPC results
    #data_pf = spio.loadmat('pfNMPC.mat', squeeze_me = True)
    #uAll_pf = data_pf['pfNMPC']['uAll']
    #xmeasureAll_pf = data_pf['pfNMPC']['xmeasureAll']
    #ObjReg_pf = data_pf['pfNMPC']['ObjReg']
    #ObjEcon_pf = data_pf['pfNMPC']['ObjEcon']
    
    nu = u0.shape[0]
    uAll = uAll.reshape(nu, MPCit,order='F').copy()
    uAll_mat = uAll_mat.reshape(nu, MPCit,order='F').copy()
    #uAll_pf = reshape(uAll_pf,(nu,MPCit))

    #Add initial control
    u0_0 = reshape(u0[:,0],(nu,1))
    uAll = hstack((u0_0, uAll))
    uAll_mat = hstack((u0_0, uAll_mat))
    #uAll_pf = append(u0[:,0], uAll_pf)

    #Add initial states
    xmeasure = reshape(xmeasure,(xmeasure.shape[0],1))
    xmeasureAll = hstack((xmeasure,xmeasureAll))
    #xmeasureAll_pf = hstack(xmeasure,xmeasureAll_pf)

    x = linspace(0, MPCit, MPCit/T)
    xi = append(0,x)

    #Plotting
    #Figure: Objective function comparison
    plt.plot(x,ObjReg, 'g', x, ObjEcon, 'b', x, ObjReg_mat, 'ro', x, ObjEcon_mat, 'k*')
    plt.title('Objective function')
    plt.xlabel('Number of MPC iteration [-]')
    plt.ylabel('Objective function [-]')
    plt.legend(['iNMPC:Full-Python', 'iNMPC:Economic-Python', 'iNMPC:Full-Matlab', 'iNMPC:Economic-Matlab'])
    plt.show()
    
    #Figure: Concentration at stage 1 (reboiler)
    plt.plot(xi,xf[0]*ones(MPCit+1),'r', xi, xmeasureAll[0,],'g', xi[0:150], xmeasureAll_mat[0,],'bo')
    plt.ylabel('Concentration [-]')
    plt.xlabel('Time [min]')
    plt.title('Distillation: Bottom Composition')
    plt.legend(['Steady-state','iNMPC-Python', 'iNMPC-Matlab'])
    plt.show()
    
    #Figure: Concentration at feed stage
    plt.plot(xi,xf[NF]*ones(MPCit+1),'r',xi,xmeasureAll[NF,],'g', xi[0:150], xmeasureAll_mat[NF,],'bo')
    plt.ylabel('Concentration [-]')
    plt.xlabel('Time [min]')
    plt.title('Distillation: Feed Composition')
    plt.legend(['Steady-state','iNMPC-Python','iNMPC-Matlab'])
    plt.show()
    
    #Figure: Concentration at stage NT (top)
    plt.plot(xi,xf[NT]*ones(MPCit+1),'r',xi,xmeasureAll[NT,],'g',xi[0:150],xmeasureAll_mat[NT,],'bo')
    plt.ylabel('Concentration [-]')
    plt.xlabel('Time [min]')
    plt.title('Distillation: Top Composition')
    plt.legend(['Steady-state','iNMPC-Python', 'iNMPC-Matlab'])
    plt.show()

    #Figure: Concentration in CSTR
    plt.plot(xi,xf[NT+1]*ones(MPCit+1),'r',xi,xmeasureAll[NT+1,],'g', xi[0:150], xmeasureAll_mat[NT+1,],'bo')
    plt.ylabel('Concentration [-]')
    plt.xlabel('Time [min]')
    plt.title('CSTR: Concentration')
    plt.legend(['Steady-state','iNMPC-Python', 'iNMPC-Matlab'])
    plt.show()
              
    #Figure: Holdup in CSTR
    plt.plot(xi,xf[2*NT-1]*ones(MPCit+1),'r',xi,xmeasureAll[2*NT-1,],'g', xi[0:150], xmeasureAll_mat[2*NT-1,],'bo')
    plt.ylabel('Holdup [-]')
    plt.xlabel('Time [min]')
    plt.title('CSTR: Holdup')
    plt.legend(['Steady-state','iNMPC-Python', 'iNMPC-Matlab'])
    plt.show()
    
    #Figure: u[0] LT control input
    plt.plot(xi, u_opt[0]*ones(MPCit+1), 'r', xi, uAll[0,],'g', xi, uAll_mat[0,],'bo')
    plt.ylabel('LT [m^3/min]')
    plt.xlabel('Time [min]')
    plt.title('Control input for LT')
    plt.legend(['Steady-state','iNMPC-Python', 'iNMPC-Matlab'])
    plt.show()

    #Figure: u[1] VB control input
    plt.plot(xi, u_opt[1]*ones(MPCit+1), 'r', xi, uAll[1,],'g',xi, uAll_mat[1,],'bo')
    plt.ylabel('VB [m^3/min]')
    plt.xlabel('Time [min]')
    plt.title('Control input for VB')
    plt.legend(['Steady-state', 'iNMPC-Python', 'iNMPC-Matlab'])
    plt.show()

    #Figure: u[2] F control input
    plt.plot(xi, u_opt[2]*ones(MPCit+1), 'r', xi, uAll[2,],'g', xi, uAll_mat[2,],'bo')
    plt.ylabel('F [kmol/min]')
    plt.xlabel('Time [min]')
    plt.title('Control input for F')
    plt.legend(['Steady-state', 'iNMPC-Python', 'iNMPC-Matlab'])
    plt.show()

    #Figure: u[3] D control input
    plt.plot(xi, u_opt[3]*ones(MPCit+1), 'r', xi, uAll[3,],'g',xi, uAll_mat[3,],'bo')
    plt.ylabel('D [kmol/min]')
    plt.xlabel('Time [min]')
    plt.title('Control input for D')
    plt.legend(['Steady-state', 'iNMPC-Python', 'iNMPC-Matlab'])
    plt.show()

    #Figure: u[4] B control input
    plt.plot(xi, u_opt[4]*ones(MPCit+1), 'r', xi, uAll[4,],'g', xi, uAll_mat[4,],'bo')
    plt.ylabel('B [kmol/min]')
    plt.xlabel('Time [min]')
    plt.title('Control input for B')
    plt.legend(['Steady-state', 'iNMPC-Python', 'iNMPC-Matlab'])
    plt.show()

