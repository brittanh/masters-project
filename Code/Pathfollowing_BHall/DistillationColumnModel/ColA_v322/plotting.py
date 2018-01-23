#!/opt/local/bin/python
# -*- encoding: ascii -*-
"""
    @purpose: Plots the results (iNMPC vs pfNMPC and MATLAB vs Python)
    @author: Brittany Hall
    @date: 08.11.2017
    @version: 0.1
    @updates:
"""
import matplotlib.pyplot as plt
import scipy.io as spio
from numpy import reshape, append, hstack, linspace, ones, transpose, vstack, shape
from params import params

def plotting(u0, xmeasure, MPCit, T):

    NT = params['dist']['NT']
    NF = params['dist']['NF']
    
    ##=================Loading in Results==================================##
    
    #Steady state data
    data = spio.loadmat('CstrDistXinit.mat', squeeze_me = True, struct_as_record=False)
    Xinit = data['Xinit']
    xf = Xinit[0:84]
    u_opt = Xinit[84:]
    
    #Loading in iNMPC Python results
    data_ideal = spio.loadmat('iNMPC_Python.mat',squeeze_me = False)
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
    
    #Loading in pfNMPC Python results
    data_pf = spio.loadmat('pfNMPC_Python.mat',squeeze_me = False)
    uAll_pf = data_pf['pf']['uAll']
    uAll_pf = uAll_pf[0,0]
    xmeasureAll_pf = data_pf['pf']['xmeasureAll']
    xmeasureAll_pf = xmeasureAll_pf[0,0]
    ObjReg_pf = data_pf['pf']['ObjReg']
    ObjReg_pf = transpose(ObjReg_pf[0,0])
    ObjEcon_pf = data_pf['pf']['ObjEcon']
    ObjEcon_pf = transpose(ObjEcon_pf[0,0])
    T_pf = data_pf['pf']['T']
    mpcit_pf = data_pf['pf']['mpciterations']
    
    #Loading in iNMPC MATLAB results
    data_iMat = spio.loadmat('iNMPC_MATLAB.mat',squeeze_me = False)
    xmeasureAll_mat = data_iMat['xmeasureAll']
    uAll_mat = data_iMat['uAll']
    ObjReg_mat = transpose(data_iMat['ObjReg'])
    ObjEcon_mat = transpose(data_iMat['ObjEcon'])

    #Loading in pfNMPC MATLAB results
    data_pfmat = spio.loadmat('pfNMPC_MATLAB.mat', squeeze_me = True)
    uAll_pfmat = data_pfmat['uAll']
    xmeasureAll_pfmat = data_pfmat['xmeasureAll']
    ObjReg_pfmat = data_pfmat['ObjReg']
    ObjEcon_pfmat = data_pfmat['ObjEcon']

    #Reshaping values
    nu = u0.shape[0]
    uAll= uAll.reshape(nu, MPCit,order='F').copy()
    uAll_mat = uAll_mat.reshape(nu, MPCit,order='F').copy()
    uAll_pf = uAll_pf.reshape(nu, MPCit,order='F').copy()
    uAll_pfmat = uAll_pfmat.reshape(nu, MPCit,order='F').copy()

    #Add initial control
    u0_0 = reshape(u0[:,0],(nu,1))
    uAll = hstack((u0_0, uAll))
    uAll_mat = hstack((u0_0, uAll_mat))
    uAll_pf = hstack((u0_0, uAll_pf))
    uAll_pfmat = hstack((u0_0, uAll_pfmat))

    #Add initial states
    xmeasure = reshape(xmeasure,(xmeasure.shape[0],1))
    xmeasureAll = hstack((xmeasure,xmeasureAll))
    xmeasureAll_pf = hstack((xmeasure,xmeasureAll_pf))
    xmeasureAll_pfmat = hstack((xmeasure,xmeasureAll_pfmat))

    x = linspace(0, MPCit, MPCit/T)
    xi = append(0,x)

    ##====================Plotting (Python vs Matlab)======================##
    #Figure: Objective function comparison (Python vs Matlab)
    plt.plot(x,ObjReg, 'g', x, ObjEcon, 'b', x, ObjReg_mat, 'ro', x, ObjEcon_mat, 'y*', x, ObjReg_pf,'kh', x, ObjEcon_pf, 'c*', x, ObjReg_pfmat, 'm^', x, ObjEcon_pfmat, 'ys')
    plt.title('Objective functions: Python vs MATLAB')
    plt.xlabel('Number of MPC iteration [-]')
    plt.ylabel('Objective function [-]')
    plt.legend(['iNMPC:Full-Python', 'iNMPC:Economic-Python', 'iNMPC:Full-Matlab', 'iNMPC:Economic-Matlab', 'pfNMPC:Full-Python', 'pfNMPC:Economic-Python', 'pfNMPC:Full-Matlab', 'pfNMPC:Economic-Matlab'])
    plt.show()
    
    #Figure: Concentration at stage 1 (reboiler)
    plt.plot(xi,xf[0]*ones(MPCit+1),'r', xi, xmeasureAll[0,],'gs', xi[0:150], xmeasureAll_mat[0,],'bo', xi, xmeasureAll_pf[0,], 'k*', xi, xmeasureAll_pfmat[0,],'c^')
    plt.ylabel('Concentration [-]')
    plt.xlabel('Time [min]')
    plt.title('Distillation: Bottom Composition')
    plt.legend(['Steady-state','iNMPC-Python', 'iNMPC-Matlab', 'pfNMPC-Python', 'pfNMPC-Matlab'])
    plt.show()
    
    #Figure: Concentration at feed stage
    plt.plot(xi,xf[NF]*ones(MPCit+1),'r',xi,xmeasureAll[NF,],'gs', xi[0:150], xmeasureAll_mat[NF,],'bo', xi, xmeasureAll_pf[NF,],'k*', xi, xmeasureAll_pfmat[NF,],'c^')
    plt.ylabel('Concentration [-]')
    plt.xlabel('Time [min]')
    plt.title('Distillation: Feed Composition')
    plt.legend(['Steady-state','iNMPC-Python','iNMPC-Matlab', 'pfNMPC-Python', 'pfNMPC-Matlab'])
    plt.show()
    
    #Figure: Concentration at stage NT (top)
    plt.plot(xi,xf[NT]*ones(MPCit+1),'r',xi,xmeasureAll[NT,],'gs',xi[0:150],xmeasureAll_mat[NT,],'bo', xi, xmeasureAll_pf[NT,],'k*', xi, xmeasureAll_pfmat[NT,],'c^')
    plt.ylabel('Concentration [-]')
    plt.xlabel('Time [min]')
    plt.title('Distillation: Top Composition')
    plt.legend(['Steady-state','iNMPC-Python', 'iNMPC-Matlab', 'pfNMPC-Python', 'pfNMPC-Matlab'])
    plt.show()
    
    #Figure: Concentration in CSTR
    plt.plot(xi,xf[NT+1]*ones(MPCit+1),'r',xi,xmeasureAll[NT+1,],'gs', xi[0:150], xmeasureAll_mat[NT+1,],'bo', xi, xmeasureAll_pf[NT+1,],'k*', xi, xmeasureAll_pfmat[NT+1,],'c^')
    plt.ylabel('Concentration [-]')
    plt.xlabel('Time [min]')
    plt.title('CSTR: Concentration')
    plt.legend(['Steady-state','iNMPC-Python', 'iNMPC-Matlab', 'pfNMPC-Python', 'pfNMPC-Matlab'])
    plt.show()
    
    #Figure: Holdup in CSTR
    plt.plot(xi,xf[2*NT-1]*ones(MPCit+1),'r',xi,xmeasureAll[2*NT-1,],'gs', xi[0:150], xmeasureAll_mat[2*NT-1,],'bo', xi, xmeasureAll_pf[2*NT-1,], 'k*', xi, xmeasureAll_pfmat[2*NT-1,],'c^')
    plt.ylabel('Holdup [-]')
    plt.xlabel('Time [min]')
    plt.title('CSTR: Holdup')
    plt.legend(['Steady-state','iNMPC-Python', 'iNMPC-Matlab', 'pfNMPC-Python', 'pfNMPC-Matlab'])
    plt.show()
    
    #Figure: u[0] LT control input
    plt.plot(xi, u_opt[0]*ones(MPCit+1), 'r', xi, uAll[0,],'gs', xi, uAll_mat[0,],'bo', xi, uAll_pf[0,],'k*', xi, uAll_pfmat[0,],'c^')
    plt.ylabel('LT [m^3/min]')
    plt.xlabel('Time [min]')
    plt.title('Control input for LT')
    plt.legend(['Steady-state','iNMPC-Python', 'iNMPC-Matlab', 'pfNMPC-Python', 'pfNMPC-Matlab'])
    plt.show()
    
    #Figure: u[1] VB control input
    plt.plot(xi, u_opt[1]*ones(MPCit+1), 'r', xi, uAll[1,],'gs',xi, uAll_mat[1,],'bo', xi, uAll_pf[1,],'k*',xi, uAll_pfmat[1,],'c^')
    plt.ylabel('VB [m^3/min]')
    plt.xlabel('Time [min]')
    plt.title('Control input for VB')
    plt.legend(['Steady-state', 'iNMPC-Python', 'iNMPC-Matlab', 'pfNMPC-Python', 'pfNMPC-Matlab'])
    plt.show()
    
    #Figure: u[2] F control input
    plt.plot(xi, u_opt[2]*ones(MPCit+1), 'r', xi, uAll[2,],'gs', xi, uAll_mat[2,],'bo', xi, uAll_pf[2,],'k*', xi, uAll_pfmat[2,],'c^')
    plt.ylabel('F [kmol/min]')
    plt.xlabel('Time [min]')
    plt.title('Control input for F')
    plt.legend(['Steady-state', 'iNMPC-Python', 'iNMPC-Matlab', 'pfNMPC-Python', 'pfNMPC-Matlab'])
    plt.show()
    
    #Figure: u[3] D control input
    plt.plot(xi, u_opt[3]*ones(MPCit+1), 'r', xi, uAll[3,],'gs',xi, uAll_mat[3,],'bo', xi, uAll_pf[3,],'k*', xi, uAll_pfmat[3,],'c^')
    plt.ylabel('D [kmol/min]')
    plt.xlabel('Time [min]')
    plt.title('Control input for D')
    plt.legend(['Steady-state', 'iNMPC-Python', 'iNMPC-Matlab', 'pfNMPC-Python', 'pfNMPC-Matlab'])
    plt.show()
    
    #Figure: u[4] B control input
    plt.plot(xi, u_opt[4]*ones(MPCit+1), 'r', xi, uAll[4,],'gs', xi, uAll_mat[4,],'bo', xi, uAll_pf[4,], 'k*', xi, uAll_pfmat[4,],'c^')
    plt.ylabel('B [kmol/min]')
    plt.xlabel('Time [min]')
    plt.title('Control input for B')
    plt.legend(['Steady-state', 'iNMPC-Python', 'iNMPC-Matlab', 'pfNMPC-Python', 'pfNMPC-Matlab'])
    plt.show()
    
    ##=============================Plotting (Python Only)====================##
    
    #Figure: Objective function comparison (iNMPC vs pfNMPC)
    plt.plot(x,ObjReg, 'r', x, ObjReg_pf,'g')
    plt.title('Objective function: iNMPC vs pfNMPC')
    plt.xlabel('Number of MPC iteration [-]')
    plt.ylabel('Objective function [-]')
    plt.legend(['Objective Function: iNMPC', 'Objective Function: pfNMPC'])
    plt.show()
    
    #Figure: Concentration at stage 1 (reboiler)
    plt.plot(xi, xmeasureAll[0,],'r', xi, xmeasureAll_pf[0,], 'g')
    plt.ylabel('Concentration [-]')
    plt.xlabel('Time [min]')
    plt.title('Distillation: Bottom Composition')
    plt.legend(['iNMPC', 'pfNMPC'])
    plt.show()
    
    #Figure: Concentration at feed stage
    plt.plot(xi,xmeasureAll[NF,],'r', xi, xmeasureAll_pf[NF,],'g')
    plt.ylabel('Concentration [-]')
    plt.xlabel('Time [min]')
    plt.title('Distillation: Feed Composition')
    plt.legend(['iNMPC','pfNMPC'])
    plt.show()
    
    #Figure: Concentration at stage NT (top)
    plt.plot(xi,xmeasureAll[NT,],'r', xi, xmeasureAll_pf[NT,],'g')
    plt.ylabel('Concentration [-]')
    plt.xlabel('Time [min]')
    plt.title('Distillation: Top Composition')
    plt.legend(['iNMPC', 'pfNMPC'])
    plt.show()

    #Figure: Concentration in CSTR
    plt.plot(xi,xmeasureAll[NT+1,],'r', xi, xmeasureAll_pf[NT+1,],'g')
    plt.ylabel('Concentration [-]')
    plt.xlabel('Time [min]')
    plt.title('CSTR: Concentration')
    plt.legend(['iNMPC', 'pfNMPC'])
    plt.show()
              
    #Figure: Holdup in CSTR
    plt.plot(xi,xmeasureAll[2*NT-1,],'r',xi, xmeasureAll_pf[2*NT-1,], 'g')
    plt.ylabel('Holdup [-]')
    plt.xlabel('Time [min]')
    plt.title('CSTR: Holdup')
    plt.legend(['iNMPC', 'pfNMPC'])
    plt.show()
    
    #Figure: u[0] LT control input
    plt.plot(xi, uAll[0,],'r', xi, uAll_pf[0,],'g')
    plt.ylabel('LT [m^3/min]')
    plt.xlabel('Time [min]')
    plt.title('Control input for LT')
    plt.legend(['iNMPC', 'pfNMPC'])
    plt.show()

    #Figure: u[1] VB control input
    plt.plot(xi, uAll[1,],'r', xi, uAll_pf[1,],'g')
    plt.ylabel('VB [m^3/min]')
    plt.xlabel('Time [min]')
    plt.title('Control input for VB')
    plt.legend([ 'iNMPC', 'pfNMPC'])
    plt.show()

    #Figure: u[2] F control input
    plt.plot(xi, uAll[2,],'r', xi, uAll_pf[2,],'g')
    plt.ylabel('F [kmol/min]')
    plt.xlabel('Time [min]')
    plt.title('Control input for F')
    plt.legend(['iNMPC', 'pfNMPC'])
    plt.show()

    #Figure: u[3] D control input
    plt.plot(xi, uAll[3,],'r', xi, uAll_pf[3,],'g')
    plt.ylabel('D [kmol/min]')
    plt.xlabel('Time [min]')
    plt.title('Control input for D')
    plt.legend(['iNMPC', 'pfNMPC'])
    plt.show()

    #Figure: u[4] B control input
    plt.plot( xi, uAll[4,],'r', xi, uAll_pf[4,], 'g')
    plt.ylabel('B [kmol/min]')
    plt.xlabel('Time [min]')
    plt.title('Control input for B')
    plt.legend(['iNMPC', 'pfNMPC'])
    plt.show()

    print "Plotting is finished\n"
