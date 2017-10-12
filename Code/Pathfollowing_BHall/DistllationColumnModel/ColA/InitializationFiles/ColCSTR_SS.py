#!/opt/local/bin/python
# -*- encoding: ascii -*-
"""
    @purpose: Steady state optimization for CSTR and distillation column A
    Creates 'CstrDistXinit.mat', 'LambdaCstrDist.mat' and 'Qmax.mat'
    @author: Brittany Hall
    @date: 11.10.2017
    @version: 0.1
    @updates:
"""

from scipy.io import savemat
from casadi import *
from numpy import append, ones, transpose, shape, abs, size, concatenate, array
from scipy.linalg import eigvals
#from sympy import Matrix
from buildmodel import *
from params import * #imports cstr and distillation column parameters
from nlp_solve import *
import time

#Unpacking parameter values
NT = params['dist']['NT']
LT = params['dist']['LT']
VB = params['dist']['VB']
F = params['dist']['F']
D = params['dist']['D']
B = params['dist']['B']


#Symbolic
x = SX.zeros(2*NT+2,1)
l = SX.zeros(2*NT+2,1)
for i in range(0,2*NT+2):
   x[i] = SX.sym('x_'+ str(i+1))
   l[i] = SX.sym('l_'+ str(i+1))

u1 = SX.sym('u1')                                                           #LT
u2 = SX.sym('u2')                                                           #VB
u3 = SX.sym('u3')                                                            #F
u4 = SX.sym('u4')                                                            #D
u5 = SX.sym('u5')                                                            #B

#Collecting states and inputs
x = vertcat(x[:])
x = vertcat(x,u1)
x = vertcat(x,u2)
x = vertcat(x,u3)
x = vertcat(x,u4)
x = vertcat(x,u5)

#Decision variables (states and controls)
Xinit = 0.5*ones((2*NT+2,1))
Uinit = vertcat(Xinit, LT)
Uinit = vertcat(Uinit, VB)
Uinit = vertcat(Uinit, F)
Uinit = vertcat(Uinit, D)
Uinit = vertcat(Uinit, B)

#Define the dynamics as equality constraints and additional inequality constraints
obj, eq, lbx, ubx, lbg, ubg = buildmodel(x,params)
prob = {'f': obj, 'x': x, 'g': eq}
options = {}
tic = time.time()
startnlp = tic
w0 = Uinit
lbw = lbx
ubw = ubx
sol = nlp_solve(prob, options, w0, lbw, ubw, lbg, ubg)
toc = time.time() - tic
elapsednlp = toc
print('IPOPT solver runtime = %f\n', elapsednlp)

u = sol['x']
lam = sol['lam_g']
lam[NT+1:end] = -1*lam[NT+1:end]

Xinit = u
#Saving steady state data to be used in dynamic optimization
savemat('CstrDistXinit1', Xinit)
savemat('LamdaCstrDist1.mat', lam)

"""
Compute Hessian and perform Greshgorin convexification
"""
xsol = u
lamda = {}
lamda['eqnonlin'] = lam
l = vertcat(l[:])
L = obj + transpose(l)*eq

Lagr = Function('Lagr', [x,l],[L])
H = Function('H', [x,l], hessian(Lagr(x,l)))
cons = Function('Cons', [x], [eq])
Jcon = Function('Jcon', [x,cons], jacobian(cons(x,cons)))

eqVal = cons(xsol)
Hx = H(xsol, lamda['eqnonlin'])
Hx = full(Hx)
Jac = Jcon(xsol)
Jac = full(Jac)

#Nullspace of the constraints and its eigenvalue
rH = transpose(Jac.nullspace())*Hx*Jac.nullspace()
eigen_RH = eigvals(rH)

#Calculating the Greshgorin convexification
def Greshgorin(H):
    numH = H.shape[0]
    Q = zeros((numH,numH))
    delta = 2.5 #with measurement noise of 1 percent
    for i in range(0,numH): #iterate every row of the Hessian
        sumRow = 0
        for j in range(0,numH):
            if j != i:
                sumRow += abs(H[i,j])
        if H[i,i] <= sumRow: #include equality
            Q[i,i] = sumRow - H[i,i] + delta
    Q = diag(Q)
    return H, Q

Hxxl, Qmax = Greshgorin(Hx)
savemat('Qmax', Qmax)

#Check at some initial point for optimization
xstat = Xinit[0:2*NT+2]
u0 = array([[2.5],[3.5],[0.6],[0.5],[0.5]])
xeval = concatenate((xstat,u0))
Jeval = Jcon(xeval)
Jeval = full(Jeval)
Hxxl = H(xeval, lam['eqnonlin'])
Hxxl = full(Hxxl)
Hconv = Hxxl + diag(Qmax)
rHe = transpose(Jeval.nullspace())*Hconv*Jeval.nullspace()
