#!/opt/local/bin/python
# -*- encoding: ascii -*-
"""
    @purpose: Objective function to be solved
    @author: Brittany Hall
    @date: 08.10.2017
    @version: 0.1
    @updates:
"""
from casadi import *
from numpy import size, append, transpose
import scipy.io as spio
from distcolcstr_prob import itPredHor
from params import *

def objective(x,y,p,N):
    
    nPrimal = size(x)
    nDual = size(y['lam_g'])
    nParam = size(p)

    #Model parameters
    NT = params['dist']['NT']          #Number of stages in distillation column
    Uf = 0.3                                             #Feed rate to CSTR F_0
    ~,state,xdot,inputs = ColCSTR_model(Uf)
    sf = Function('sf',[state,inputs],[xdot])

    global nx nu nk d tf ns
    h = tf/nk

    #preparing collocation matrices
    ~,C,D,d = collocationSetup()

    #NLP variable vector
    V = {}                                #Decision variables (control + state)
    obj = 0                                                 #Objective function
    cons = {}
    delta_time = 1
    alpha = 1
    beta = 1
    gamma = 1

    #Loading initial states and controls
    data = spio.loadmat('CstrDistXinit.mat', squeeze_me = True)
    Xinit = data['Xinit']
    xf = Xinit[0:84]
    u_opt = Xinit[85:89]

    #Initial conditions
    X0 = MX.sym('X0', nx)
    V = append(V,X0)
    cons = append(cons, X0 - x[0:nx,0])
    cons_x0 = X0 - x[0:nx,0]

   #Formulating the NLP
   Xk = X0
   data = spio.loadmat('Qmax.mat', squeeze_me = True)
   param['Qmax'] = data['Qmax']
   ssoftc = 0
   for i in range(0,N):
       obj, cons, V, Xk, params, ssoftc = itPredHorizon(Xk, V, cons, obj, params, i, ssoftc)
                               
   V = vertcat(V[:])
   con = vertcat(cons[:])

   #Objective function and constraint functions
   f = Function('f', [V], [obj], chr('V'), chr('objective'))
   c = Function('c', [V], [cons], chr('V'), chr('constraint'))
   cx0 = Function('cx0',[X0], [cons_x0], chr('X0'), chr('constraint'))
                               
   #Constructing Lagrangian
   lag_expr = obj + transpose(y['lam_g'])*cons
   g = gradient(f)
   lagr = Function('lagr', [V], [lag_expr], chr('V'), chr('lag_expr'))
   H = Function('H',hessian(lagr,['V','lag_expr']))
   Hobj = hessian(f,['V','objective'])
   J = jacobian(c,['V','constraint'])
   Jp = jacovian(cx0,['X0','constraint'])
                               
   f = f(x)
   g = g(x)
   H = H(x)
   Lxp = H[0:nPrimal,0:nParam]
   J = J(x)
   Jtemp = zeros((nDual,nParam))
   cp = Jp(x[0:nParam])
   Jtemp[0:nParam, 0:nParam] = full(cp)
   cp = spares(Jtemp)
   cst = c(x)
   
   #Evaluation of objective function used for Greshgorin bound
   Hobj = Hobj(x)
   Hobj = sparse(Hobj)
   f = full(f)
   g = sparse(g)
   H = sparse(H)
   Lxp = sparse(Lxp)
   J = sparse(J)
   cp = sparse(cp)
   cst = full(cst)
   
   #Equality constraint
   Jeq = J
   dpe = cp
   
   return f, g, H, Lxp, cst, J, cp, Jeq, dpe, Hobj
                               
                               
                   

