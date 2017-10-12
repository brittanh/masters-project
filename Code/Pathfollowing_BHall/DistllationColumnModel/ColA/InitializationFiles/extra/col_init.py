#!/opt/local/bin/python
# -*- encoding: ascii -*-
"""
    @purpose: Generates steady state data for distillation column.
              Solves data in col_init.npy.
    @author: Brittany Hall
    @date: 06.10.2017
    @version: 0.1
    @updates:
    """
from numpy import *
from scikits.odes import ode
from col_model import *

def col_init():
    """
    Simulating stabilized LV-model for 20000 min to get initial steady state
    values for the distillation column
    """

    #Setting up stiff ode solver
#    ode15s = integrate.ode(f)
#    ode15s.set_integrator('vode', method='bdf', order=15, nsteps=3000)
#    t,x = ode15s.integrate('cola_main',[0 20000],0.5*transpose(ones((1,82))))
    t0 = 1                                                        #initial time
    x0 = 0.5*transpose(ones((1,82)))                             #intial states
    solution = ode('cvode',col_LV, old_api=False).solve(linspace(t0,20000),x0)
    t = solution.values.t
    x = solution.values.x[:,0]
    lengthx = size(x)
    Xinit = transpose(x[lengthx[0],:])
    
    print(x)
    raw_input()
    #Nominal inputs
    LT=2.70629                                                          #Reflux
    VB=3.20629                                                          #Boilup
    D=0.5                                                           #Distillate
    B=0.5                                                              #Bottoms
    F=1.0                                                             #Feedrate
    zF=0.5                                                    #Feed composition
    qF=1.0                                                #Feed liquid fraction
    Uinit = array([[LT], [VB], [D], [B], [F],[zF],[qF]])
    #Saving data to file
    save('col_init.npy',array([[lengthx],[Xinit],[Uinit]]))
    #d=load('col_init.npy')
    #print(Uinit==d)

col_init()
