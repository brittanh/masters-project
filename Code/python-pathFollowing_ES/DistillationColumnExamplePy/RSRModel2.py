__author__ = 'Vlad Minasides'
__date__ = 26.11

import numpy as np
import casadi as cd
import casadi.tools as cdtl


def idealRSR(structDaeModel=False, **kwargs):
    """ This is a muticomponent RSR model with a CSTR reactor,
        ideal distillation column with liquid hydrodynamics based on Francis Weir formula
        and distillate D recycled with a purge P split  """
    modelParameters = dict(kwargs["modelParameters"])

    assert(modelParameters["includeLiquidDynamics"] in [0,1])
    incLD = modelParameters["includeLiquidDynamics"]

    # Declaration of design parameters
    # Number of trays
    assert(isinstance(modelParameters["NT"], int))
    NT = modelParameters["NT"]
    # Number of components
    assert(isinstance(modelParameters["NC"], int))
    NC = modelParameters["NC"]
    # Feed tray number
    assert(isinstance(modelParameters["NF"], int))
    NF = modelParameters["NF"]
    # Relative volatilites
    assert(isinstance(modelParameters["alpha"], list) and len(modelParameters["alpha"])==NC-1)
    alpha = np.array(modelParameters["alpha"])
    # Caloric factor (p. 119)
    assert(0.0 <= modelParameters["q"] <= 1.0)
    q = modelParameters["q"]

    # Parameters for Franci's Weir Formula L(i) = K*(Mi -Mow)^1.5
    # Default values [1.625 , 1.177]
    Kuf = modelParameters["Kuf"]                          # Constant above feed
    Kbf = modelParameters["Kbf"]                          # Constant below feed

    # Need to change this value in case the flow rates change order of magnitude
    # Wier formula Li=Kf*(Mi-Muw)^1.5
    Muw = modelParameters["Muw"]                # Liquid holdup under weir (kmol)



#    B = 0; D = 0; F = 0; F_0 = 460.0/60.;
#    T_R = 340; z_F0 = [0.9];
    # Declaration of states
    states = cdtl.struct_symSX([
        cdtl.entry("x", shape=(NT+1,NC-1)),       # Liquid composition
        cdtl.entry("M", shape=(NT+1,1))])        # Molar holdup
    x, M = states[...]

    # Declaration of inputs
#    for key , value in inputs.
    if kwargs.has_key("controlLevels"):
        controlLevels = kwargs["controlLevels"]
        if controlLevels[0]:
            nx = kwargs["nx"]
            z0 = kwargs["z0"]

            inputs = cdtl.struct_symSX(["L_T", "V_B", "F_0", "z_F0"])
            L_T, V_B, F_0, z_F0 = inputs[...] #, F_0, T_R, P
            B = z0[nx+4] + controlLevels[1]*(M[0] - z0[nx/2])
            D = z0[nx+3] + controlLevels[2]*(M[nx/2-2] - z0[nx-2])
            F = z0[nx+5] + controlLevels[3]*(M[nx/2-1] - z0[nx-1])
        else:
            inputs = cdtl.struct_symSX(["L_T", "V_B","D", "B", "F", "F_0", "z_F0"]) #, "F_0", "T_R", "P"
            L_T, V_B, D, B, F, F_0, z_F0 = inputs[...] #, F_0, T_R, P
    else:
        inputs = cdtl.struct_symSX(["L_T", "V_B","D", "B", "F", "F_0", "z_F0"]) #, "F_0", "T_R", "P"
        L_T, V_B, D, B, F, F_0, z_F0 = inputs[...] #, F_0, T_R, P

    t = cd.SX.sym("t")
    y = cd.SX.sym("y", NT-1, NC-1)              # Vapor composition
    Li = cd.SX.sym("Li", NT, 1)                 # Liquid flow
    Vi = cd.SX.sym("Vi", NT, 1)                 # Vapor flow


    dMdt = cd.SX.sym("dMdt", NT+1, 1)             # Total molar holdup
    dMxdt = cd.SX.sym("dMxdt", NT+1, NC-1)        # Componentwise molar holdup
    dxdt = cd.SX.sym("dxdt", NT+1, NC-1)          # Rate of change of composition

#    dMrdt = cd.SX.sym("dMdt", 1, 1)             # Total CSTR componentwise molar holdup


    # Vapor flows (Constant, no dynamics) V[NT-1] skipped since D given
    for i in range(NT-1):
        if i<NF:
            Vi[i] = V_B
        else:
            Vi[i] = V_B # + (1-q)*F
    Vi[NT-1] = np.inf

    # Liquid flows (Wier formula) L[0] not required since B given
    Li[0] = np.inf
    for i in range(1, NT-1):
        if i<=NF:
            Li[i] = Kbf*(M[i] - Muw)**1.5*incLD + (1-incLD)*L_T
        else:
            Li[i] = Kuf*(M[i] - Muw)**1.5*incLD + (1-incLD)*L_T
    # Top tray liquid in
    Li[NT-1] = L_T


    # Vapor liquid equilibrium
    # (multicomponent ideal VLE, Stichlmair-Fair, 'Distillation', p. 36, 1998)
    for i in range(NT-1):
        for j in range(NC - 1):
            y[i,j] = x[i,j]*alpha[j]/(1 + cd.mtimes(alpha[:] - 1.0,x[i,:]))

#    z_F = x[NT,:]
#==============================================================================
#     # Molar balances
#==============================================================================

    # Level control
#    D =  498./60. - (M[NT-1] - 212.5)*1.1
#    B =  460./60. - (M[0] - 147.5)*1.1

    # Partial Reboiler (for simplicityy one can use L(0) = B)
    dMdt[0]     = Li[1] - Vi[0] - B
    for j in range(NC - 1):
        dMxdt[0,j]  = Li[1]*x[1,j] - Vi[0]*y[0,j] - B*x[0,j]

    # Stripping section
    for i in range(1,NF+1):
        dMdt[i]  = Li[i+1] - Li[i]+ Vi[i-1] - Vi[i]
        for j in range(NC-1):
            dMxdt[i,j] = Li[i+1]*x[i+1,j] - Li[i]*x[i,j] + Vi[i-1]*y[i-1,j] - Vi[i]*y[i,j]

    # Correction for the feed stage:
    dMdt[NF] += F
    for j in range(NC - 1):
        dMxdt[NF,j] +=  F*x[NT,j]

    # Enrichment section
    for i in range(NF+1,NT-1):
        dMdt[i]  = Li[i+1] - Li[i]+ Vi[i-1] - Vi[i]
        for j in range(NC - 1):
            dMxdt[i,j] = Li[i+1]*x[i+1,j] - Li[i]*x[i,j] + Vi[i-1]*y[i-1,j] - Vi[i]*y[i,j]

    # Total Condenser (for simplicityy one can use V(NT-1) = D )
    dMdt[NT-1] = Vi[NT-2] - Li[NT-1] - D
    for j in range(NC - 1):
        dMxdt[NT-1,j] = Vi[NT-2]*y[NT-2,j] - Li[NT-1]*x[NT-1,j]  - D*x[NT-1,j]

    # CSTR model
    k1 = 0.341/60.0 #1e8*np.exp(-6.67e4/(8.14*T_R))
    dMdt[NT] = F_0 + D - F
    for j in range (NC-1):
        dMxdt[NT,j] = F_0*z_F0[j] + D*x[NT-1,j] - F*x[NT,j] - k1*M[NT]*x[NT,j]

    for i in range(NT+1):
        for j in range(NC - 1):
            dxdt[i,j] = (dMxdt[i,j]-x[i,j]*dMdt[i])/M[i]

#==============================================================================
#   Construct sumbolic expressions
#==============================================================================

    # Construct helpers

    # Construct an SX function

    if structDaeModel:

        # Get the state structure for rhs assignment
        # Assign right hand side
        rhs = cdtl.struct_SX(states)
        rhs = cd.vertcat(dxdt,dMdt)
        print states
        print inputs

        return t, states, inputs, rhs

    else:
        lhs = cd.vertcat(x,M)
        rhs = cd.vertcat(dxdt,dMdt)

        return t, states, inputs, rhs

if __name__ == "__main__":

    ## Column model parameters
    NC = 2; NF = 12; NT = 22; q = 1.0; alpha = [2.0]
    includeLiquidDynamics = 1;  # Liquid dynamics Weir Formula L(i) = K*(Mi -Muw)^1.5
    scaleM_i = 1.0
    Kuf = [1.177*scaleM_i**1.5]; Kbf = [1.625*scaleM_i**1.5];  Muw = [0.25/scaleM_i] ;

    modelParameters = {
    "NT"    : NT,
    "NC"    : NC,
    "NF"    : NF,
    "alpha" : alpha,
    "q"     : q,
    "Kuf"   : Kuf,
    "Kbf"   : Kbf,
    "Muw"   : Muw,
    "includeLiquidDynamics" : includeLiquidDynamics,
    }

    t, lhs, inputs, rhs = idealRSR(modelParameters=modelParameters)

    print t, lhs, inputs, rhs
    print lhs["x"][0]

    temp_rhs = lhs["x"][0]

    z0 = np.load("./nominalValuesRSR_1010_from_sim.npy")

    # print z0[:46]
    x0 = np.array((z0[:46]), ndmin = 2)
    p = np.array((z0[46:]), ndmin=2)

    # Construct DAE model
    DAEModel = cd.SXFunction(cd.daeIn(t=t, x = lhs, p = inputs), cd.daeOut(ode = rhs))

    DAEModel.init()

    tf = 2000
    tsteps = 4000

    #! Simulate over Tf seconds, 100 timesteps.
    ts = np.linspace(0,tf,tsteps)

    integrator = cd.Integrator("cvodes", DAEModel)

    sim=cd.Simulator(integrator,ts)
    sim.init()

    sim.setInput(x0.T, "x0")
    sim.setInput(p.T, "p")

    sim.evaluate()
    x_results = np.array(sim.getOutput())

    print x_results







