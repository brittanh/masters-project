#  _________________________________________________________________________
#
#  Pyomo: Python Optimization Modeling Objects
#  Copyright (c) 2014 Sandia Corporation.
#  Under the terms of Contract DE-AC04-94AL85000 with Sandia Corporation,
#  the U.S. Government retains certain rights in this software.
#  This software is distributed under the BSD License.
#  _________________________________________________________________________

"""
This is the class to design and test iNMPC controllers
"""

from __future__ import division

__all__ = ['IdealNMPC']

from pyomo.core.base import *
from pyomo.dae import *
from pyomo.dae.plugins.colloc import Collocation_Discretization_Transformation
from pyomo.dae.plugins.finitedifference import Finite_Difference_Transformation
from pyomo.opt import SolverFactory
from pyomo.opt import SolverStatus, TerminationCondition
from decimal import Decimal
import time

## =============================================================================
## NEW IMPORTS FOR (N)MPC
## =============================================================================
from Disturbance import *
from InitialCondition import *
from ManipulatedVar import *

# New import from plotting
from DynamicPlot import *
    
class IdealNMPC():
    """
    The definition of an iNMPC class. 
    ideal Nonlinear Model Predictive Control    
    
    Public attributes:
    *disturbances: set of all the disturbances
    *manipulated_var: set of all the manipulated variables
    *state_var: set of all the state variables
    *initial_conditions: set of all the intial conditions
    	
    The default configuration is to discretize everything using the 
    LAGRANGE-RADAU representation with three collocation points and 
    prediction_horizon finite elements. It also considers the distrubances to 
    be constant in the unless it is specified the oposite later.
    """
    def __init__(self, model=None, model_time=None, prediction_horizon=0, 
    	    control_horizon=0, simulation_periods=0, 
    	    default_discretization=True):    	
	
	if model == None:
		raise TypeError("A pyomo model must be specified")
#	if type(model) is not ConcreteModel:
#		raise TypeError("The component specified using the 'model'" + \
#			" keyword must be a concrete Pyomo model")		
	self._model = model
	
	if model_time is None:
            raise TypeError("A differential set must be specified for the time")
        if model_time.type() is not ContinuousSet:
            raise TypeError("The component specified using the 'model_time'" + \
            	    " keyword must be a differential set")
        dt = model.find_component(model_time.cname(True))
        
        if dt is None:
            raise ValueError("DifferentialSet '%s' is not a valid component" + \
            	    " of the discretized model instance" %(diffset.cname(True)))
        self._model_time = dt        	
        
        if type(prediction_horizon) is not int or int(prediction_horizon) <= 0:
            raise TypeError("The prediction horizon must be a positive" + \
            	    " integer greater than zero")
        self._prediction_horizon = prediction_horizon
        
        if type(control_horizon) is not int or int(control_horizon) <= 0:
            raise TypeError("The control horizon must be a positive integer" + \
            	    " greater than zero")
        if control_horizon > prediction_horizon:
            raise ValueError("The prediction horizon must be greater or" + \
            	    " equal to the control horizon")
        self._control_horizon = control_horizon
    	
    	if type(simulation_periods) is not int or int(simulation_periods) <= 0:
            raise TypeError("The number of simulation periods must be a" + \
            	    " positive integer greater than zero")
        self._simulation_periods = simulation_periods    	    
    	    	
    	self._sampling_time = max(self._model_time)#/self._prediction_horizon    	        
	self._model_time.add(self._sampling_time*self._prediction_horizon)
	self._model_time.discard(self._sampling_time)    		
    	
    	self.points_in_fe = 0
    	self._discretize_method = ''
    	
    	self._real_time = []
    	self.disturbances = {}    	
    	self.manipulated_var = {}
    	self.state_var = {}
    	self.initial_conditions = {}    	
    	self.data_type_classification() 	 	
	
    	if default_discretization:    	    	    	
            self.discretize_model(discretization_type='OrthogonalCollocation',
            	    wrt=self._model_time, ncp=3, scheme='LAGRANGE-RADAU')                   
 		
    def discretize_model(self, **kwds):
    	"""
    	This function applies the discretization scheme over the dynamic model
    	defined previously.
    	
    	Keyword Arguments:
    	discretization_type The desiered type of discretization to apply
    		      	    FiniteDifferences or OrthogonalCollocation
        ncp           The desired number of collocation points over each 
                      finite element.
        wrt           Indicates which ContinuousSet the transformation 
                      should be applied to. If this keyword argument is not
                      specified then the same scheme will be applied to all
                      ContinuousSets.
        scheme        Indicates which finite difference method to apply. 
                      Options are LAGRANGE-RADAU, LAGRANGE-LEGENDRE, or 
                      HERMITE-CUBIC for OrthogonalCollocation. BACKWARD, CENTRAL
                      or FORWARD for FiniteDifferences.
    	"""    	
    	tmpdisc_type = kwds.pop('discretization_type','OrthogonalCollocation')
	#tmpnfe = kwds.pop('nfe',10)
	tmpnfe = self._prediction_horizon
        tmpncp = kwds.pop('ncp',3)        
        tmpds = kwds.pop('wrt',None)
        
        if (tmpdisc_type == 'OrthogonalCollocation'):
            tmpscheme = kwds.pop('scheme','LAGRANGE-RADAU')
            _scheme_name = tmpscheme.upper()
            self._discretize_method = Collocation_Discretization_Transformation()
            self._discretize_method.apply(self._model, wrt=tmpds, 
            	    nfe=tmpnfe, ncp=tmpncp, scheme=tmpscheme, inplace=True)
            self.points_in_fe = tmpncp
        if (tmpdisc_type == 'FiniteDifferences'):
            tmpscheme = kwds.pop('scheme','BACKWARD')
            _scheme_name = tmpscheme.upper()
            self._discretize_method = Finite_Difference_Transformation()
            self._discretize_method.apply(self._model, wrt=tmpds, 
            	    nfe=tmpnfe, scheme=tmpscheme, inplace=True)
            self.points_in_fe = 1    	
    	self.set_RealTime()    	 
    	# Modify objective based on the new continuous set
	for obj in self._model.component_objects(Objective, active=True):	    
	    obj_rule = obj.rule
	    obj_sense = obj.sense
	    obj_expr = obj.expr
	    obj._active = False
	self._model.NEW_OBJ = Objective(rule=obj_rule, expr=obj_expr, 
					sense=obj_sense)   	

    def data_type_classification(self):
    	"""
    	This function classify all the data types in the pyomo model in 
    	dictionaries that relate the name of the object with the object itself
    	"""
    	for obj in self._model.component_objects(Var, active=True):
    	    if (type(obj) is ManipulatedVar):
    	        self.manipulated_var[cname(obj)] = obj
    	    elif isinstance(obj, Var):
    	        self.state_var[cname(obj)] = obj
    	for obj in self._model.component_objects(Param, active=True):
    	    if (type(obj) is Disturbance):
    	        self.disturbances[cname(obj)] = obj
    	    elif (type(obj) is InitialCondition):
    	    	self.initial_conditions[cname(obj)] = obj       

    def set_PieceWiseStepManipulatedVars(self):    
    	"""
    	This function uses the reduce_collocation_points method from the dae
    	module to transform the manipulated variables in piece wise step
    	functions
    	"""    	
    	if (type(self._discretize_method) == Collocation_Discretization_Transformation): 
    	    for ManVar in self.manipulated_var.itervalues():    	       	    
    	        self._discretize_method.reduce_collocation_points(self._model, 
    	        	var=ManVar, ncp=1, diffset=self._model_time)

    def set_ControlHorizonConstraints(self):
    	"""
    	This function adds the constraints of the control horizon over all the
    	manipulated variables. In other words, it set their value constant
    	after the control horizon length.
    	"""
        for ManVar in self.manipulated_var.itervalues():
            cons_name = ManVar.name + "_ControlHorizon_constraint"
            initial_time = min(self._model_time)
            final_time = max(self._model_time)    
            control_time = self._control_horizon*self._sampling_time           
            def _ControlHorizonConstraint(m, l):
                j_tmp = sorted(ManVar.keys())
                if (l == initial_time):
                    return ManVar[l] - ManVar[j_tmp[j_tmp.index(l)+1]] == 0
                elif (l >= control_time and 
                      Decimal(str(l))%Decimal(str(self._sampling_time)) == 0 and 
                      l < final_time):                   
                    return ManVar[l] - ManVar[j_tmp[j_tmp.index(l)+1]] == 0
                else:
                    return Constraint.Skip
            self._model.add_component(cons_name, Constraint(self._model_time,
            	    rule = _ControlHorizonConstraint))
        
    #def add_InitialConditions(self):
    #    """
    #    Add the initial condition equations and artificial variables 
    #    to the model
    #    """
    #    for initC in self.initial_conditions.itervalues():
    #        var_name = initC.name + '_Artificial'            
    #        initC_sets = initC.get_sets()
    #        cs = self._model.dum3, self._model.dum2,
    #        self._model.add_component(var_name, Var(cs))        
    #	for obj in self._model.component_objects(Var, active=True):
    #       print obj
    #       for i in obj:
    #           print i
        
    def set_RealTime(self):
    	"""
    	This function sets the real time of the NMPC controller based on the 
    	sampling time and number of simulation periods
    	"""
    	self._real_time = [tmp_time + (period*self._sampling_time) \
    	    for period in range(self._simulation_periods) \
    	    for tmp_time in self._model_time if tmp_time <= self._sampling_time]
    	    
    def get_StVarSetNumber(self, a_StVar):
	"""
	This function returns the number of sets of a state variable and if
	the time set is included in it or not
	"""
	number_of_set = a_StVar.dim()	
	time_in_sets = False
	if (number_of_set == 1 and a_StVar._index == self._model_time):
	    time_in_sets = True
	elif (number_of_set >= 2 and 
		self._model_time in a_StVar._implicit_subsets):    
	    time_in_sets = True	       
        return number_of_set, time_in_sets
        
    def write_StVar(self, StVar, model_time):
        """
        This functions returns a tuple of two string where the first one are the
	names of the state variables and the second one the values of these 
	variables for the time specified        
        """
        var_titles = ''
        var_values = ''
        var_set_num, var_has_time = self.get_StVarSetNumber(StVar)
        if (var_set_num == 0):
            var_titles = var_titles + str(StVar) + ' '
            var_values = var_values + str(StVar.value) + ' '
        elif (var_set_num == 1):
            if (var_has_time):
                var_titles = var_titles + str(StVar) + ' '
                var_values = var_values + str(StVar[model_time].value) + ' '	
            else:
                for indx in StVar._index:
                    var_titles = var_titles + str(StVar) + \
                    	'[' + str(indx) + ']' + ' '
                    var_values = var_values + str(StVar[indx].value) + ' ' 	              	
        elif (var_set_num >= 2):     
            if (var_has_time):
                time_index = StVar._implicit_subsets.index(self._model_time)		  	                                
                for indx in StVar._index:
                    indx_l = list(indx)
                    if (indx_l[time_index] == model_time): 
                    	indx_l[time_index] = 'Time'
                        var_titles = var_titles + str(StVar) + '[' + \
                            ', '.join(map(str, indx_l)).replace(' ', '') + \
                            ']' + ' '	    
                    	indx_l[time_index] = model_time  
                        var_values = var_values + \
                        	str(StVar[tuple(indx_l)].value) + ' '      
            else:
                for indx in StVar._index:
                    var_titles = var_titles + str(StVar) + '[' + \
                        ', '.join(map(str, indx)).replace(' ', '') + ']' + ' '	
                    var_values = var_values + str(StVar[indx].value) + ' '    
        return var_titles, var_values                             
    	    
    def solve(self, solver_name='ipopt', options=[], file_name='untitled', 
    	    dynamic_plot=False):
    	"""
    	This function solves the iNMPC over all the simulation periods stablished
    	"""    	
    	#Check the validity of the input data
    	for ManVar in self.manipulated_var.itervalues():
            ManVar.check_ManipulatedVar(self._model_time)
        for dist in self.disturbances.itervalues():
            dist.check_disturbance(self._model_time)
        for init_cond in self.initial_conditions.itervalues():
            init_cond.check_InitialCondition(self._model_time)
    	#It is neccesary to add the reamining components to the model    
        self.set_PieceWiseStepManipulatedVars()
        self.set_ControlHorizonConstraints()
        #Create solver
        opt = SolverFactory(solver_name)
        #Solver options
        for op in options:
            opt.options[op[0]] = op[1]
        #Open file to save results    	    	
	profile_results_list = []	    	
    	computational_results = open(file_name+'_iNMPC_computational_results.txt', 'w')
    	#Write titles in profile file
    	ans = 'Time '
    	for dis in self.disturbances.itervalues():
            ans = ans + str(dis) + ' '
        for ManVar in self.manipulated_var.itervalues():
            ans = ans + str(ManVar) + ' '
        for StVar in self.state_var.itervalues():  
            ans = ans + self.write_StVar(StVar, min(self._model_time))[0]            
	ans = ans + '\n'
	number_disturbances = len(self.disturbances)
	number_ManipulatedVar = len(self.manipulated_var)
	number_StateVar = (len(ans.split(' ')) - len(self.disturbances) \
		         - len(self.manipulated_var) - 2)
	profile_results_list.append(ans)	
    	computational_results.write('Periods TerminationCondition' +\
    		' SimulationTime[secs] ObjectiveFunction \n')
    	#Create plots
    	if (dynamic_plot):    	                   
            ydata_plot = {}
            dynamic_plots = {}
            plot_names = ['Disturbances', 'Manipulated Variables',
                'State Variables']
            plot_legends = {'Disturbances':\
            	    ans.split(' ')[1:number_disturbances+1],
            	    'Manipulated Variables':\
            	    ans.split(' ')[number_disturbances+1:number_disturbances+\
            	    			number_ManipulatedVar+1],
            	    'State Variables':\
            	    ans.split(' ')[number_disturbances+number_ManipulatedVar+1:\
            			   number_disturbances+number_ManipulatedVar+\
            			   number_StateVar+1] }
	    plot_total_variables = {'Disturbances':number_disturbances,
            	    'Manipulated Variables':number_ManipulatedVar,
            	    'State Variables':number_StateVar }            	    
            for i in plot_names:
                dynamic_plots[i] = DynamicPlot(i, 'Time', '', 
                	plot_legends[i] , plot_total_variables[i])                
        #Solution
        for period in range(self._simulation_periods):            
            current_time = period*self._sampling_time
            #Update disturbances
            for dist in self.disturbances.itervalues():            	
                dist.update_disturbance(current_time)
            #Read initial conditions
            if (period == 0):
                pass # Do nothing for now
            else:
                for init_cond in self.initial_conditions.itervalues():
                    init_cond.update_InitialCondition(self._sampling_time)            
            #Solution
	    initial_time = time.time()	
            NMPC_results = opt.solve(self._model, tee=False)            
            final_time = time.time()                   
            print '------------------------------------------------------------'
            print 'Simulation period: ' + str(period+1)
            print 'Termination condition: ' +\
                str(NMPC_results.solver.termination_condition)
            print 'Objective function value: ' + str(value(self._model.NEW_OBJ))
            print 'Solution time: ' + str(final_time - initial_time) + ' secs'            
            #print 'Solution time: ' + str(NMPC_results.solver.Time) + ' secs'                        
            #Write computational results
            computational_results.write(str(period+1) + ' ' +\
            	    str(NMPC_results.solver.termination_condition) +\
            	    ' ' + str(final_time - initial_time) +\
            	    ' ' + str(value(self._model.NEW_OBJ)) + ' \n')
            #Write profile results in files
            for i in range(0, self.points_in_fe+1):                                    
                 t_m = self._real_time[i]
                 t_r = self._real_time[i] + current_time
                 ans = str(t_r) + ' '
                 for dis in self.disturbances.itervalues():
                     ans = ans + str(dis.get_realValue(t_r)) + ' '                     
                 for ManVar in self.manipulated_var.itervalues():
                     ans = ans + str(ManVar[t_m].value) + ' '  
                 for StVar in self.state_var.itervalues():
                     ans = ans + self.write_StVar(StVar, t_m)[1]      		                             
            	 ans = ans + '\n'                    	 
            	 profile_results_list.append(ans)            	        	 
            #Dynamic plot            
            if (dynamic_plot):      
            	#Reults matrix    		         	
            	results = [map(float,rr.split(' ')[:-1]) \
            		   for rr in profile_results_list \
            		   if profile_results_list.index(rr)>0]
            	results_transpose = [[row[i] for row in results] \
            				for i in range(len(results[0]))]            	            	            	    		        
		# Plots
		t_plot = results_transpose[0]
		ydata_plot['Disturbances'] = results_transpose[1:\
						number_disturbances+1]
		ydata_plot['Manipulated Variables'] = \
			results_transpose[number_disturbances+1:\
				number_disturbances+1+number_ManipulatedVar]
		ydata_plot['State Variables'] = \
			results_transpose[number_disturbances+1+\
				number_ManipulatedVar:\
				number_disturbances+1+number_ManipulatedVar+\
				number_StateVar]		
		for i in plot_names:
		    dynamic_plots[i].update_plot(t_plot, ydata_plot[i]) 			    		    
        #Write and close files
        profile_results = open(file_name+'_iNMPC_profiles_results.txt', 'w')
        profile_results.writelines(profile_results_list)
        profile_results.close()
        computational_results.close()          
        for i in plot_names:
                dynamic_plots[i].display_final_plot()
