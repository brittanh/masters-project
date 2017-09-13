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

# Imports for iNMPC
from pyomo.core.base import *
from pyomo.core.base.var import SimpleVar
from pyomo.dae import *
from pyomo.dae.plugins.colloc import Collocation_Discretization_Transformation
from pyomo.dae.plugins.finitedifference import Finite_Difference_Transformation
from pyomo.opt import SolverFactory
from pyomo.opt import SolverStatus, TerminationCondition
from decimal import Decimal

# Imports for nl writing
from pyomo.opt import ProblemFormat
from pyomo.core.base import SymbolMap
import re
import time
from subprocess import call
from pyomo.opt.plugins.sol import ResultsReader_sol

## =============================================================================
## NEW IMPORTS FOR (N)MPC
## =============================================================================
from Disturbance import *
from InitialCondition import *
from ManipulatedVar import *

# New import from plotting
from DynamicPlot import *
        
    
def find_comp(name, pyomo_model):
    """ 
    A function that returns a component of a pyomo model given its name as a
    string including the corresponding indexes
    """
    name_tmp = re.search("(.*)\[(.*)\]",name)   
    if name_tmp is None:
        component = pyomo_model.find_component(name)
        pyomo_object = component
        ix = ''
    else:
        component = pyomo_model.find_component(name_tmp.groups()[0] + '[**]')
        ix = name_tmp.groups()[1]
	for index in component:          	    
	    ind = ''
	    try: 
	        for i in range(len(index)):
		    if i == len(index) - 1:
			ind += str(index[i])
		    else:		   	    	   	    
			ind += str(index[i]) + ','	    	    	    	        	    	    	    
	    except TypeError:
		ind = str(index)
	    if ind == ix:    	    	    
		pyomo_object = component[index]   		    
    return pyomo_object    

class nl_map_iNMPC():
    """
    The definition of a class that creates and modify the nl file of an 
    optimization problem
    """
    
    def __init__ (self, pyomo_model, file_name='untitled'):
        """
	It initialize the objects of this class creating the .nl, .col, .row 
	files associated with the pyomo model. 
	It also creates the symbolic map and dictionaries that relates pyomo's 
	object with nl names
	It creates a dictionary between the artificial constraints and the 
	corresponding variable
	"""
	# Important atributes
	self.name = file_name	
	self.nl_change = False
	self.x_numvar = None
	self.d_numvar = None
	self.cons_num = None
	self.x_line = None
	self.d_line = None
	self.calls2update_Parameters = 0	
	self.r_line = None	
	self.call2update_InitialGuess = 0
	self.var_dict = {}
	self.num_var_lo = 0
	self.num_var_up = 0
	self.dual_dict = {}
	# Write .nl .col .row files			
	t00 = time.time()
	pyomo_model.write(file_name +'.nl', format = ProblemFormat.nl, 
		io_options = {'symbolic_solver_labels':True,
				'file_determinism':3})
	self.nl_writting_time = time.time() - t00
	#Creating mapping and dictionaries	
	t00 = time.time()
	with open(file_name +'.col') as f:
	    varnames = [line.strip() for line in f.readlines()]
	with open(file_name +'.row') as f:
	    connames = [line.strip() for line in f.readlines()]
	objname = connames.pop()				
	#Create symbol map identical to that generated by the nl file interface
	self.symbol_map = SymbolMap()				
	symbol_list = []	
	self.pyomo2nl = {} 
	self.nl2pyomo = {} 
	self.pyomoName2pyomoObject = {}
	for cntr, name in enumerate(varnames):
	    symbol = 'v' + str(cntr)
	    component = find_comp(name, pyomo_model)
	    symbol_list.append((component, symbol))   
	    self.pyomo2nl[name] = 'V' + str(cntr) 	
	    self.nl2pyomo['V' + str(cntr)] = name
	    self.pyomoName2pyomoObject[name] = component
	for cntr, name in enumerate(connames):
	    symbol = 'c' + str(cntr)
	    component = find_comp(name, pyomo_model)
	    symbol_list.append((component, symbol))   
	    self.pyomo2nl[name] = 'C' + str(cntr)
	    self.nl2pyomo['C' + str(cntr)] = name
	    self.pyomoName2pyomoObject[name] = component
	symbol = 'o'+str(0)
	component = pyomo_model.find_component(objname)
	symbol_list.append((component, symbol))		
	self.symbol_map.addSymbols(symbol_list) 
	# Create additional dictionaries
	self.update_InitialGuess(pyomo_model, var_class = True)
	self.var_classification()
	#Dictionary for artificial constraint and variable	
	self.artificial_classification()
	self.nl_dictWritting_time = time.time() - t00
	
    def artificial_classification(self):
    	"""
    	It creates a dictionary between the artificial variables and its
    	corresponding artificial constraint
    	"""
    	self.var2artcons = {}
	self.artcons2var = {}
	id2var = {}
	nl_file = open(self.name + '.nl','r')
	nl_lines = nl_file.readlines()
	nl_file.close()
	var_indicator = re.search("(S4 )(\d\d*)( sens_state_0)", 
					''.join(nl_lines))  
	var_num = int(var_indicator.groups()[1])                                        
	var_line = nl_lines.index(''.join(var_indicator.groups()) + '\n')			
	cons_indicator = re.search("(S5 )(\d\d*)( sens_init_constr)", 
					''.join(nl_lines))			
	cons_num = int(cons_indicator.groups()[1])
	cons_line = nl_lines.index(''.join(cons_indicator.groups()) + '\n')			
	if var_num == cons_num:
	    for line_number in range(var_line + 1, var_line + var_num + 1, 1):
	        line = nl_lines[line_number].strip().split()					
		var_name = self.nl2pyomo['V' + str(line[0])]
		id_var = line[1]					
		id2var[id_var] = var_name
	    for line_number in range(cons_line + 1, cons_line + cons_num + 1,1):
		line = nl_lines[line_number].strip().split()					
		cons_name = self.nl2pyomo['C' + str(line[0])]
		var_name = id2var[line[1]]
		self.var2artcons[var_name] = cons_name
		self.artcons2var[cons_name] = var_name				
		nl_lines[line_number] = line[0] + ' ' + '1' + '\n'										
	else:
	    raise ValueError("Error in the suffixes " +\
				"denfition within the model")	
	# Overwrite nl file
	nl_file = open(self.name + '.nl', 'w')
	nl_file.seek(0,0)		
	nl_file.writelines(nl_lines)			
	nl_file.close()		 
            			
    def var_classification(self):
    	"""
    	It creates a dictionary between the variable code (V###) and a list
    	corresponding to its current value, its lower bound and its upper bound
    	"""
    	#Variables with lower and upper bounds    	    	
        nl_file = open(self.name + '.nl','r')
	nl_lines = nl_file.readlines()
	nl_file.close()
	#self.x_numvar = re.search("(x)(\d\d*)", ''.join(nl_lines))
	#self.x_line = nl_lines.index('x' + self.x_numvar.groups()[1] + "\n")
	self.b_line = nl_lines.index("b\n")		
        for i in range(int(self.x_numvar.groups()[1])):
            var_number = i
            var_value = float(nl_lines[self.x_line + i + 1].split()[1])
            var_bounds = nl_lines[self.b_line + i + 1].split()
            if var_bounds[0] == '0':
                self.var_dict['V' + str(i)] = [var_value, 
                		float(var_bounds[1]), float(var_bounds[2])]
                self.num_var_lo += 1
                self.num_var_up += 1		
            elif var_bounds[0] == '1':
                self.var_dict['V' + str(i)] = [var_value, 
                		None, float(var_bounds[1])]                
                self.num_var_up += 1
            elif var_bounds[0] == '2':
            	self.var_dict['V' + str(i)] = [var_value, 
            			float(var_bounds[1]), None]
            	self.num_var_lo += 1                
            elif var_bounds[0] == '3':
            	self.var_dict['V' + str(i)] = [var_value, 
            			None, None] 
	       
    def update_InitialGuess(self, pyomo_model_with_results, var_class=False, update_duals=False):
        """
        Update the initial values (guess) of all the primal variables 
        in the nl file
        It uses the values set in the correcponding pyomo model
        It can also initialize the dual variables
        """
	pmwr = pyomo_model_with_results		
	# Read nl file as lines
	nl_file = open(self.name + '.nl', 'r')		
	nl_lines = nl_file.readlines()	
	nl_file.close()
	if self.call2update_InitialGuess == 0 or self.nl_change:
	# Identify number of variables (x primal, d dual)
    	    self.x_numvar = re.search("(x)(\d\d*)", ''.join(nl_lines))					
    	    self.d_numvar = re.search("(d)(\d\d*)", ''.join(nl_lines))
    	    self.cons_num = int(nl_lines[1].strip().split()[1])
    	    self.var_num = int(nl_lines[1].strip().split()[0])			
    	    # Identify index of x in list	
    	    self.x_line = nl_lines.index('x' + self.x_numvar.groups()[1] + "\n")
    	    self.nl_change = False
        # Read and modify lines after xnn
        if int(self.x_numvar.groups()[1]) < self.var_num or var_class:
            self.r_line = nl_lines.index("r\n")	
	    nl_lines_1 = nl_lines[:self.x_line+1]
	    nl_lines_2 = nl_lines[self.r_line:]
	    nl_lines_1.pop(-1)	    
	    nl_lines_1.append('x' + str(self.var_num) + "\n")
	    for variable_number in range(0, self.var_num, 1):
	        pyomo_var_name = self.nl2pyomo['V' + str(variable_number)]			
	        pyomo_var_value = 0
	        nl_lines_1.append(str(variable_number) + ' ' +
	        			str(pyomo_var_value) + "\n")
	    del nl_lines
	    nl_lines = nl_lines_1 + nl_lines_2	
	    self.x_numvar = re.search("(x)(\d\d*)", ''.join(nl_lines))
	    self.x_line = nl_lines.index('x' + self.x_numvar.groups()[1] + "\n")
	    self.nl_change = True
    	else:			
	    for line_number in range(self.x_line + 1, 
	    	    	self.x_line + int(self.x_numvar.groups()[1]) + 1, 1):			
	        line = nl_lines[line_number].strip().split()			
	        pyomo_var_name = self.nl2pyomo['V' + line[0]]			
	        pyomo_var_value = self.pyomoName2pyomoObject[pyomo_var_name].value
		if pyomo_var_value is None:
		    pyomo_var_value = 0
	        value_bounds = self.var_dict['V' + line[0]]
	        self.var_dict['V' + line[0]] = [pyomo_var_value, 
	        				value_bounds[1], value_bounds[2]] 
	        nl_lines[line_number] = line[0] + ' ' + str(pyomo_var_value) + "\n"			
	# Check the update duals option
	if update_duals:
	    # Check if d is in the nl file			
	    if self.d_numvar == None:			
	        nl_lines_1 = nl_lines[:self.x_line]
		nl_lines_2 = nl_lines[self.x_line:]
		nl_lines_1.append('d' + str(self.cons_num) + "\n")
		for constraint_number in range(0, self.cons_num, 1):
		    pyomo_cons_name = self.nl2pyomo['C' + str(constraint_number)]			
		    pyomo_dual_value = pmwr.dual[self.pyomoName2pyomoObject[pyomo_cons_name]]
		    self.dual_dict['C' + str(constraint_number)] = pyomo_dual_value
		    nl_lines_1.append(str(constraint_number) + ' ' +
		    	    			str(pyomo_dual_value) + "\n")
		    
		del nl_lines
		nl_lines = nl_lines_1 + nl_lines_2	
		self.d_numvar = re.search("(d)(\d\d*)", ''.join(nl_lines))
		self.d_line = nl_lines.index('d' + self.d_numvar.groups()[1] + "\n")
		self.nl_change = True
	    # If d is already in the nl file
	    else:		
		if self.call2update_InitialGuess == 0: 
		    self.d_line = nl_lines.index('d' + self.d_numvar.groups()[1] + "\n")
	        # Read and modify lines after dnn
	        for line_number in range(self.d_line+1, 
	        	self.d_line + int(self.d_numvar.groups()[1]) + 1, 1):
		    line = nl_lines[line_number].strip().split()
		    pyomo_cons_name = self.nl2pyomo['C' + line[0]]					
		    pyomo_dual_value = pmwr.dual[self.pyomoName2pyomoObject[pyomo_cons_name]]
		    self.dual_dict['C' + line[0]] = pyomo_dual_value
		    nl_lines[line_number] = line[0] + ' ' + str(pyomo_dual_value) + "\n"
	        self.nl_change = False
        # Overwrite nl file
	nl_file = open(self.name + '.nl', 'w')
	nl_file.seek(0,0)		
	nl_file.writelines(nl_lines)			
	nl_file.close()	
	self.call2update_InitialGuess += 1	
	
    def update_Parameters(self, dict_NewValues):
        """
	It updates the parameters of the model that are defined as artificial 
	constraints (e.g. initial conditions and disturbances)
	"""					
	# Read nl file
	nl_file = open(self.name+'.nl','r')
	nl_lines = nl_file.readlines()
	nl_file.close()	 
	# Modify values
	if self.calls2update_Parameters == 0 or self.nl_change:
	    self.InitialParam2line = {} # Dictionary to save lines modified
	    self.r_line = nl_lines.index("r\n")
	    if self.cons_num is None:
	        self.cons_num = nl_lines[1].strip().split()[1]
	    for line_number in range(self.r_line+1,self.r_line+1+int(self.cons_num)):
	        C_ID = line_number - (self.r_line+1)
		cons_name = self.nl2pyomo['C'+str(C_ID)]
		if cons_name in self.artcons2var:
		    var_name = self.artcons2var[cons_name]		    
		    nl_lines[line_number] = '4 ' + str(dict_NewValues[var_name][2]) + "\n"					
		    self.InitialParam2line[var_name] = line_number
	else:
	    for var_name in self.InitialParam2line.keys():
	        line_number = self.InitialParam2line[var_name]				
		nl_lines[line_number] = '4 ' + str(dict_NewValues[var_name][-1]) + "\n"
	# Overwrite nl file
	nl_file = open(self.name+'.nl', 'w')
	nl_file.seek(0,0)		
	nl_file.writelines(nl_lines)			
	nl_file.close()	
	self.calls2update_Parameters += 1
	
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
	# Dictionary for updating the parameters
    	self.gen_dictUpdateParameters() 	    	

    def data_type_classification(self):
    	"""
    	This function classify all the data types in the pyomo model in 
    	dictionaries that relate the name of the object with the object itself    	
    	"""
    	for obj in self._model.component_objects(Var, active=True):
    	    if (type(obj) is ManipulatedVar):
    	        self.manipulated_var[cname(obj)] = obj
    	    elif (type(obj) is SimpleInitialCondition or type(obj) is IndexedInitialCondition):
    	        self.initial_conditions[cname(obj)] = obj  
    	    elif (type(obj) is Disturbance):
    	        self.disturbances[cname(obj)] = obj
    	    elif isinstance(obj, Var):
    	        self.state_var[cname(obj)] = obj    	
    	        
    def gen_dictUpdateParameters(self):
        """
        A function that generates the corresponding dictionaries between the 
        artificial variables and their new value
        It is called after the model discretization
        Keys: name with indexes as a string
        Values: list with the obejct, the index, and its current value
        """        
        self.dictUpdateParameters = {}
        # Disturbances
        for dis in self.disturbances.itervalues():
            for t in self._model_time:
                name = dis.name + '[' + str(t) + ']'
		name = name.replace("'", "")
                values = dis.forecast(self._model_time, dis.rule)
                self.dictUpdateParameters[name] = [dis, t, values[t]]
        # Initial conditions
	for initCond in self.initial_conditions.itervalues():
	    if type(initCond) == SimpleInitialCondition:
	        name = initCond.name
		name = name.replace("'", "")
	        value = initCond.value	        
	        self.dictUpdateParameters[name] = [initCond, None, value]
	    else:
	        for indx in initCond._index:
		    if ',' in str(indx):		        	
	                name = initCond.name + '[' + \
	            		''.join(str(indx)[1:-1].split(' ')) + ']'
			name = name.replace("'", "")
	            else:
	            	name = initCond.name + '[' + str(indx) + ']'	            
			name = name.replace("'", "")
	            value = initCond[indx].value
	            self.dictUpdateParameters[name] = [initCond, indx, value]
	
    def update_dictUpdateParameeters(self, current_time):
    	"""
    	It updates the dictionary that updates the parameters in the nl file
    	with the new values
    	"""
    	time_set = [current_time + t for t in self._model_time]
    	for item_name in self.dictUpdateParameters:
    	    item = self.dictUpdateParameters[item_name]
    	    if type(item[0]) == Disturbance:    	            	    	     
                value = item[0].forecast(time_set, item[0].rule)[item[1] + current_time]                                
            else:
            	#print item[1], type(item[1])    
                if item[1] == None:
                    value = item[0].variable[self._sampling_time].value			                    
            	elif type(item[1]) == tuple:            	    	
            	    var_indx = item[0].var_indx(item[1], self._sampling_time)
            	    value = item[0].variable[var_indx].value
              	else:
              	    var_sets = item[0].variable._implicit_subsets	    	    
		    other_set = item[0]._index
		    other_set_position = var_sets.index(other_set)	      		      	
		    var_idx = [0]*2
		    var_idx[other_set_position] = item[1]
		    var_idx[1-other_set_position] = self._sampling_time
		    value = item[0].variable[tuple(var_idx)].value	
            self.dictUpdateParameters[item_name] = [item[0], item[1], value]	    
    
    def suffixes_declaration(self):
        """
        A function that declare all the suffixes neccesaries to access the nl
        file or for the asNMPC
        *It must be call after discretization
        """
        self._model.dual = Suffix(direction=Suffix.IMPORT)
	self._model.ipopt_zL_out = Suffix(direction=Suffix.IMPORT)
	self._model.ipopt_zU_out = Suffix(direction=Suffix.IMPORT)
	self._model.sens_state_0 = Suffix(direction=Suffix.EXPORT) 
	self._model.sens_state_1 = Suffix(direction=Suffix.EXPORT) 
	self._model.sens_state_value_1 = Suffix(direction=Suffix.EXPORT) 
	self._model.sens_init_constr = Suffix(direction=Suffix.EXPORT) 
	self._model.sens_sol_state_1 = Suffix(direction=Suffix.IMPORT) 
	self._model.sens_sol_state_1_z_L = Suffix(direction=Suffix.IMPORT)
	self._model.sens_sol_state_1_z_U = Suffix(direction=Suffix.IMPORT)
	
    def suffixes_definition(self):
        """
        A function that defines all the suffixes for the artificial variables,
        artificial constraints, dual variables ...
        *It must be called after discretization
        """
        suffix_enum = 1
        #Disturbances
        for dis in self.disturbances.itervalues():
            #Creation of artificial constraint for disturbances	
            art_cons_name = dis.name + '_artificial_constraint'
            dis_value = dis.forecast([t for t in self._model_time], dis.rule)
	    def _art_constraint(m, i):	    	
	    	dis[i] = dis_value[i]
	        return dis[i] == dis_value[i]
	    self._model.add_component(art_cons_name, Constraint(self._model_time,
	        rule = _art_constraint))	    
	    art_const = find_comp(art_cons_name, self._model)
            #Creation of suffixes for disturbances
            for indx in dis._index:
                self._model.sens_state_0[dis[indx]] = suffix_enum
		self._model.sens_state_1[dis[indx]] = suffix_enum
		self._model.sens_state_value_1[dis[indx]] = 0
		self._model.sens_init_constr[art_const[indx]] = suffix_enum
		suffix_enum += 1
        #Initial conditions
	for initCond in self.initial_conditions.itervalues():
	    #Creation of artificial constraints for initial conditions
	    art_cons_name = initCond.name + '_artificial_constraint'
	    #Simple Initial condition	    
	    if type(initCond) == SimpleInitialCondition:	    	  
	    	initCond.__class__ = SimpleVar    			   	
	        def _art_constraint(m):	    		    	
	            return initCond == initCond.value
	        self._model.add_component(art_cons_name, 
	        		Constraint(rule = _art_constraint))       	        
	        #Creation of suffixes for initial conditions
	        art_const = find_comp(art_cons_name, self._model)
	        self._model.sens_state_0[initCond] = suffix_enum
		self._model.sens_state_1[initCond] = suffix_enum
		self._model.sens_state_value_1[initCond] = 0
		self._model.sens_init_constr[art_const] = suffix_enum
		suffix_enum += 1
	    #Indexed Initial condition
    	    else:	        
    	    	all_index = initCond._index   	       	    		    		
		def _art_constraint(m, *args):	   		    	
	            return initCond[args] == initCond[args].value
	        self._model.add_component(art_cons_name, Constraint(all_index, 
	        				rule = _art_constraint))	    	            
	        #Creation of suffixes for disturbances       
	        art_const = find_comp(art_cons_name, self._model)
	        for indx in all_index:	        	
		    self._model.sens_state_0[initCond[indx]] = suffix_enum
		    self._model.sens_state_1[initCond[indx]] = suffix_enum
		    self._model.sens_state_value_1[initCond[indx]] = 0
		    self._model.sens_init_constr[art_const[indx]] = suffix_enum
		    suffix_enum += 1		    

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
        #Declare suffixes
        self.suffixes_declaration()        
        #Define suffixes
        self.suffixes_definition()
        #Create the nl file
        iNMPC_nl_map = nl_map_iNMPC(self._model, file_name = file_name)      
        sol_Reader = ResultsReader_sol()                
        #Solver options
        opt_str = 'print_level 0 wantsol 1'
        for op in options:
            opt_str += ' ' + str(op[0]) + ' ' + str(op[1])  	                       
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
	    #Current time
            current_time = period*self._sampling_time
	    #Update initial guess
	    if period == 0:
	        iNMPC_nl_map.update_InitialGuess(self._model, update_duals=False)
            #Update disturbances and initial conditions
            if period > 0:
                self.update_dictUpdateParameeters(current_time)  
                iNMPC_nl_map.update_Parameters(self.dictUpdateParameters)
            #Solution
	    initial_time = time.time()	
	    call([solver_name,file_name+'.nl',opt_str])	    
            #NMPC_results = opt.solve(self._model, tee=False)            
            final_time = time.time()                  
            #Read solution            
            NMPC_results = sol_Reader(file_name+'.sol', None, None, 
            	    suffixes=['dual', 'ipopt_zU_out', 'ipopt_zL_out', 
            	    	'ipopt_zU_in', 'ipopt_zL_in', 'sens_state_0', 
            	    	'sens_state_1', 'sens_state_value_1'])            	
            NMPC_results._smap = iNMPC_nl_map.symbol_map		
            self._model.solutions.load_from(NMPC_results,
            	    				delete_symbol_map=False)
            iNMPC_nl_map.update_InitialGuess(self._model, update_duals=True)
            print '------------------------------------------------------------'
            print 'Simulation period: ' + str(period+1)
            print 'Termination condition: ' +\
                str(NMPC_results.solver.termination_condition)
            print 'Objective function value: ' + str(value(self._model.NEW_OBJ))
            print 'Solution time: ' + str(final_time - initial_time) + ' secs'
            print '------------------------------------------------------------'
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
