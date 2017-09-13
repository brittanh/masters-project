#  _________________________________________________________________________
#
#  Pyomo: Python Optimization Modeling Objects
#  Copyright (c) 2014 Sandia Corporation.
#  Under the terms of Contract DE-AC04-94AL85000 with Sandia Corporation,
#  the U.S. Government retains certain rights in this software.
#  This software is distributed under the BSD License.
#  _________________________________________________________________________

"""
This defines the initial condition for the NMPC
"""

__all__ = ['InitialCondition', 'SimpleInitialCondition', 'IndexedInitialCondition']

from pyomo.core import *
from pyomo.core.base.var import Var, SimpleVar, _VarData, _VarDataWithDomain
from pyomo.dae import *

def InitialCondition(sVar, *arg, **kwds):
    """
    A initial condition for a IVP definition, which is defined as a mutable 
    parameter that is modify recurrently in the (N)MPC scheme      
    """
    # Check validity of variable
    if (not isinstance(sVar,Var)):
        raise TypeError(
            "%s is not a variable. Can only take the derivative " \
            "of a Var component." % (sVar))               
    if (sVar.dim() == 0):
        raise IndexError("The variable %s is not indexed by any " \
       	    "ContinuousSets. A derivative may only be taken with " \
       	    "respect to a continuous domain" % (sVar))
    # Return the element depending in the number of sets            
    if sVar.dim() == 1:        	
        return SimpleInitialCondition(sVar, *arg, **kwds)
    else:
        return IndexedInitialCondition(sVar, *arg, **kwds)                

class SimpleInitialCondition(_VarDataWithDomain, Var):
    """
    It defines a simple IC from SimpleVar
    """
        
    def __init__(self, sVar, *args, **kwds):
    	self.variable = sVar
	self.time_set = kwds.pop('time_set',None)  
	self.time_pos = 0
	self.num_indexes = None
	arg = ()		
        _VarDataWithDomain.__init__(self, self)
        Var.__init__(self, *args, **kwds)
        #self.__class__ = SimpleVar        
	
    def check_InitialCondition(self, time_set):
        """
        It checks all the requirements for the initial condition definition
        """
        if (self.variable is None):
            raise TypeError("No variable has been assigned to the " \
            	    "initial condition '%s'"%(cname(self)))    	            
        if (not isinstance(self.variable,Var) ):#or 
        	#type(self.variable) is ManipulatedVar):
            raise DAE_Error("%s is not a free variable." % (self.variable))            
        try:
            num_contset = len(self.variable._contset)
        except:
            self.variable._contset = {} 
            self.variable._derivative = {}
            if (self.variable.dim() == 0):
                num_contset = 0
            elif (self.variable.dim() == 1):
                sidx_sets = self.variable._index
                if sidx_sets.type() is ContinuousSet:
                    self.variable._contset[sidx_sets] = 0
            else:
                sidx_sets = self.variable._implicit_subsets
                for i,s in enumerate(sidx_sets):
                    if (s.type() is ContinuousSet):
                        self.variable._contset[s] = i
            num_contset = len(self.variable._contset)
        if (num_contset == 0):
            raise DAE_Error("The variable %s is not indexed by any " \
            	    	"ContinuousSets. The variable associated with the " \
            	    	"initial condition must be indexed in a continuous " \
            	    	"domain" % (self.variable))            
        for contset in self.variable._contset:            
            if contset.name != cname(time_set):
                raise TypeError("The variable '%s' for the initial " \
                	"condition '%s' isn't indexed in the continous set "\
            	        "used in the (N)MPC"%(cname(self.variable),self.name))        
        if (self.num_indexes == 2):            
            if self._index not in self.variable._implicit_subsets:
                raise IndexError("The index of the parameter '%s' is not " \
                	"contained in the indexes of the " \
                	"variable '%s'"%(self.name,self.variable))
        elif (self.num_indexes > 2):
            if len(self._implicit_subsets)+1 != self.num_indexes:
                raise IndexError("The number of indexes of the parameter " \
                	"'%s' do not match the number of indexes of the " \
                	"variable '%s'"%(self.name,self.variable))
            for nset in self._implicit_subsets:
                if nset not in self.variable._implicit_subsets:
            	    raise IndexError("The index of the parameter '%s' is not " \
            	    	    "contained in the indexes of the " \
            	    	    "variable '%s'"%(nset.name,self.variable))
            	    
    #def update_InitialCondition(self, sampling_time):
    #    """
    #    This updates the value of the intial condition according to the 
    #    previous solution of the optimal control problem
    #    """                   
    #	# No indexes in the parameter and the variable only indexed in time
    #	if (self.num_indexes is None): # Number of indexes 1
    #	    self[None] = self.variable[sampling_time].value	         	    
	
class IndexedInitialCondition(Var):
    """
    It defines an indexed IC from Var
    """
    def __init__(self, sVar, *args, **kwds):
        """
        Initialization
        """
        self.variable = sVar
	self.time_set = kwds.pop('time_set',None)  
	self.time_pos = 0
	self.num_indexes = None	
        sidx_sets = list(self.variable._implicit_subsets)            
        for i in range(len(sidx_sets)):
            if (sidx_sets[i].name == self.time_set.name):
                self.time_pos = i
        sidx_sets.remove(self.time_set)       
        arg = tuple(sidx_sets)
        try:
	    self.num_indexes = len(self.variable._implicit_subsets)
	except:
	    pass	
	Var.__init__(self, *arg, **kwds)	
	
    def check_InitialCondition(self, time_set):
        """
        It checks all the requirements for the initial condition definition
        """
        if (self.variable is None):
            raise TypeError("No variable has been assigned to the " \
            	    "initial condition '%s'"%(cname(self)))    	            
        if (not isinstance(self.variable,Var) ):#or 
        	#type(self.variable) is ManipulatedVar):
            raise DAE_Error("%s is not a free variable." % (self.variable))            
        try:
            num_contset = len(self.variable._contset)
        except:
            self.variable._contset = {} 
            self.variable._derivative = {}
            if (self.variable.dim() == 0):
                num_contset = 0
            elif (self.variable.dim() == 1):
                sidx_sets = self.variable._index
                if sidx_sets.type() is ContinuousSet:
                    self.variable._contset[sidx_sets] = 0
            else:
                sidx_sets = self.variable._implicit_subsets
                for i,s in enumerate(sidx_sets):
                    if (s.type() is ContinuousSet):
                        self.variable._contset[s] = i
            num_contset = len(self.variable._contset)
        if (num_contset == 0):
            raise DAE_Error("The variable %s is not indexed by any " \
            	    	"ContinuousSets. The variable associated with the " \
            	    	"initial condition must be indexed in a continuous " \
            	    	"domain" % (self.variable))            
        for contset in self.variable._contset:            
            if contset.name != cname(time_set):
                raise TypeError("The variable '%s' for the initial " \
                	"condition '%s' isn't indexed in the continous set "\
            	        "used in the (N)MPC"%(cname(self.variable),self.name))        
        if (self.num_indexes == 2):            
            if self._index not in self.variable._implicit_subsets:
                raise IndexError("The index of the parameter '%s' is not " \
                	"contained in the indexes of the " \
                	"variable '%s'"%(self.name,self.variable))
        elif (self.num_indexes > 2):
            if len(self._implicit_subsets)+1 != self.num_indexes:
                raise IndexError("The number of indexes of the parameter " \
                	"'%s' do not match the number of indexes of the " \
                	"variable '%s'"%(self.name,self.variable))
            for nset in self._implicit_subsets:
                if nset not in self.variable._implicit_subsets:
            	    raise IndexError("The index of the parameter '%s' is not " \
            	    	    "contained in the indexes of the " \
            	    	    "variable '%s'"%(nset.name,self.variable))
            	    
    def var_indx(self, param_index, time_value):
        """
        This function receives the list of indexes of a parameter and sort them
        according to the indexes of a variable
        """
        param_list = list(param_index)    
        param_list.insert(self.time_pos, time_value)
        var_list = param_list                
        return tuple(var_list) 
     
    #def update_InitialCondition(self, sampling_time):
    #    """
    #    This updates the value of the intial condition according to the 
    #    previous solution of the optimal control problem
    #    """                   
    #	# No indexes in the parameter and the variable only indexed in time
    #	if (self.num_indexes is None): # Number of indexes 1
    #	    self[None] = self.variable[sampling_time].value
    #	# One index in the parameter		
    #	elif (self.num_indexes == 2):		
    #	    var_sets = self.variable._implicit_subsets	    	    
    #	    other_set = self._index
    #	    other_set_position = var_sets.index(other_set)	      
    #	    for i in other_set:	    	
    #	    	if other_set_position == 1:
    #	    	    var_idx = sampling_time,i	    	    
    #	    	elif other_set_position == 0:
    #	    	    var_idx = i,sampling_time
    #	        self[i] = self.variable[var_idx].value	        
    #	# Case for more than one index in the parameter
    #	elif (self.num_indexes > 2):
    #	    all_indexes = self._index
    #	    for idx in all_indexes:
    #	        idx_var = self.var_indx(idx, sampling_time)   		        
    #	        self[idx] = self.variable[idx_var].value            	    
