#  _________________________________________________________________________
#
#  Pyomo: Python Optimization Modeling Objects
#  Copyright (c) 2014 Sandia Corporation.
#  Under the terms of Contract DE-AC04-94AL85000 with Sandia Corporation,
#  the U.S. Government retains certain rights in this software.
#  This software is distributed under the BSD License.
#  _________________________________________________________________________

"""
This defines the manipulated variables for the NMPC
"""

__all__ = ['ManipulatedVar']

from pyomo.core.base.var import Var, _VarData
from pyomo.core import *
from pyomo.dae import *

class ManipulatedVar(Var):
    """
    A manipulated variable definition, which is defined as a variable
    """
    def __init__(self, *arg, **kwds):    	
        Var.__init__(self, *arg, **kwds)        
        self.num_args = len(arg)        
        
    def check_ManipulatedVar(self, time_set):
        """
        It checks all the requirement for a Manipulated variable
        """
        if self.num_args > 1:
            raise ValueError("The ManipulatedVar '%s' must " \
            	    "receive only one continuous set as argument"%(cname(self)))
        if self.num_args == 0:
            raise ValueError("The ManipulatedVar '%s' must receive " \
            	    "one continuous set as argument"%(cname(self)))        
        if type(self._index) is not ContinuousSet:
            raise ValueError("The argument for the ManipulatedVar '%s' must "\
            	    "be a continuos set"%(cname(self)))            
        if cname(self._index) != cname(time_set):
            raise TypeError("The set defined for the ManipulatedVar " \
            	    "'%s' does not match with the continous time set used " \
            	    "in the (N)MPC"%(cname(self)))
        if self._implicit_subsets is None:
            self.time_pos = 0
        else:
            sidx_sets = list(self._implicit_subsets)       
            for i in range(len(sidx_sets)):
                if (sidx_sets[i].name == time_set.name):
                    self.time_pos = i    
            
    def var_indx(self, param_index, time_value):
        """
        This function receives the list of indexes of a parameter and sort them
        according to the indexes of a variable
        """
        param_list = list(param_index)    
        param_list[self.time_pos] = time_value                      
        return tuple(param_list)            
