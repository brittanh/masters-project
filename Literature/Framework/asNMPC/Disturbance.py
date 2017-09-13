#  _________________________________________________________________________
#
#  Pyomo: Python Optimization Modeling Objects
#  Copyright (c) 2014 Sandia Corporation.
#  Under the terms of Contract DE-AC04-94AL85000 with Sandia Corporation,
#  the U.S. Government retains certain rights in this software.
#  This software is distributed under the BSD License.
#  _________________________________________________________________________


"""
This defines the class of disturbances in the NMPC
"""

__all__ = ['Disturbance']

from pyomo.environ import *
from pyomo.core.base.var import Var, _VarData
from pyomo.dae import *

def _constant_forecast(time_set, forecast_function):
    """
    This is the function update when the forecast is constant
    """    
    ans = {}    
    for t_indx in time_set:    	
        ans[t_indx] = forecast_function(min(time_set))             
    return ans

def _perfect_forecast(time_set, forecast_function):
    """
    This is the function update when the forecast is constant
    """
    ans = {}
    for t_indx in time_set:
        ans[t_indx] = forecast_function(t_indx)        
    return ans

class Disturbance(Var):
    """
    A disturbance definition, which is defined as a mutable parameter that is 
    modify recurrently in the (N)MPC scheme      
    """
    def __init__(self, *args, **kwds):
    	# All forecasting types
    	self.all_forecasts = {
    		'CONSTANT': _constant_forecast,
    		'PERFECT': _perfect_forecast}    	
    	# Specific kwd for this class
    	tmpforecast = kwds.pop('forecast','CONSTANT')
    	self.forecast_name = tmpforecast.upper()    	
    	self.forecast = self.all_forecasts.get(self.forecast_name,None)
    	if (self.forecast is None):
    	    raise ValueError("Unknown forecasting method '%s' specified using" \
    	    	    " the 'forecast' keyword. Valid forecasting methods" \
    	    	     "are 'CONSTANT' or 'PERFECT'"%(tmpforecast))    	
    	self.rule = kwds.pop('rule',None)    	   
        self.num_args = len(args)    	
    	# Code copy from the __init__ method of the Var class    	
        Var.__init__(self, *args, **kwds)
                                   
           
    def check_disturbance(self, time_set):
        """
        It checks all the requirements for the disturbance definition
        """
        if (self.rule is None):
            raise TypeError("No rule has been defined for the disturbance '%s'"\
            	    %(cname(self)))        
        if (self.num_args > 1):
            raise ValueError("The disturbance '%s' must receive only one" \
            	    " continuous set as argument"%(cname(self)))
        if (self.num_args == 0):
            raise ValueError("The disturbance '%s' must receive one " \
            	    "continuous set as argument"%(cname(self)))        
        if (type(self._index) is not ContinuousSet):
            raise ValueError("The argument for the disturbance '%s' must "\
            	    "be a continuos set"%(cname(self)))            
        if (cname(self._index) != cname(time_set)):
            raise TypeError("The set defined for the disturbance '%s' does "\
            	    "not match with the continous time set used in" \
            	    "the (N)MPC"%(cname(self)))
        
    #def update_disturbance(self, current_time):
    #	"""
    # 	This updates the value of the disturbance according to the function
    #	defined by the modeler and the type of forecasting used    	
    #	"""                    
    #    # Time representation
    #    dt = self._index
    #    time = [current_time + t for t in dt]        
    #    function = self.rule        
    #	# Forecast                  
    #	disturbance_values = self.forecast(time, function)        
    #    # Assigment
    #    for t in dt:            
    #        self[t] = disturbance_values[t+current_time]              
            
    def get_realValue(self, current_time):
        """
    	This returns the real value of the disturbance for the current time  	
    	"""      
    	dt = self._index
        time = [current_time + t for t in dt]   
    	disturbance_values = _perfect_forecast(time, self.rule)
    	return disturbance_values[current_time]    	              	                   
