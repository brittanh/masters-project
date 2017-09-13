#  _________________________________________________________________________
#
#  Pyomo: Python Optimization Modeling Objects
#  Copyright (c) 2014 Sandia Corporation.
#  Under the terms of Contract DE-AC04-94AL85000 with Sandia Corporation,
#  the U.S. Government retains certain rights in this software.
#  This software is distributed under the BSD License.
#  _________________________________________________________________________

"""
Dynamic plotting of a set of variables
"""

__all__ = ['DynamicPlot']

import matplotlib.pyplot as plt

class DynamicPlot():
    """
    This is a class for dynamic plotting
    """    
    def __init__(self, name, x_label, y_label, legend, total_variables):
	plt.ion()    	        
        self.name = name
        self.legend = legend
        self.total_variables = total_variables
        self.fig, self.ax = plt.subplots()
        plt.xlabel(x_label)
        plt.ylabel(y_label)
        plt.title(self.name)        
        self.lines = []
        for i in range(total_variables):
            ld, = self.ax.plot([], [])
            self.lines.append(ld)      
        plt.legend(self.legend)
        self.ax.set_autoscaley_on(True)
    
    def update_plot(self, xdata_input, ydata_input):
        """
        This function update the plot based on new values
        """
        for i in range(self.total_variables):
            self.lines[i].set_xdata(xdata_input)		
            self.lines[i].set_ydata(ydata_input[i])
	self.ax.set_xlim(min(xdata_input), max(xdata_input))
	self.ax.relim()
	self.ax.autoscale_view()		
	self.fig.canvas.draw()
	self.fig.canvas.flush_events()
	
    def display_final_plot(self):
    	"""
    	It desactives the interactive plotting and shows the final plot 
    	"""
    	plt.ioff()
    	plt.show()
