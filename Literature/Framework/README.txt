This framework has only been tested using Pyomo version 4.1.10519. 
Copy the three folders (asNMPC, iNMPC, and iNMPC_nl_interface) to Pyomo's directory (e.g. /usr/local/lib/python2.7/dist-packages/pyomo/)
For this framework to work properly some of the Pyomo's core files have to be modify:

	Directory			File		Modifications
	.../pyomo/dae			misc.py		*In the function update_conset_indexed_component add the method update_expressions
							*Add the function update_expression
							*Add the function get_index_information
							*Add the function _get_idx
	.../pyomo/dae/plugins		colloc.py	*Add the function reduce_collocation_points
							*Add the function _interpolation_coeffs
							*At the begining of the script import get_index_information from pyomo.dae.misc 
	../pyomo/environ		__init__.py	*Add the names iNMPC, asNMPC, and iNMPC_nl_interface to the corresponding list to load these modules


The detailed modifications can be found in the file "PyomoModifications_V4.1.10519"