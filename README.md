# masters-project
## Brittany Hall's Specialization Project Autumn 2017 ##

### Norwegian University of Science and Technology (Trondheim, Norway)

Based on work done by Eka Suwartadi: [Implementation of path-following algorithm for LICQ case study](https://github.com/detu/licq-path-following).
Eka's work was utilised throughout this project to create an equivalent but more universal algorithm in Python.

To use the code the following elements are required to be installed:

* Python 2.7
* [numpy](http://www.numpy.org/)
* [casadi](https://github.com/casadi/casadi/wiki) (version v3.2.0)
* 
## Purpose
This code implements economic NMPC using a path-following algorithm.
## Getting Started
The scripts work as follows:

* Steady state optimization 
	* Edit `params.py` to have desired column and CSTR parameters
	* Run `Col_CSTR_SS.py`
		* Requires `params.py` and `nlp_solve.py`
		* Creates files needed for dynamic optimization:
			* `CstrDistXinit.mat` : contains steady state optimal inputs
			* `LambdaCstrDist.mat` : contains steady state optimal Lagrange multipliers 
			*  `Qmax.mat` : Greshgorin convexification objective function weights
* Dynamic Optimization
## Distillation Column Model ##

The distillation column model used was developed by Sigurd Skogestad and is referred to as "Column A".

Details of the model and the original files can be found on his [website](http://folk.ntnu.no/skoge/book/matlab_m/cola/cola.html).

## CSTR Model ##

Simple first order reaction A ->B 