# masters-project
Brittany Hall's Specialization Project Autumn 2017

Based on work done by Eka Suwartadi (https://github.com/detu/licq-path-following).

Used his work and converted it to use Python.
The code requires:
** Python
** Pyomo
** 

This code is an implementation of economic NMPC scheme.
The workflow is:
1. steady-state optimization (gives the set-point)
2. Dynamic optimization (MPC controller)

The steady-state code is located at \models\ folder and is named distACstrSS.py

