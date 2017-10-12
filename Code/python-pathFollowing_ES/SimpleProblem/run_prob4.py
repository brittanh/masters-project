'''
Created on 18. nov. 2015

@author: suwartad
'''

# -*- coding: utf-8 -*-
import numpy as np
import matplotlib.pyplot as plt
from jpredictor import *
from prob4 import prob4, objective

dual    = np.array([[0.735758882342878], [23.729329433526743]])
xstart  = np.array([[1.000000000000000], [0.367879441171442]])
primal,dual,t,x1,x2 = jpredictor(prob4, objective, np.array([[8.0],[1.0]]), np.array([[1.0],[-4.0]]), xstart, dual, 0.1, np.array([]), np.array([]), 0)


xsolution = np.array([x1[-1],x2[-1]])
print "primal solution: " 
print xsolution

# ploting
plt.subplot(2,1,1)
plt.plot(t,x1,'--rs');
plt.xlabel('t (parameter)'); 
plt.ylabel('x1');
plt.title('x1 vs. t');

plt.subplot(2,1,2)
plt.plot(t,x2,'--rs');
plt.xlabel('t (parameter)');
plt.ylabel('x2');
plt.title('x2 vs. t');

plt.show()

'''
# ploting 
numInfo = numel(info);
t  = np.zeros([1,numInfo]);
x1 = np.zeros([1,numInfo]);
x2 = np.zeros([1,numInfo]);

# need to change this part!
# STRUCT and Dictionary in Python baca ini http://docs.scipy.org/doc/scipy/reference/tutorial/io.html
for i=1:numInfo
    t(1,i)  = info(i).t;
    x1(1,i) = info(i).x(1,:);
    x2(1,i) = info(i).x(2,:);
end

#xsolution = [x1(end) x2(end)];
xsolution = [x1[-1] x2([-1]];

C(:,1) = {'LineWidth'; 2};
C(:,2) = {'MarkerEdgeColor'; 'k'};
C(:,3) = {'MarkerFaceColor'; 'g'};


#figure(1)
plt.subplot(2,1,1)
plt.plot(t,x1,'--rs', C{:});
plt.xlabel('t (parameter)'); 
plt.ylabel('x1');
plt.title('x1 vs. t');

#figure(2);
plt.subplot(2,1,2)
plt.plot(t,x2,'--rs', C{:});
plt.xlabel('t (parameter)');
plt.ylabel('x2');
plt.title('x2 vs. t');
'''
