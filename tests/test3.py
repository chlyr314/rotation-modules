# test 3: Given an ARRAY of rotation vectors and the
# corresponding array of rotations around them, plot the 
# rotation path of a vector v0
import sys
import os

pwd = os.getcwd()
sys.path.append(pwd + '/../utils')

import numpy as np
from rotplot import plotrot as pltr
import matplotlib.pyplot as plt

# get unit vectors
s1 = np.array([0,1,0])
s2 = np.array([0,0,1])
s3 = np.array([1,0,0])

a = np.asarray([s1, s2, s3])

# define rotations about unit vectors (in degrees)
theta1 = 90
theta2 = 90
theta3 = 90

theta = [theta1,theta2, theta3]

# Set an initial vector
v0 = np.array([0,0,1])*1/np.sqrt(1)

# Discretize theta in n part
n = 20

# initialize figure handle for 3d plot
fig = plt.figure()
ax = fig.add_subplot(111,projection='3d')

for i,th in enumerate(theta):
    s = a[i]
    v0 = pltr(s,th,v0,n,ax)


plt.show(ax)
