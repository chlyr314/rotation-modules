# Second test file: Given an axis of rotation <s>, an angle of 
# rotation <theta> and an initial vector <v0>, it plots the path
# of the rotation. Can also plot 1st and second order approximations 
# of the rotation.

import numpy as np
from rotplot import plotrot as pltr
import matplotlib.pyplot as plt

# Axis of rotation
s = np.array([1, 1, 1])*1/np.sqrt(3)
#s = np.array([0, 1, 0])

# Compound rotation around s (in degrees)
theta = 360

# Set an initial vector
v0 = np.array([1, 0, 1])*1/np.sqrt(2)

# set number of steps
n=30

# initialize figure handle for 3d plot
fig = plt.figure()
ax = fig.add_subplot(111, projection='3d')

# plot
v = pltr(s, theta, v0, n, ax)
#######
theta = 90
col = 'k'
dflag = True
n=20
v = pltr(s, theta, v0, n, ax, col, dflag)
#######
plt.show(ax)
