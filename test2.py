# Second test file: Given an axis of rotation <s>, an angle of 
# rotation <theta> and an initial vector <v0>, it plots the path
# of the rotation. Can also plot 1st and second order approximations 
# of the rotation.

import numpy as np
from rots import get_rotmat
from rots import get_axial
from rotplot import plotrot as pltr
import matplotlib.pyplot as plt

# Axis of rotation
#s = np.array([1,1,1])*1/np.sqrt(3)
s = np.array([0,1,0])

# Compound rotation around s (in degrees)
theta = 90

# Set an initial vector
v0 = np.array([1,0,1])*1/np.sqrt(2)

# Discretize theta in n part
n = 15

# Get rotation matrix and axial matrix for dtheta
dtheta = theta/n
r,S= get_rotmat(s,dtheta)

fig = plt.figure()
ax = fig.add_subplot(111,projection='3d')

# plot
v0 = pltr(r,v0,n,S,ax,d1flag=False,d2flag=False)

plt.show(ax)
