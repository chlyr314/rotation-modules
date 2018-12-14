# This file showcases that the compount axis of rotation for 2 sequential
# rotations is not the vector addition of the axes for the first and second
# rotation. This can be seen from the fact that it is not in their span
import sys
import os

pwd = os.getcwd()
sys.path.append(pwd + '/../utils')

import numpy as np
from rots import get_rotmat
from rotplot import plotrot as pltr
import matplotlib.pyplot as plt


# Plot incremental rotation

# Axis of rotation
s = np.array([0, 1, 0])

# Compound rotation around s (in degrees)
theta = 45

# Set an initial vector
v0 = np.array([1, 0, 1])*1/np.sqrt(2)

# Discretize theta in n part
n = 15

# Get rotation matrix and axial matrix for dtheta
dtheta = theta/n

fig = plt.figure()
ax = fig.add_subplot(111, projection='3d')

# plot
v0 = pltr(s, theta, v0, n, ax, d2flag=False)

# Second rotation
s = np.array([0, 0, 1])
v0 = pltr(s, theta, v0, n, ax, d2flag=False)

# rotate initial vector v0 to v and compare results
sz = np.array([1, 0, 1])*1/np.sqrt(2)
v1 = np.array([1, 1, 0])*1/np.sqrt(2)

# compound axis of rot
t = np.cross(v0, v1)
s = t/np.linalg.norm(t)

# get angle between vectors v0-v1
sintheta = np.dot(v0, v1)/(np.linalg.norm(v0)*np.linalg.norm(v1))
theta = 2*np.arcsin(sintheta)*180/np.pi
dtheta = theta/n
r, S = get_rotmat(s, dtheta)

v0 = pltr(s, theta, v0, n, ax, col='b')

plt.show()
