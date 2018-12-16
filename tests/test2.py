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
from rotplot import minRotMat
import matplotlib.pyplot as plt


# Axis of rotation
s1 = np.array([0, 1, 0])

# Second rotation
s2 = np.array([0, 0, 1])

# Rotation angle (in degrees) around s1 and s2
theta = 45

# Set an initial vector
v0 = np.array([1, 0, 1])*1/np.sqrt(2)

# Decompose theta in n (~small) rotations
n = 15

fig = plt.figure()
ax = fig.add_subplot(111, projection='3d')

# plot
v1 = pltr(s1, theta, v0, n, ax, d2flag=False)

v2 = pltr(s2, theta, v1, n, ax, d2flag=False)

# get compound axis of rotation, rotation angle and rotation matrix
s3,theta_c,R = minRotMat(v0,v2)

v0 = pltr(s3, theta_c, v0, n, ax, col='b')

plt.show()
