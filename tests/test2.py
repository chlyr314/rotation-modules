# Second test file: Given an axis of rotation <s>, an angle of 
# rotation <theta> and an initial vector <v0>, it plots the path
# of the rotation. Can also plot 1st and second order approximations 
# of the rotation.
import sys
import os

pwd = os.getcwd()
sys.path.append(pwd + '/../utils')

from mpl_toolkits.mplot3d import Axes3D
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
n=20

# initialize figure handle for 3d plot
fig = plt.figure()
ax = fig.add_subplot(111, projection='3d')

# plot
v = pltr(s, theta, v0, n, ax)

# Plot linear approx
theta = 200
col = 'k'
dflag = True
d2flag = True
n=10
v = pltr(s, theta, v0, n, ax, col, dflag,d2flag)

# plot sphere
u = np.linspace(0, 2 * np.pi, 100)
v = np.linspace(0, np.pi, 100)
x = np.outer(np.cos(u), np.sin(v))
y = np.outer(np.sin(u), np.sin(v))
z = np.outer(np.ones(np.size(u)), np.cos(v))

# Plot the surface
ax.plot_surface(x, y, z, color='b',alpha = 0.3)

plt.title('n=10,angle=60', fontsize = 20)
plt.show(ax)
