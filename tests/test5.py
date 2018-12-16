# Routine that plots path of vector v0 as it rotates around axis <s>.
# You can make it so that <s> also rotates around another fixed vector
# at the same time
# TODO: add transformation from non-additive to additive components!
import sys
import os

pwd = os.getcwd()
sys.path.append(pwd + '/../utils')

from rots import *
import numpy as np
from rotplot import plotrot as pltr
from rotplot import pvec
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D

# INPUT =================

# vector to be rotated
v0 = np.array([1,0,0])

# axis of rotation - it will also rotate
s = np.array([1, 1, 1])*1/np.sqrt(3)

total_rotation = 570.0     # total rotation around s

# rotation of axis of rotation around z
phi = 40;
e3 = np.array([0,0,1])
k = 30 # number of increments
iters = [j for j in range(1,k+1)] # generate sequence with number of iterations

# END OF INPUT ============

dphi = phi/k   # increment for rotation of s around e3
Rs,a = get_rotmat(e3,dphi)

dtheta = total_rotation/k  # increment for rotation around s

fig = plt.figure()
ax = fig.add_subplot(111, projection='3d')
pvec(v0,ax,'k',[0,0,0],'--')
R = np.identity(3)

history = np.zeros((3,k))

for i in range(1,k+1):

    v0 = pltr(s, dtheta, v0, 1, ax, d2flag=False)
    dRv,dAv = get_rotmat(s,dtheta)

    R = np.matmul(dRv,R)   # compound rotation

    s = np.matmul(Rs,s)   # rotate axis of rotation
    compound_axis,compound_theta = get_axial(R)

    # compound axial vector extracted from TOTAL rotation vector
    # we expect discontinuities
    compound_axis = compound_axis*np.pi*compound_theta/180.0 #total axial vector
    history[0,i-1] = compound_axis[0]
    history[1,i-1] = compound_axis[1]
    history[2,i-1] = compound_axis[2]



# plot

pvec(v0,ax,'k',[0,0,0],'--') # plot rotated vector

# plot sphere
u = np.linspace(0, 2 * np.pi, 25)
v = np.linspace(0, np.pi, 25)
x = np.outer(np.cos(u), np.sin(v))
y = np.outer(np.sin(u), np.sin(v))
z = np.outer(np.ones(np.size(u)), np.cos(v))

# Plot the surface
ax.plot_surface(x, y, z, color='b',alpha = 0.3)

# Plot history of components of rotation vector
fig2 = plt.figure()

theta_x = plt.plot(dtheta*np.array(iters),history[0,:],'b-',label = 'theta_x')
theta_y = plt.plot(dtheta*np.array(iters),history[1,:],'r--',label = 'theta_y')
theta_z = plt.plot(dtheta*np.array(iters),history[2,:],'k-.',label = 'theta_z')
plt.axhline(y=0,color='k',ls='-',lw=0.5)
plt.legend()
plt.xlabel('Rotation angle (degrees)', fontsize = 18)
plt.ylabel('Components of axial vector', fontsize = 14)
plt.title('Rotation about follower axis, theta =570', fontsize = 20)
plt.show()
