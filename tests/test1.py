# Input a unit vector representing the axis of rotation and the 
# range of angles of rotation around it. Loop over the range and 
# plot a) the norm of the residual vector for the whole range, where 
# the residual vector is given as r = (R-Ra)*v and b) the components of
# the exact axial vector that we give as input vs the components of the 
# axial vector as is retrieved from the get_axial() method. Singularity should
# be demonstrated in the retrieved axial vector compronents when Ra is close
# to 180 degrees. (total lagrangian formulation problem etc..)
# TODO: add quaternion extraction for comparison
import sys
import os
import numpy as np
from matplotlib import pyplot as plt

pwd = os.getcwd()
sys.path.append(pwd + '/../utils')

from rots import *

# Define any unit vector by specifying first 2 components
nx, ny = 0.3, 0.4
nz = np.sqrt(1.0-nx*nx-ny*ny)

# Define a unit vector representing the axis of rotation and 
# a rotation angle around it (in degrees)
e = np.array([1, 1, 1])*1/np.sqrt(3)

# Define range of angle of rotations
theta_min = 1.0
theta_max = 720.0
dtheta = 1.0

# END OF INPUT ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# Generate thetas
numpoints = (theta_max-theta_min)/dtheta + 1
thetas = np.linspace(theta_min, theta_max, numpoints)

num = len(thetas)

axials = np.zeros((3, num)) # array that stores exact axial
R_axials = np.zeros((3, num)) # array that stores axials retrieved from R

# define a vector and a container to store residual norm between
# R*v and R_a*v. Should be zero!
v = np.array([1,3,4])
residual = np.zeros(num)

for i, theta in enumerate(thetas):

    # get axial vector
    s = e*np.pi*theta/180.0

    # get rotation matrix
    R, _, _ = get_rotmat(s)

    # reverse now: get axial vector from rotation matrix
    s_dot = get_axial(R,1)

    # get rotation matrix from retrieved axial vector
    Ra,_,_ = get_rotmat(s_dot)

    # store exact and retrieved axial vectors 
    axials[:,i] = s
    R_axials[:,i] = s_dot

    # store norm of residual vector
    res = np.linalg.norm(np.matmul(R-Ra,v)) # should be zero!
    residual[i] = res

# Plot residual norm for residual vector r = (R-Ra)*v.
# Should be zero all the way
plt.figure(1)
plt.plot(thetas,residual)

# Plot components of the axial vector as we give it and as it is retrieved
# from the get_axial() algorithm. Singularity in component form should occur
# every (2*k*pi + pi).
fig, ax = plt.subplots(3, 1, sharex = True)
ids = ['$\\theta_x$','$\\theta_y$','$\\theta_z$']
for j in range(3):
    ax[j].plot(thetas,axials[j,:], thetas, R_axials[j,:], 'b--', linewidth=1)
    ax[j].set(ylabel=ids[j])

fig.suptitle('Components of exact and retrieved axial vector\n$e = \
        [{0}, {1}, {2}]$'.format(np.round(e[0],4),np.round(e[1],4),np.round(e[2],4)))
plt.xlabel('Angle of rotation(degrees)')
plt.xticks(np.arange(0, theta_max+60, step=60))

plt.show(fig)
plt.show(1)

################################
#q, M = vec2qua(s,theta)
#
#print('exponential map: \n', r)
#print('Determinant of the matrix is: ', np.linalg.det(r))
#
#print('Quaternion: \n', M)
#print('Determinant of the matrix is: ', np.linalg.det(M))


