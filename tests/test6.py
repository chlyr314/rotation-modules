# In this test the effect of mapping the increment (or variation) of the
# non-additive axial vector to the additive increment of the rotation vector
# so that the current axial vector can be updated by simple addition. Then
# the resulting compound rotation tensor is extracted and its effect is 
# compared with the exact rotation tensor <Rtot> by determining the norm
# of a residual vector r = (R-Rtot)*v, where v is any vector. The results are
# then plotted by varying both dtheta and phi_e (see comments below).
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


# INPUT ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# Determine axial vector
e0 = np.array([0.0,1.0,0.0]) # axis of rotation (unit vector)
theta = 0.1 # current configuration

s0 = (theta*np.pi/180.0)*e0 # axial vector

# Define the measure of variation of e0. By this we determine the 
# incremental axial vector (or the variation of the axial vector) by
# rotating the initial axis of rotation e0 by a small amount, determined 
# by phi_e (See second loop)
dummy = np.array([1.0,0.0,0.0]) # rotate e0 around axis x - dont touch that

# Define a dummy vector (different from e0 though) so that we can have a 
# measure of the accuracy of the Tangent map operator. This vector will be
# The measure of accuracy will be the norm of (R_exact _ R_tangent)*v
v = np.array([0.0, 0.0, 1.0])

# Define a min and a max value for incremental rotation angle (in degrees)
# and a step
min_dtheta = 1
max_dtheta = 90
step = 1

# Determine which levels of variation (in degrees) to plot for second figure
phi_lvls = [10,20,45,90,125]

# END OF INPUT ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

#  rotation angle around incremental axis
dtheta_vals = list(range(min_dtheta,max_dtheta, step))
phi_vals = list(range(2,180,1)) # measure of variation. 

resvals = np.zeros((len(dtheta_vals), len(phi_vals)), dtype=float)

tan_axial = np.zeros((3,len(dtheta_vals), len(phi_vals)), dtype=float)
tot_axial = np.zeros((3,len(dtheta_vals), len(phi_vals)), dtype=float)

# Get rotation matrix corresponding to the axial vector
R0, A0, T = get_rotmat(e0,theta)

for index,dtheta in enumerate(dtheta_vals): # dtheta is the incremental
                                            # rotation angle
    i = 0

    for phi_e in phi_vals:      # <phi_e> is the measure of variation of
                                # the incremental axis of rotation

        # Define an incremental (non-additive) axial vector
        # we define it by rotating e0 by a small ammount around an axis
        R_d, _, _ = get_rotmat(dummy,phi_e) 

        # Get variation of axis of rotation 
        de = np.matmul(R_d,e0)

        # Incremental axial vector (non-additive)
        ds = (dtheta*np.pi/180.0)*de

        DR, _, _ = get_rotmat(de, dtheta) # incremental rotation tensor

        # Get exact total rotation tensor and corresponding compound
        # axial vector
        Rtot = np.matmul(DR,R0)
        e_tot,theta_tot = get_axial(Rtot)
        tot_ax = (theta_tot*np.pi/180.0)*e_tot # "exact" axial vector

        # s0 + ds will not give the axial vector that corresponds to Rtot,
        # because ds is not additive. We need to apply a transformation 
        # in order to get the additive increment to the axial vector, say a_ds
        tan_ax, R, _ = updateAxialVector(s0,ds,T) # approximate axial vector

        # See what the norm of a residual vector does by increasing phi_e or 
        # dtheta. If tangent map was exact, then res_norm should --> 0
        r = np.matmul(R-Rtot,v)
        res_norm = np.linalg.norm(r)
        resvals[index,i] = res_norm

        # Store components of "exact" and approximated total axial vector
        tan_axial[0,index,i],tan_axial[1,index,i],tan_axial[2,index,i] = \
                tan_ax[0],tan_ax[1],tan_ax[2]
        tot_axial[0,index,i],tot_axial[1,index,i],tot_axial[2,index,i] = \
                tot_ax[0],tot_ax[1],tot_ax[2]

        i = i + 1

# Plot the phi_e vs res_norm curves for a range of incremental thetas
# This plot shows how the residual norm changes with respect to the variation
# of the incremental axis of rotation, keeping the incremental rotation around
# it constant. Plot for several dtheta values.
plt.figure(1)
for j in range(len(dtheta_vals)):
    plt.plot(phi_vals, resvals[j][:])

plt.xlabel('Variation measure (in angles)')
plt.ylabel('Norm of residual vector')
subtit = '\nRange of del_theta:[{0}-{1},{2}]'.format(min_dtheta,max_dtheta,step)
plt.title('Accuracy of tangent map of incremental axial vector'+ subtit)


# Plot norm of residual at different phi_e levels, for all dthetas
# This plot shows how the residual norm changes with respect to increasing
# levels of the incremental rotation dtheta, keeping the incremental axis of
# rotation constant. Plot for 4 different measures of the incremental axis
plt.figure(2)
if (len(dtheta_vals) != 1):
    ind = [ind_val for ind_val,k in  enumerate(phi_vals) if k in phi_lvls]
    plt.plot(dtheta_vals, resvals[:,ind])

    plt.xlabel('Incremental rotation angle(dtheta)')
    plt.ylabel('Norm of residual vector')
    plt.title('Change in accuracy at different variation levels')
    labels = [str(val)+' degrees' for val in phi_lvls]
    plt.legend(labels)

# Plot axial vector max error components vs variation of axial vector measure
# (tot_axial_i - tan_axial_i) vs phi_e. For each phi_e value, we calculate at
# which dtheta the error in the axial vector components is max and then plot
# the max_error vs phi_e for each component (3 plots)
maxPhiX = list(range(len(phi_vals)))
maxPhiY = list(range(len(phi_vals)))
maxPhiZ = list(range(len(phi_vals)))
for m in range(len(phi_vals)):
    maxPhiX[m] = np.absolute(np.amax(tot_axial[0,:,m]-tan_axial[0,:,m]))
    maxPhiY[m] = np.absolute(np.amax(tot_axial[1,:,m]-tan_axial[1,:,m]))
    maxPhiZ[m] = np.absolute(np.amax(tot_axial[2,:,m]-tan_axial[2,:,m]))

fig, ax = plt.subplots(3, 1, sharex=True)
fig.subplots_adjust(hspace=0.5)
ax[0].plot(phi_vals,maxPhiX)
ax[1].plot(phi_vals,maxPhiY)
ax[2].plot(phi_vals,maxPhiZ)
ax[0].set(title='max_err vs variation of incremental axial vector',ylabel='θx')
ax[1].set(ylabel='θy')
ax[2].set(ylabel='θz')
ax[0].yaxis.tick_right()
ax[1].yaxis.tick_right()
ax[2].yaxis.tick_right()

# Plot axial vector max error components vs incremental angle of rotation
# (dtheta). For each dtheta, we calculate at which phi_e the error in the 
# axial vector components is max and then plot the max_error vs dtheta 
# for each component (3 plots)
if (len(dtheta_vals) != 1):
    maxdThetaX = list(range(len(dtheta_vals)))
    maxdThetaY = list(range(len(dtheta_vals)))
    maxdThetaZ = list(range(len(dtheta_vals)))
    for m in range(len(dtheta_vals)):
        maxdThetaX[m] = np.absolute(np.amax(tot_axial[0,m,:]-tan_axial[0,m,:]))
        maxdThetaY[m] = np.absolute(np.amax(tot_axial[1,m,:]-tan_axial[1,m,:]))
        maxdThetaZ[m] = np.absolute(np.amax(tot_axial[2,m,:]-tan_axial[2,m,:]))

    fig2, ax2 = plt.subplots(3, 1, sharex=True)
    fig.subplots_adjust(hspace=0.5)
    ax2[0].plot(dtheta_vals,maxdThetaX)
    ax2[1].plot(dtheta_vals,maxdThetaY)
    ax2[2].plot(dtheta_vals,maxdThetaZ)
    ax2[0].set(title='max_err vs incremental angle of rotation',ylabel='θx')
    ax2[1].set(ylabel='θy')
    ax2[2].set(ylabel='θz')
    ax2[0].yaxis.tick_right()
    ax2[1].yaxis.tick_right()
    ax2[2].yaxis.tick_right()

plt.show(1)
plt.show(fig)
plt.show(fig2)
