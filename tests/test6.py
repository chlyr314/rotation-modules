# In this test the effect of mapping the increment (or variation) of the
# non-additive axial vector to the additive increment of the rotation vector
# so that the current axial vector can be updated by simple addition. Then
# the resulting compound rotation tensor is extracted and its effect is 
# compared with the exact rotation tensor <Rtot> by determining the norm
# of a residual vector r = (R-Rtot)*v, where v is any vector. The results are
# then plotted by varying both dtheta and phi_e (see comments below).
# TODO : add features to plot and plot total axial vector components
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
e0 = np.array([0.0,1.0,0.0])
theta = 45.0

s0 = (theta*np.pi/180.0)*e0 # axial vector

# incremental rotation angle
dtheta = 20.0

# Define the measure of variation of e0. By this we determine the 
# incremental axial vector (or the variation of the axial vector) by
# rotating the initial axis of rotation e0 by a small amount, determined 
# by phi_e
phi_e = 5.0   # degrees
dummy = np.array([1.0,0.0,0.0]) # rotate e0 around axis x - dont touch that


# END OF INPUT ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


dtheta_vals = list(range(2,92,2)) #  rotation angle around incremental axis
phi_vals = list(range(2,400,2)) # measure of variation. 


for dtheta in dtheta_vals:

    resvals = np.zeros(len(phi_vals), dtype=float)
    i = 0

    for phi_e in phi_vals:

        # Get rotation matrix corresponding to the axial vector
        R0, A0, T = get_rotmat(e0,theta)

        # Define an incremental (non-additive) axial vector
        # we define it by rotating e0 by a small ammount around an axis
        R_d, _, _ = get_rotmat(dummy,phi_e)

        # Get variation of axis of rotation
        de = np.matmul(R_d,e0)

        # Incremental axial vector (non-additive)
        ds = (dtheta*np.pi/180.0)*de

        DR, _, _ = get_rotmat(de, dtheta) # incremental rotation tensor

        # Get exact total rotation tensor
        Rtot = np.matmul(DR,R0)

        # s0 + ds will not give the axial vector that corresponds to Rtot, because
        # ds is not additive. We need to apply a transformation in order to get the
        # additive increment to the axial vector, say a_ds

        _, R, _ = updateAxialVector(s0,ds,T)

        # See what the norm of a residual vector does by increasing phi_e or dtheta
        v = np.array([0.0, 0.0, 1.0])
        r = np.matmul(R-Rtot,v)
        res_norm = np.linalg.norm(r)
        resvals[i] = res_norm
        i = i +1

    plt.plot(phi_vals, resvals)

plt.show()
