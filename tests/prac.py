# Input vector axis of rotation and rotation around it to retrieve 
# rotation matrix with a) exponential map algorithm and b) quaternions
# and compare
import sys
import os
import numpy as np

pwd = os.getcwd()
sys.path.append(pwd + '/../utils')

from rots import *

s = np.array([1, 1, 1])*1/np.sqrt(3)
theta = 70

r, A =get_rotmat(s, theta)

q, M = vec2qua(s,theta)

print('exponential map: \n', r)
print('Determinant of the matrix is: ', np.linalg.det(r))

print('Quaternion: \n', M)
print('Determinant of the matrix is: ', np.linalg.det(M))

#ax = get_axial(r)
#print('get_axial gives: ')
#print(ax[0])
#print(np.round(ax[1], 4))

