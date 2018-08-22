# Input vector axis of rotation and rotation around it to retrieve 
# rotation matrix with a) exponential map algorith and b) quaternions
# and compare

from rots import qua2mat
from rots import get_rotmat
from rots import get_axial
import numpy as np

#s = np.array([0.128842, 0.412293, 0.901894])
s = np.array([1, 1, 1])*1/np.sqrt(3)
theta = 350

r, A =get_rotmat(s, theta)
M = qua2mat(s,theta)

print('exponential map: \n', r)
print('Determinant of the matrix is: ', np.linalg.det(r))

print('Quaternion: \n', M)
print('Determinant of the matrix is: ', np.linalg.det(M))

#ax = get_axial(r)
#print('get_axial gives: ')
#print(ax[0])
#print(np.round(ax[1], 4))

