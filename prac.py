from rots import get_rotmat
from rots import get_axial
import numpy as np

#s = np.array([0.128842, 0.412293, 0.901894])
s = np.array([1, 1, 1])*1/np.sqrt(3)
theta = 0.1
r, A =get_rotmat(s, theta)
print(r)
print('Determinant of the matrix is: ', np.linalg.det(r))
print('initial vector s: ')
print(s)
print(theta)

ax = get_axial(r)
print('get_axial gives: ')
print(ax[0])
print(np.round(ax[1], 4))

