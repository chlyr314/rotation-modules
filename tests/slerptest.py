import sys
import os

pwd = os.getcwd()
sys.path.append(pwd + '/../utils')

import numpy as np
from rots import *
from rotplot import *
import matplotlib.pyplot as plt

s1 = np.array([1,0,0])
theta1 = 20

s2 = np.array([1,1,1])/np.sqrt(3)
theta2 = 120

q1, R1 = vec2qua(s1,theta1)
q2, R2 = vec2qua(s2,theta2)

fig = plt.figure()
ax = fig.add_subplot(111,projection='3d')

ref = [0, 5, 5]

for i in range(3):

    v = R1[:,i]
    pvec(v,ax,'k',ref)

for t in range(5):

    ti = (t+1)/5
    ref = [t+1, 5, 5]
    q = slerp(q1,q2,ti)
    R = qua2mat(q)
    print(np.linalg.det(R))
    for i in range(3):

        v = R[:,i]
        pvec(v,ax,'r',ref)

#r = 10
#ax.set_xlim(0,r)
#ax.set_ylim(0,r)
#ax.set_zlim(0,r)
plt.show()

