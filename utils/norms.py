import numpy as np
import time
import sys
from rots import get_rotmat
from rots import get_axial

def skew(a):
        A = np.array([[0,-a[2],a[1]],[a[2],0,-a[0]],[-a[1],a[0],0]])
        return A
def getparams(x,theta,s):

    w = np.linalg.norm(s)
    theta = theta*np.pi/180

    if x == 1:
        g = 1
        a, b = np.sin(theta)/w, 2*(np.sin(theta/2))**2/w**2
    elif x == 2:
        g = 1/w
        a, b = np.sin(theta), 2*(np.sin(theta/2))**2
    elif x == 3:
        g = np.tan(theta/2)/w
        a, b = 2*(np.cos(theta/2))**2, 2*(np.cos(theta/2))**2
    elif x == 4:
        g = np.sin(theta/2)/w
        a, b = 2*np.cos(0.5*theta), 2
    else:
        g = theta/w
        a, b = np.sin(theta)/theta, 2*(np.sin(theta/2))**2/theta**2
    return g,a,b

def rotator(s,theta,gamma=5):

    # norm of pseudovector
    w = np.linalg.norm(s)

    # get normalization parameters

#####################################
s = np.array([1,1,1])
theta = 180.5
x = 5

start_time = time.time()

g, a, b= getparams(x,theta,s)
W = skew(s)
A = g*W
I = np.identity(3)
R = I + a*A + b*np.matmul(A,A)

print(R)

n,a = get_axial(R)
print(a)

print("--- %s seconds ---" % (time.time() - start_time))

