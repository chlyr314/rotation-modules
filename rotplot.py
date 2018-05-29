import numpy as np
from mpl_toolkits import mplot3d
from numpy import *
from matplotlib import pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from matplotlib.patches import FancyArrowPatch
from mpl_toolkits.mplot3d import proj3d

class Arrow3D(FancyArrowPatch):
    def __init__(self, xs, ys, zs, *args, **kwargs):
        FancyArrowPatch.__init__(self, (0,0), (0,0), *args, **kwargs)
        self._verts3d = xs, ys, zs

    def draw(self, renderer):
        xs3d, ys3d, zs3d = self._verts3d
        xs, ys, zs = proj3d.proj_transform(xs3d, ys3d, zs3d, renderer.M)
        self.set_positions((xs[0],ys[0]),(xs[1],ys[1]))
        FancyArrowPatch.draw(self, renderer)

# function to plot marker in 3D
def pvec(p,ax,c='k',lins='-'):

    ax.plot3D([0,p[0]],[0,p[1]],[0,p[2]],linestyle = lins, color=c,
            linewidth=1)
    a = Arrow3D([0,p[0]],[0,p[1]],[0,p[2]],mutation_scale=20,lw=0,
                arrowstyle="-|>", color=c)
    ax.add_artist(a)

# function to connect points p1 and p2
def conp(p1,p2,ax,c='k',lins='-'):

    ax.plot3D([p1[0],p2[0]],[p1[1],p2[1]],[p1[2],p2[2]],linestyle = lins,
            color=c)


# function for plotting rotation path
def plotrot(r,v0,n,S,ax,d2flag=False):


    # plot initial vector and axis of rotation vector
    s = np.array([-S[1,2],S[0,2],-S[0,1]]) # axis of rot
    pvec(s/np.linalg.norm(s),ax,'r','--')
    pvec(v0,ax,'k','--')
    vin = v0

    # plot rotation path
    for i in range(n):
        ax.scatter3D(vin[0],vin[1],vin[2],c='k')
        v = np.matmul(r,vin)
        conp(vin,v,ax)
        vin = v

    # plot rotated vector
    pvec(v,ax,'k')
    vout = v
    # plot path of first order approx to rotation R ~ I + S
    vin = v0
    for j in range(n):
        cv = np.matmul(S,vin)
        v1 = vin+cv
        conp(vin,v1,ax,'b','--')
        ax.scatter3D(v1[0],v1[1],v1[2],c='b')
        vin = v1

    # plot path of second order approx to rotation
    if d2flag:
        vin = v0
        for k in range(n+1):
            cv = np.matmul(S,vin)
            cv2 = 0.5*np.matmul(S,cv)
            v1 = vin+cv+cv2
            conp(vin,v1,ax,'r','--')
            ax.scatter3D(v1[0],v1[1],v1[2],c='r')
            vin = v1
    ub = np.linalg.norm(v0)
    ax.set_xlabel('x')
    ax.set_ylabel('y')
    ax.set_zlabel('z')
    ax.set_xlim([-1,ub])
    ax.set_ylim([-1,ub])
    ax.set_zlim([-1,ub])
    return vout
