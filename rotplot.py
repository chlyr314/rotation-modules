import numpy as np
from mpl_toolkits import mplot3d
from numpy import *
from matplotlib import pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from matplotlib.patches import FancyArrowPatch
from mpl_toolkits.mplot3d import proj3d


# function that plots arrowheads on tips of vector lines
class Arrow3D(FancyArrowPatch):
    def __init__(self, xs, ys, zs, *args, **kwargs):
        FancyArrowPatch.__init__(self, (0,0), (0,0), *args, **kwargs)
        self._verts3d = xs, ys, zs

    def draw(self, renderer):
        xs3d, ys3d, zs3d = self._verts3d
        xs, ys, zs = proj3d.proj_transform(xs3d, ys3d, zs3d, renderer.M)
        self.set_positions((xs[0],ys[0]),(xs[1],ys[1]))
        FancyArrowPatch.draw(self, renderer)

# function to plot vector in 3D
def pvec(p,ax,c='k',lins='-'):
    """ 
    Plots vector p with arrow head. 
    INPUTS : 1) p : the end points of the vector
                 2) ax : figure handle
                 3) c : vector color (default is black)
                 4) lins : linestyle (default is continuous)
    """
    ax.plot3D([0,p[0]],[0,p[1]],[0,p[2]],linestyle = lins, color=c,
            linewidth=1)
    a = Arrow3D([0,p[0]],[0,p[1]],[0,p[2]],mutation_scale=20,lw=0,
                arrowstyle="-|>", color=c)
    ax.add_artist(a)

# function to connect points p1 and p2
def conp(p1,p2,ax,c='k',lins='-'):
    """ 
    Connects points p1 and p2 with lines in 3D
    INPUTS : 1) p1 : first point
             2) p2 : second point
             3) c : line color (default is black)
             4) lins : linestyle (default is continuous)
    """
    ax.plot3D([p1[0],p2[0]],[p1[1],p2[1]],[p1[2],p2[2]],linestyle = lins,
            color=c)


# function for plotting rotation path
def plotrot(r,v0,n,S,ax,d1flag=False,d2flag=False):

    """
    Function that plots rotation path from point v0 to point v.
    It will discretize the path in n steps and will apply the same 
    rotation sequentially : v = R v0 = dR*dR*dR....*dR*v0 where dR 
    is the rotation tensor that corresponds to dtheta = theta/n. This
    is done just to get the rotation path that results from the effect of
    the total rotation tensor R.

    INPUTS: 1) r : rotation matrix that corresponds to dtheta = theta/n
            2) v0 : vector to be rotated
            3) n  : number of steps. Total rotation theta = n*dtheta
            4) S  : Skew-symmetric matrix associated with dR. Unit axis of
                    rotation is extracted from it
            5) ax : current figure handle
            6) d1flag  : flag to plot 1st order approximation path
            7) d2flag  : flag to plot 2nd order approximation path
    
    OUTPUT: v : position of v0 after the effect of total rotation

    """
    # plot initial vector and axis of rotation vector
    s = np.array([-S[1,2],S[0,2],-S[0,1]]) # axis of rot
    pvec(s/np.linalg.norm(s),ax,'r','--')
    pvec(v0,ax,'k','--')

    vin = v0

    # plot rotation path
    for i in range(n):
        ax.scatter3D(vin[0],vin[1],vin[2],c='k',alpha=0.5,s=10)
        v = np.matmul(r,vin)
        conp(vin,v,ax)
        vin = v

    # plot rotated vector
    pvec(v,ax,'k')

    # output vector
    vout = v

    # plot path of first order approx to rotation R ~ I + S
    if d1flag:
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
    ax.set_xlim([-ub,ub])
    ax.set_ylim([-ub,ub])
    ax.set_zlim([-ub,ub])

    return vout
