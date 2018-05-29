import numpy as np
import time
import sys

def get_rotmat(a,theta = None,angtype = 'degrees'):
    """ Returns rotation tensor from axial rotation vector

        Input: 1) <a>       is the axial vector (axis of rotation).

               2) <theta>   is the angle of rotation around <a>, assumed
                            positive for counter-clockwise rotation. If value
                            is None, then the total rotation  is the norm and 
                            the input vector <a> needs to be decomposed into 
                            a unit vector and a scalar multiple.
                            If theta is specified,then <a>
                            needs only be normalized. Default value is None.

               3) <angtype> is the unit for the input rotation. Default is
                            degrees. Options: a)'degrees'
                                              b)'rads'

        Output: 1) Rotation tensor <R>
                2) Skew matrix <S> 
    """

    # Get norm of axial vector
    norms = np.linalg.norm(a)

    # Decompose input into unit vector (axis) and average rotation around it
    # if rotation is given in the form of components of <a>

    # Check arguments and act accordingly
    if theta is None:
        theta = norms # Extract rotation around axial vector = norm
    else:
        try:
            dum = theta*np.pi/180.
            if angtype is 'degrees':
                theta = dum
            else:
                print('<Angtype> argument assumed "rads".')
        except TypeError:
            sys.exit('ERROR: wrong type of argument <theta>.')

    # If axis of rotation vector not unit, make it unit
    if norms != 1.0:
        a = a/norms

    # If rotation is zero, then rotation tensor is the identity tensor,
    # else, use exponential map to return rotation tensor [R]
    if theta == 0.0:
        print('No rotation. [R] = [I] \n')
        return np.identity(3)

    else:

        # Get skew symmetric axial tensor associated with a
        A = np.array([[0,-a[2],a[1]],[a[2],0,-a[0]],[-a[1],a[0],0]])

        # Outer product tensor
        O = np.outer(a,a)
        A2=O-np.identity(3)

        # Get rotation tensor
        R = np.cos(theta)*np.identity(3) + (1-np.cos(theta))*O + np.sin(theta)*A
        S = theta*A
        return R,S




def get_axial(R):
    """ Returns axial unit vector and angle of (counterclockwise) rotation

       Input:  1) <R>   -    is the rotation matrix

       Output: 1) <n>   -    axial unit vector with components the rotations
                             with respect to the three axes
               2) <a>   -    angle of rotation in degrees

       Description:      Perform check to see if determinant of input is 1.
                         Determine the cos<a> from the invariant relation.
                         Then, a) if cos<a>=1 --> <a>=0 and R=I
                               b) if cos<a>=-1 --> <a>=180 and <n> is
                                  determined from the nonzero column of
                                  matrix (R+I)
                               c) if 0< <a> <180, then determine a from
                                  skew-symmetric part of R using, making
                                  use of the Euler-Rodriguez formula.
                                  The direction of <n> is uniquely determined
                                  from the components of the axial tensor.

                         CHECK ISSUES:  Both the determinand and the cos<a>
                                        that are used in checks are rounded
                                        to the third decimal.
    """
    if np.round(np.linalg.det(R),3) != 1.0:
        sys.exit('Determinant of rotation matrix is not 1. Exiting')

    # Step 1 - Calculate cos of angle from tr(R) relation
    cos = 0.5*(np.trace(R)-1)

    # Step 2 - Cases : a =0, a= 180, 0 < a < 180
    if np.round(cos,8) == 1:                            # No Rotation, a = 0
        a = 0
        n=np.array([1,0,0])
    elif np.round(cos,8) == -1:                   # Angle of rotation a = 180

        a=np.pi

        print('to pire ws 180!!')
        # Columns of R+I are parallel to axis of rotation!
        M = R + np.identity(3)

        # Pick the first non-zero column vector of M = R + I
        i = 0
        normn=np.linalg.norm(M[0:3,i])

        while normn==0:
            i = i+1
            normn = np.linalg.norm(M[0:3,i])

        # Normalize it so it becomes unit
        n = 1/normn*(M[0:3,i])

    else:

        # Get angle of rotation from arccos - unique for 0<a<180
        a = np.arccos((np.trace(R)-1)/2)

        # Calculate sina and components of axial unit vector from skew(R)
        sin = np.sqrt(1-cos**2)
        par = 1/(2*sin)

        a1 = par*(R[2,1]-R[1,2])
        a2 = par*(R[0,2]-R[2,0])
        a3 = par*(R[1,0]-R[0,1])

        n=np.array([a1,a2,a3])

    a = 180*a/np.pi
    return n,a



