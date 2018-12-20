import numpy as np
import time
import sys

def get_rotmat(a, theta = None, angtype = 'degrees'):
    """ Returns rotation tensor, Spin(skew) tensor
        and tangent operator from axial rotation vector.
        The tangent operator maps the non-additive increment
        of the axial vector to an additive increment to the
        current axial vector

        Input: 1) <a>       is the axial vector (axis of rotation).

               2) <theta>   is the angle of rotation around <a>, assumed
                            positive for counter-clockwise rotation. If value
                            is None, then the total rotation  is the norm and
                            the input vector <a> needs to be decomposed into
                            a unit vector and a scalar multiple.
                            If theta is specified, then <a>
                            needs only be normalized. Default value is None.

               3) <angtype> is the unit for the input rotation. Default is
                            degrees. Options: a)'degrees'
                                              b)'rads'

        Output: 1) Rotation tensor <R>
                2) Skew matrix <S>
                3) Tangent Operator <T>
        ######################################################################
        Changelog: - Created
                   - Added Tangent Map Operator (Krenk 3.47)
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
        A = np.array([[0, -a[2], a[1]], [a[2], 0, -a[0]], [-a[1], a[0], 0]])

        # Outer product tensor
        O = np.outer(a, a)
        A2=O-np.identity(3)

        # Get rotation tensor
        R = np.cos(theta)*np.identity(3) + (1-np.cos(theta))*O + np.sin(theta)*A
        S = theta*A

        # Get tangent operator T. T maps an increment to the axial vector
        # (which is non-additive) to the incremental rotation vector, which
        # is additive
        cot = 1/np.tan(0.5*theta)
        T = O + 0.5*theta*cot*(np.identity(3)-O) - 0.5*theta*A

        return R, S, T


def updateAxialVector(s0, ds, T0):
    """ Returns the updated axial vector and the corresponding
        rotation tensor

        *****************************************************
        Input:  1) <s0>  -    current axial vector

                2) <ds>  -    variation of axial vector (non-additive)

                3) <T0>  -    Rotational tangent operator
        *****************************************************
        Output: 1) <s>   -    updated axial vector

                2) <R>   -    updated rotation tensor

                3) <T>   -    updated tangent operator

    """

    # Additive incremental rotation vector
    a_ds = np.matmul(T0, ds)

    # Total updated axial vector
    s = s0 + a_ds

    # Total (compound) rotation
    theta = np.linalg.norm(s)*180.0/np.pi

    # Compound axis of rotation
    e = s/np.linalg.norm(s)

    # Compound rotation tensor
    R, _, T = get_rotmat(e, theta)

    return s, R, T



def get_axial(R,s=None):
    """ Returns axial unit vector and angle of (counterclockwise) rotation

       Input:  1) <R>   -    is the rotation matrix

               2) <s>   -    flag to return axial vector if specified.
                             If not specified, it returns unit vector (axis)
                             and angle of rotation (in degrees)

       Output: 1) <e>   -    axial vector if <s> is specified

                                   OR

               1) <n>   -    axial unit vector with components the rotations
                             with respect to the three axes
               2) <a>   -    angle of rotation in degrees
       if <s> is not specified

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
    if np.round(np.linalg.det(R), 3) != 1.0:
        sys.exit('Determinant of rotation matrix is not 1. Exiting')

    # Step 1 - Calculate cos of angle from tr(R) relation
    cos = 0.5*(np.trace(R)-1)

    # Step 2 - Cases : a =0, a= 180, 0 < a < 180
    if np.round(cos, 8) == 1:                            # No Rotation, a = 0
        a = 0
        n=np.array([1, 0, 0])
    elif np.round(cos, 8) == -1:                   # Angle of rotation a = 180

        a=np.pi

        print('to pire ws 180!!')
        # Columns of R+I are parallel to axis of rotation!
        M = R + np.identity(3)

        # Pick the first non-zero column vector of M = R + I
        i = 0
        normn=np.linalg.norm(M[0:3, i])

        while normn==0:
            i = i+1
            normn = np.linalg.norm(M[0:3, i])

        # Normalize it so it becomes unit
        n = 1/normn*(M[0:3, i])

    else:

        # Get angle of rotation from arccos - unique for 0<a<180
        a = np.arccos((np.trace(R)-1)/2)

        # Calculate sina and components of axial unit vector from skew(R)
        sin = np.sqrt(1-cos**2)
        par = 1/(2*sin)

        a1 = par*(R[2, 1]-R[1, 2])
        a2 = par*(R[0, 2]-R[2, 0])
        a3 = par*(R[1, 0]-R[0, 1])

        n=np.array([a1, a2, a3])

    if s:
        return n*a # returns axial vector
    else:
        a = 180*a/np.pi # returns unit vector < e > and angle of rotation
        return n,a      # in degrees


def vec2qua(n, theta):
    """
    INPUT:
           <n> : vector representing the axis of rotation
           <theta> : scalar representing the rotation in degrees around <n>

    OUTPUT:
           <q> : unit quaternion representing the rotation
           <M> : 3x3 matrix associated with this rotation

    Description: Function that accepts a vector <n> and a rotation value,
    <theta> - in degrees! - and returns the corresponding unit quaternion,
    along with the rotation matrix associated with <n> and <theta>
    """
    # Normalize vector and convert degrees to rad
    n = n/np.linalg.norm(n)
    th = theta*np.pi/180

    # Quaternion components. Re(q) = w
    w, x, y, z = (np.cos(th/2), n[0]*np.sin(th/2), n[1]*np.sin(th/2),
                     n[2]*np.sin(th/2))

    # Quaternion form q = [w, x, y, z]
    q = np.array([w, x, y, z])

    # get rotation matrix from quaternion
    M = qua2mat(q)

    return q, M


def qua2mat(q):
    """
    INPUT: <q> : unit quaternion

    OUTPUT: <M> : Rotation matrix associated with unit quaternion

    Description: Function that returns the rotation matrix associated with
                 the unit quaternion <q>
    """
    w, x, y, z = q[0],q[1],q[2],q[3]

    # Matrix Components
    m11, m12, m13 = 1-2*y**2-2*z**2, 2*x*y-2*w*z, 2*x*z+2*w*y
    m21, m22, m23 = 2*x*y+2*w*z, 1-2*x**2-2*z**2, 2*z*y-2*w*x
    m31, m32, m33 = 2*x*z-2*w*y, 2*y*z+2*w*x, 1-2*x**2-2*y**2

    M = np.array([[m11, m12, m13], [m21, m22, m23], [m31, m32, m33]])

    return M


def slerp(q1,q2,t):
    """
    INPUT: <q1> : unit quaternion at start
           <q2> : unit quaternion at end
            <t> : parameter 0<t<1. For t=0, q=q1. For t=1, q=q2.

    OUTPUT: <q> : interpolated quaternion in between

    Description: Function that performs spherical linear interpolation of
                 quaternions q1 and q2.
    """

    # Angle between q1 and q2
    cosz = q1.dot(q2)
    z = np.arccos(cosz)

    # slerp coefficients
    a = np.sin((1-t)*z)/np.sin(z)
    b = np.sin(t*z)/np.sin(z)

    q = a*q1 + b*q2

    return q
