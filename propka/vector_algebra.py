"""
Vector calculations
===================

Vector algebra for PROPKA.
"""
import logging
import math
from typing import Optional, Protocol, Union


_LOGGER = logging.getLogger(__name__)


class _XYZ(Protocol):
    """
    Protocol for types which have x/y/z attributes, like Vector or Atom.
    """
    x: float
    y: float
    z: float


class Vector:
    """Vector"""

    x: float
    y: float
    z: float

    def __init__(self,
                 xi: float = 0.0,
                 yi: float = 0.0,
                 zi: float = 0.0,
                 atom1: Optional[_XYZ] = None,
                 atom2: Optional[_XYZ] = None):
        """Initialize vector.

        Args:
            xi:  default x-coordinate
            yi:  default y-coordinate
            zi:  default z-coordinate
            atom1:  atom center used to define default coordinate
            atom2:  two atom centers used to define vector
        """
        self.x = xi
        self.y = yi
        self.z = zi

        if atom1:
            # make vector pointing to atom1
            self.x = atom1.x
            self.y = atom1.y
            self.z = atom1.z

            if atom2:
                # make inter-atomic vector (atom1 -> atom2)
                self.x = atom2.x - self.x
                self.y = atom2.y - self.y
                self.z = atom2.z - self.z

    def __add__(self, other: _XYZ):
        return Vector(self.x + other.x,
                      self.y + other.y,
                      self.z + other.z)

    def __sub__(self, other: _XYZ):
        return Vector(self.x - other.x,
                      self.y - other.y,
                      self.z - other.z)

    def __mul__(self, other: Union["Vector", "Matrix4x4", float]):
        """Dot product, scalar and matrix multiplication."""
        if isinstance(other, Vector):
            return self.x * other.x + self.y * other.y + self.z * other.z
        elif isinstance(other, Matrix4x4):
            return Vector(
                xi=other.a11*self.x + other.a12*self.y + other.a13*self.z
                + other.a14*1.0,
                yi=other.a21*self.x + other.a22*self.y + other.a23*self.z
                + other.a24*1.0,
                zi=other.a31*self.x + other.a32*self.y + other.a33*self.z
                + other.a34*1.0
                )
        elif type(other) in [int, float]:
            return Vector(self.x * other, self.y * other, self.z * other)
        raise TypeError(f'{type(other)} not supported')

    def __rmul__(self, other):
        return self.__mul__(other)

    def __pow__(self, other: _XYZ):
        """Cross product."""
        return Vector(self.y * other.z - self.z * other.y,
                      self.z * other.x - self.x * other.z,
                      self.x * other.y - self.y * other.x)

    def __neg__(self):
        res = Vector(xi=-self.x,
                     yi=-self.y,
                     zi=-self.z)
        return res

    def sq_length(self):
        """Return vector squared-length"""
        return self.x * self.x + self.y * self.y + self.z * self.z

    def length(self) -> float:
        """Return vector length."""
        return math.sqrt(self.sq_length())

    def __str__(self):
        return '{0:>10.4f} {1:>10.4f} {2:>10.4f}'.format(
            self.x, self.y, self.z)

    def __repr__(self):
        return '<vector>'

    def orthogonal(self):
        """ Returns a vector orthogonal to self """
        res = Vector(self.y, -self.x, 0)
        if abs(self.y) < abs(self.z):
            res = Vector(self.z, 0, -self.x)
        return res

    def rescale(self, new_length: float):
        """ Rescale vector to new length while preserving direction """
        frac = new_length/(self.length())
        res = Vector(xi=self.x*frac, yi=self.y*frac, zi=self.z*frac)
        return res


class Matrix4x4:
    """A 4-by-4 matrix class."""

    def __init__(self,
                 a11i=0.0, a12i=0.0, a13i=0.0, a14i=0.0,
                 a21i=0.0, a22i=0.0, a23i=0.0, a24i=0.0,
                 a31i=0.0, a32i=0.0, a33i=0.0, a34i=0.0,
                 a41i=0.0, a42i=0.0, a43i=0.0, a44i=0.0):
        """Initialize with matrix elements."""
        # Row 1
        self.a11 = a11i
        self.a12 = a12i
        self.a13 = a13i
        self.a14 = a14i
        # Row 2
        self.a21 = a21i
        self.a22 = a22i
        self.a23 = a23i
        self.a24 = a24i
        # Row 3
        self.a31 = a31i
        self.a32 = a32i
        self.a33 = a33i
        self.a34 = a34i
        # Row 4
        self.a41 = a41i
        self.a42 = a42i
        self.a43 = a43i
        self.a44 = a44i


def angle(avec: Vector, bvec: Vector) -> float:
    """Get the angle between two vectors.

    Args:
        avec:  vector 1
        bvec:  vector 2
    Returns:
        angle in radians
    """
    dot = avec * bvec
    return math.acos(dot / (avec.length() * bvec.length()))


def angle_degrees(avec: Vector, bvec: Vector) -> float:
    """Get the angle between two vectors in degrees.

    Args:
        avec:  vector 1
        bvec:  vector 2
    Returns:
        angle in degrees
    """
    return math.degrees(angle(avec, bvec))


def signed_angle_around_axis(avec: Vector, bvec: Vector, axis: Vector) -> float:
    """Get signed angle of two vectors around axis in radians.

    Args:
        avec:  vector 1
        bvec:  vector 2
        axis:  axis
    Returns:
        angle in radians
    """
    norma = avec**axis
    normb = bvec**axis
    ang = angle(norma, normb)
    dot_ = bvec*(avec**axis)
    if dot_ < 0:
        ang = -ang
    return ang


def rotate_vector_around_an_axis(theta: float, axis: Vector, vec: Vector) -> Vector:
    """Rotate vector around an axis.

    Args:
        theta:  rotation angle (in radians)
        axis:  axis for rotation
        vec:  vector to rotate
    Returns:
        rotated vector
    """
    gamma = 0.0
    if axis.y != 0:
        if axis.x != 0:
            gamma = -axis.x/abs(axis.x)*math.asin(
                axis.y/(math.sqrt(axis.x*axis.x + axis.y*axis.y)))
        else:
            gamma = math.pi/2.0
        rot_z = rotate_atoms_around_z_axis(gamma)
        vec = rot_z * vec
        axis = rot_z * axis
    beta = 0.0
    if axis.x != 0:
        beta = -axis.x/abs(axis.x)*math.acos(
            axis.z/math.sqrt(axis.x*axis.x + axis.z*axis.z))
        rot_y = rotate_atoms_around_y_axis(beta)
        vec = rot_y * vec
        axis = rot_y * axis
    rot_z = rotate_atoms_around_z_axis(theta)
    vec = rot_z * vec
    rot_y = rotate_atoms_around_y_axis(-beta)
    vec = rot_y * vec
    rot_z = rotate_atoms_around_z_axis(-gamma)
    vec = rot_z * vec
    return vec


def rotate_atoms_around_z_axis(theta: float) -> Matrix4x4:
    """Get rotation matrix for z-axis.

    Args:
        theta:  angle of rotation (radians)
    Returns:
        rotation matrix
    """
    return Matrix4x4(
        a11i=math.cos(theta),
        a12i=-math.sin(theta),
        a13i=0.0,
        a14i=0.0,
        a21i=math.sin(theta),
        a22i=math.cos(theta),
        a23i=0.0,
        a24i=0.0,
        a31i=0.0,
        a32i=0.0,
        a33i=1.0,
        a34i=0.0,
        a41i=0.0,
        a42i=0.0,
        a43i=0.0,
        a44i=1.0
        )


def rotate_atoms_around_y_axis(theta: float) -> Matrix4x4:
    """Get rotation matrix for y-axis.

    Args:
        theta:  angle of rotation (radians)
    Returns:
        rotation matrix
    """
    return Matrix4x4(
        a11i=math.cos(theta),
        a12i=0.0,
        a13i=math.sin(theta),
        a14i=0.0,
        a21i=0.0,
        a22i=1.0,
        a23i=0.0,
        a24i=0.0,
        a31i=-math.sin(theta),
        a32i=0.0,
        a33i=math.cos(theta),
        a34i=0.0,
        a41i=0.0,
        a42i=0.0,
        a43i=0.0,
        a44i=1.0
        )
