"""Vector algebra for PROPKA."""
import math
from propka.lib import info, get_sorted_configurations


class Vector:
    """Vector"""

    def __init__(self, xi=0.0, yi=0.0, zi=0.0, atom1=None, atom2=None):
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

    def __add__(self, other):
        return Vector(self.x + other.x,
                      self.y + other.y,
                      self.z + other.z)

    def __sub__(self, other):
        return Vector(self.x - other.x,
                      self.y - other.y,
                      self.z - other.z)

    def __mul__(self, other):
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
        else:
            info('%s not supported' % type(other))
            raise TypeError

    def __rmul__(self, other):
        return self.__mul__(other)

    def __pow__(self, other):
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

    def length(self):
        """Return vector length."""
        return math.sqrt(self.sq_length())

    def __str__(self):
        return '%10.4f %10.4f %10.4f'%(self.x, self.y, self.z)

    def __repr__(self):
        return '<vector>'

    def orthogonal(self):
        """ Returns a vector orthogonal to self """
        res = Vector(self.y, -self.x, 0)
        if abs(self.y) < abs(self.z):
            res = Vector(self.z, 0, -self.x)
        return res

    def rescale(self, new_length):
        """ Rescale vector to new length while preserving direction """
        frac = new_length/(self.length())
        res = Vector(xi=self.x*frac,
                     yi=self.y*frac,
                     zi=self.z*frac)
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


def angle(avec, bvec):
    """Get the angle between two vectors.

    Args:
        avec:  vector 1
        bvec:  vector 2
    Returns:
        angle in radians
    """
    dot = avec * bvec
    return math.acos(dot / (avec.length() * bvec.length()))


def angle_degrees(avec, bvec):
    """Get the angle between two vectors in degrees.

    Args:
        avec:  vector 1
        bvec:  vector 2
    Returns:
        angle in degrees
    """
    return math.degrees(angle(avec, bvec))


def signed_angle_around_axis(avec, bvec, axis):
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


def rotate_vector_around_an_axis(theta, axis, vec):
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


def rotate_atoms_around_z_axis(theta):
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


def rotate_atoms_around_y_axis(theta):
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


class MultiVector:
    """Collection of vectors for multiple configurations of atoms.

    TODO - this class does not appear to be used or covered by tests
    """

    def __init__(self, atom1=None, atom2=None):
        """Initialize with atom configurations.

        Args:
            atom1:  first atom to define vector
            atom2:  second atom to define vector
        """
        self.vectors = []
        self.keys = []
        self.result = None
        # store vectors for all configurations of atoms
        if atom1 is not None:
            self.keys = get_sorted_configurations(atom1.configurations.keys())
            if atom2 is not None:
                keys2 = get_sorted_configurations(atom2.configurations.keys())
                if self.keys != keys2:
                    str_ = ('Cannot make multi vector: Atomic configurations '
                            'mismatch for\n   %s\n   %s\n' % (atom1, atom2))
                    raise KeyError(str_)
            for key in self.keys:
                atom1.setConfiguration(key)
                if atom2 != 0:
                    atom2.setConfiguration(key)
                vec = Vector(atom1=atom1, atom2=atom2)
                self.vectors.append(vec)

    def __getattribute__(self, name):
        try:
            return object.__getattribute__(self, name)
        except AttributeError:
            return self.do_job(name)

    def __str__(self):
        res = ''
        for i, key in enumerate(self.keys):
            res += '%s %s\n' % (key, self.vectors[i])
        return res

    def do_job(self, job):
        """Append vectors to configuration.

        Args:
            job:  name of function to apply to vectors
        Returns:
            TODO - figure out what this is
        """
        self.result = MultiVector()
        for i, vector in enumerate(self.vectors):
            func = getattr(vector, job)
            self.result.vectors.append(func())
            self.result.keys.append(self.keys[i])
        return self.get_result

    @property
    def get_result(self):
        """Return the latest result."""
        return self.result

    def generic_operation(self, operation, other):
        """Perform a generic operation between two MultiVector objects.

        Args:
            operation:  operation to perform (string)
            other:  other MultiVector object
        """
        if self.keys != other.keys:
            raise 'Incompatible keys'
        self.result = MultiVector()
        for i in range(len(self.vectors)):
            self.result.vectors.append(
                # TODO - eliminate eval() or entire class
                eval('self.vectors[%d] %s other.vectors[%d]'
                     % (i, operation, i)))
            self.result.keys.append(self.keys[i])

    def __add__(self, other):
        self.generic_operation('+', other)
        return self.result

    def __sub__(self, other):
        self.generic_operation('-', other)
        return self.result

    def __mul__(self, other):
        self.generic_operation('*', other)
        return self.result

    def __pow__(self, other):
        self.generic_operation('**', other)
        return self.result

    @staticmethod
    def generic_self_operation(_):
        """TODO - delete this."""
        return

    def __neg__(self):
        self.generic_operation('*', -1.0)
        return self.result

    def rescale(self, new_length):
        """Rescale multi-vector to new length.

        Args:
            new_length:  new length for multi-vector
        Result:
            MultiVector object
        """
        self.result = MultiVector()
        for i, vector in enumerate(self.vectors):
            self.result.vectors.append(vector.rescale(new_length))
            self.result.keys.append(self.keys[i])
        return self.res


def rotate_multi_vector_around_an_axis(theta, axis, vec):
    """Rotate a multi-vector around an axis.

    NOTE - both axis ans v must be MultiVectors.

    Args:
        theta:  angle (in radians)
        axis:  multi-vector axis
        vec:  multi-vector vector
    """
    if axis.keys != vec.keys:
        raise 'Incompatible keys in rotate MultiVector'
    res = MultiVector()
    for i, key in enumerate(vec.keys):
        res.vectors.append(rotate_vector_around_an_axis(
            theta, axis.vectors[i], vec.vectors[i]))
        res.keys.append(key)
    return res
