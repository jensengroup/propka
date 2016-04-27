from __future__ import division
from __future__ import print_function
import math
from propka.lib import info

class vector:
    """ Vector """
    def __init__(self, xi=0.0, yi=0.0, zi=0.0, atom1 = 0, atom2 = 0):
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

        return

    def __add__(self, other):
        return vector(self.x + other.x,
                      self.y + other.y,
                      self.z + other.z)    

    def __sub__(self, other):
        return vector(self.x - other.x,
                      self.y - other.y,
                      self.z - other.z)
    
    def __mul__(self, other):
        """ Dot product, scalar and matrix multiplication """

        if isinstance(other,vector):
            return self.x * other.x + self.y * other.y + self.z * other.z
        elif isinstance(other, matrix4x4):
            return vector(
                xi = other.a11*self.x + other.a12*self.y + other.a13*self.z + other.a14*1.0,
                yi = other.a21*self.x + other.a22*self.y + other.a23*self.z + other.a24*1.0,
                zi = other.a31*self.x + other.a32*self.y + other.a33*self.z + other.a34*1.0
                )
        elif type(other) in [int, float]:
            return vector(self.x * other, self.y * other, self.z * other)
        else:
            info('%s not supported' % type(other))
            raise TypeError

    def __rmul__(self,other):
       return self.__mul__(other)
        

    def __pow__(self, other):
        """ Cross product """
        return vector(self.y * other.z - self.z * other.y,
                      self.z * other.x - self.x * other.z,
                      self.x * other.y - self.y * other.x)
        

    def __neg__(self):
        res = vector(xi = -self.x,
                     yi = -self.y,
                     zi = -self.z)
        return res

    def sq_length(self):
        return self.x * self.x + self.y * self.y + self.z * self.z

    def length(self):
        return math.sqrt(self.sq_length())

    def __str__(self):
        return '%10.4f %10.4f %10.4f'%(self.x, self.y, self.z)

    def __repr__(self):
        return '<vector>'

    def orthogonal(self):
        """ Returns a vector orthogonal to self """
        res = vector(self.y, -self.x, 0)
        
        if abs(self.y) < abs(self.z):
            res = vector(self.z, 0, -self.x)

        return res

    def rescale(self, new_length):
        """ Rescale vector to new length while preserving direction """
        frac = new_length/(self.length())
        res = vector(xi = self.x*frac,
                     yi = self.y*frac,
                     zi = self.z*frac)
        return res


class matrix4x4:
    def __init__(self,
                 a11i=0.0, a12i=0.0, a13i=0.0, a14i=0.0,
                 a21i=0.0, a22i=0.0, a23i=0.0, a24i=0.0,
                 a31i=0.0, a32i=0.0, a33i=0.0, a34i=0.0,
                 a41i=0.0, a42i=0.0, a43i=0.0, a44i=0.0):

        self.a11 = a11i
        self.a12 = a12i
        self.a13 = a13i
        self.a14 = a14i

        self.a21 = a21i
        self.a22 = a22i
        self.a23 = a23i
        self.a24 = a24i

        self.a31 = a31i
        self.a32 = a32i
        self.a33 = a33i
        self.a34 = a34i

        self.a41 = a41i
        self.a42 = a42i
        self.a43 = a43i
        self.a44 = a44i

        return





# methods working on vectors


def angle(a, b):
    dot = a * b
    return math.acos(dot / (a.length() * b.length()))


def angle_degrees(a,b):
    return math.degrees(angle(a, b))


def signed_angle_around_axis(a,b, axis):
    na = a**axis
    nb = b**axis

    v = angle(na,nb)

    d = b*(a**axis)

    if d < 0:
        v =-v

    return v

def signed_angle_degrees(a,b):
    return 180/math.pi * signed_angle(a, b)

        
def rotate_vector_around_an_axis(theta, axis, v):
        #print "# 1. rotate space about the z-axis so that the rotation axis lies in the xz-plane"
        gamma = 0.0
        if axis.y != 0:
            if axis.x != 0:
                gamma = -axis.x/abs(axis.x)*math.asin(axis.y/(math.sqrt(axis.x*axis.x + axis.y*axis.y)))
            else:
                gamma = math.pi/2.0

            Rz = rotate_atoms_around_z_axis(gamma)
            v = Rz * v
            axis = Rz * axis
            
        #print "# 2. rotate space about the y-axis so that the rotation axis lies along the z-axis"
        beta = 0.0
        if axis.x != 0:
            beta = -axis.x/abs(axis.x)*math.acos(axis.z/math.sqrt(axis.x*axis.x + axis.z*axis.z))
            Ry = rotate_atoms_around_y_axis(beta)
            v = Ry * v
            axis = Ry *axis

        #print "# 3. perform the desired rotation by theta about the z-axis"
        Rz = rotate_atoms_around_z_axis(theta)
        v = Rz * v

        #print "# 4. apply the inverse of step 2."
        Ry = rotate_atoms_around_y_axis(-beta)
        v = Ry * v
        
        #print "# 5. apply the inverse of step 1."
        Rz = rotate_atoms_around_z_axis(-gamma)
        v = Rz * v
        
        return v

def rotate_atoms_around_z_axis(angle):
    Rz = matrix4x4(
        a11i = math.cos(angle), a12i = -math.sin(angle), a13i = 0.0, a14i = 0.0,
        a21i = math.sin(angle), a22i =  math.cos(angle), a23i = 0.0, a24i = 0.0,
        a31i = 0.0            , a32i = 0.0             , a33i = 1.0, a34i = 0.0,
        a41i = 0.0            , a42i = 0.0             , a43i = 0.0, a44i = 1.0
        )
    
    return Rz


def rotate_atoms_around_y_axis(angle):
    Ry = matrix4x4(
        a11i =  math.cos(angle), a12i = 0.0, a13i = math.sin(angle), a14i = 0.0,
        a21i =  0.0            , a22i = 1.0, a23i = 0.0            , a24i = 0.0,
        a31i = -math.sin(angle), a32i = 0.0, a33i = math.cos(angle), a34i = 0.0,
        a41i =  0.0            , a42i = 0.0, a43i = 0.0            , a44i = 1.0
        )
        
    return Ry



class multi_vector:
    def __init__(self, atom1=0, atom2=0):
        self.vectors = []
        self.keys = []

        # store vectors for all configurations of atoms
        if atom1!=0:
            self.keys = lib.get_sorted_configurations(atom1.configurations.keys())
            if atom2!=0:
                keys2 = lib.get_sorted_configurations(atom2.configurations.keys())
                if self.keys != keys2:
                    raise 'Cannot make multi vector: Atomic configurations mismatch for\n   %s\n   %s\n'%(atom1,atom2)
            for key in self.keys:
                atom1.setConfiguration(key)
                if atom2!=0:
                    atom2.setConfiguration(key)
                v = vector(atom1=atom1, atom2=atom2)
                self.vectors.append(v)
                #info(key,v)
        return
    
    def __getattribute__(self,name):
        try:
            return object.__getattribute__(self, name)
        except AttributeError:
            return self.do_job(name)

    def __str__(self):
        res = ''
        for i in range(len(self.keys)):
            res += '%s %s\n'%(self.keys[i], self.vectors[i])
        return res


    def do_job(self, job):
        #info(job)
        self.res = multi_vector()
        for i in range(len(self.vectors)):
            self.res.vectors.append(eval('self.vectors[%d].%s()'%(i,job)))
            self.res.keys.append(self.keys[i])
        return self.get_result

    def get_result(self):
        return self.res
    
    def generic_operation(self, operation, other):
        if self.keys != other.keys:
            raise 'Incompatable keys'

        self.res = multi_vector()
        for i in range(len(self.vectors)):
            self.res.vectors.append(eval('self.vectors[%d] %s other.vectors[%d]'%(i,operation,i)))
            self.res.keys.append(self.keys[i])
        return

    def __add__(self, other):
        self.generic_operation('+',other)
        return self.res

    def __sub__(self, other):
        self.generic_operation('-',other)
        return self.res

    def __mul__(self, other):
        self.generic_operation('*',other)
        return self.res

    def __pow__(self, other):
        self.generic_operation('**',other)
        return self.res

    def generic_self_operation(self, operation):
        return

    def __neg__(self):
        self.generic_operation('*',-1.0)
        return self.res

    def rescale(self, new_length):
        self.res = multi_vector()
        for i in range(len(self.vectors)):
            self.res.vectors.append(self.vectors[i].rescale(new_length))
            self.res.keys.append(self.keys[i])
        return self.res


def rotate_multi_vector_around_an_axis(theta, axis, v):
    """ both axis ans v must be multi_vectors """
    
    if axis.keys != v.keys:
        raise 'Incompatible keys in rotate multi_vector'
    
    res = multi_vector()
    for i in range(len(v.keys)):
        res.vectors.append(rotate_vector_around_an_axis(theta, axis.vectors[i], v.vectors[i]))
        res.keys.append(v.keys[i])
        
    return res
