import propka.vector_algebra as m

import math
import pytest
from pytest import approx

RADIANS_90 = math.pi / 2


def assert_vector_equal(v1: m.Vector, v2: m.Vector):
    assert isinstance(v1, m.Vector)
    assert isinstance(v2, m.Vector)
    assert v1.x == approx(v2.x)
    assert v1.y == approx(v2.y)
    assert v1.z == approx(v2.z)


def _matrix4x4_tolist(self: m.Matrix4x4) -> list:
    return [
        self.a11, self.a12, self.a13, self.a14, self.a21, self.a22, self.a23, self.a24,
        self.a31, self.a32, self.a33, self.a34, self.a41, self.a42, self.a43, self.a44
    ]

def assert_matrix4x4_equal(m1: m.Matrix4x4, m2: m.Matrix4x4):
    assert isinstance(m1, m.Matrix4x4)
    assert isinstance(m2, m.Matrix4x4)
    assert _matrix4x4_tolist(m1) == approx(_matrix4x4_tolist(m2))


def test_Vector__init():
    v = m.Vector()
    assert v.x == 0.0
    assert v.y == 0.0
    assert v.z == 0.0
    v = m.Vector(12, 34, 56)
    assert v.x == 12
    assert v.y == 34
    assert v.z == 56
    v1 = m.Vector(atom1=v)
    assert v1.x == 12
    assert v1.y == 34
    assert v1.z == 56
    v2 = m.Vector(5, 4, 3)
    v3 = m.Vector(atom1=v2, atom2=v1)
    assert v3.x == 7
    assert v3.y == 30
    assert v3.z == 53


def test_Vector__plusminus():
    v1 = m.Vector(1, 2, 3)
    v2 = m.Vector(4, 5, 6)
    v3 = v1 + v2
    assert v3.x == 5
    assert v3.y == 7
    assert v3.z == 9
    v3 = v1 - v2
    assert v3.x == -3
    assert v3.y == -3
    assert v3.z == -3
    v4 = -v1
    assert v4.x == -1
    assert v4.y == -2
    assert v4.z == -3


def test_Vector__mul__number():
    v1 = m.Vector(1, 2, 3)
    assert_vector_equal(v1 * 2, m.Vector(2, 4, 6))


def test_Vector__mul__Vector():
    v1 = m.Vector(1, 2, 3)
    v2 = m.Vector(4, 5, 6)
    with pytest.deprecated_call():
        assert v1 * v2 == 32
    assert v1.dot(v2) == 32
    with pytest.raises(TypeError):
        v1 @ v2  # type: ignore


def test_Vector__mul__Matrix4x4():
    v1 = m.Vector(1, 2, 3)
    assert_vector_equal(m.Matrix4x4() @ v1, m.Vector())
    m2 = m.Matrix4x4(0, 1, 0, 0, 20, 0, 0, 0, 0, 0, 300, 0, 0, 0, 0, 1)
    with pytest.deprecated_call():
        assert_vector_equal(v1 * m2, m.Vector(2, 20, 900))
    with pytest.deprecated_call():
        assert_vector_equal(m2 * v1, m.Vector(2, 20, 900))
    assert_vector_equal(m2 @ v1, m.Vector(2, 20, 900))
    with pytest.raises(TypeError):
        v1 @ m2  # type: ignore


def test_Vector__cross():
    v1 = m.Vector(1, 2, 3)
    v2 = m.Vector(4, 5, 6)
    with pytest.deprecated_call():
        assert_vector_equal(v1**v2, m.Vector(-3, 6, -3))
    assert_vector_equal(v1.cross(v2), m.Vector(-3, 6, -3))
    assert_vector_equal(v2.cross(v1), m.Vector(3, -6, 3))


def test_Vector__length():
    v1 = m.Vector(1, 2, 3)
    assert v1.length() == 14**0.5
    assert v1.sq_length() == 14


def test_Vector__orthogonal():
    v1 = m.Vector(1, 2, 3)
    assert v1.dot(v1.orthogonal()) == 0


def test_Vector__rescale():
    v1 = m.Vector(1, 2, 3)
    v2 = v1.rescale(4)
    assert v2.length() == 4
    assert v2.x / v1.x == approx(4 / 14**0.5)
    assert v2.y / v1.y == approx(4 / 14**0.5)
    assert v2.z / v1.z == approx(4 / 14**0.5)


def test_angle():
    v1 = m.Vector(0, 0, 1)
    v2 = m.Vector(0, 2, 0)
    assert m.angle(v1, v2) == RADIANS_90


def test_angle_degrees():
    v1 = m.Vector(0, 0, 3)
    v2 = m.Vector(5, 0, 0)
    assert m.angle_degrees(v1, v2) == 90


def test_signed_angle_around_axis():
    v1 = m.Vector(0, 0, 3)
    v2 = m.Vector(5, 0, 0)
    v3 = m.Vector(0, 1, 0)
    assert m.signed_angle_around_axis(v1, v2, v3) == -RADIANS_90
    v1 = m.Vector(0, 2, 3)
    v2 = m.Vector(5, 4, 0)
    assert m.signed_angle_around_axis(v1, v2, v3) == -RADIANS_90


def test_rotate_vector_around_an_axis():
    v1 = m.Vector(0, 0, 3)
    v2 = m.Vector(3, 0, 0)
    v3 = m.Vector(0, -1, 0)
    v4 = m.rotate_vector_around_an_axis(RADIANS_90, v3, v2)
    assert_vector_equal(v4, v1)


def test_rotate_atoms_around_z_axis():
    m_rot = m.rotate_atoms_around_z_axis(-RADIANS_90)
    assert_matrix4x4_equal(m_rot,
                           m.Matrix4x4(0, 1, 0, 0, -1, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 1))


def test_rotate_atoms_around_y_axis():
    m_rot = m.rotate_atoms_around_y_axis(RADIANS_90)
    assert_matrix4x4_equal(m_rot,
                           m.Matrix4x4(0, 0, 1, 0, 0, 1, 0, 0, -1, 0, 0, 0, 0, 0, 0, 1))
