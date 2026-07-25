"""Tests for ligand geometry helpers."""

from propka.ligand import are_atoms_planar
from propka.vector_algebra import Vector


def test_are_atoms_planar_with_degenerate_normal():
    atoms = [
        Vector(0.0, 0.0, 0.0),
        Vector(1.0, 0.0, 0.0),
        Vector(1.0, 0.0, 0.0),
        Vector(0.0, 1.0, 0.0),
    ]

    assert are_atoms_planar(atoms) is False


def test_are_atoms_planar_with_duplicate_test_atom():
    atoms = [
        Vector(0.0, 0.0, 0.0),
        Vector(1.0, 0.0, 0.0),
        Vector(0.0, 1.0, 0.0),
        Vector(0.0, 0.0, 0.0),
    ]

    assert are_atoms_planar(atoms) is False
