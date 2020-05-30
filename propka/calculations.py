"""PROPKA calculations."""
import math
import propka.protonate
import propka.bonds
from propka.lib import warning, info


# TODO - this file should be broken into three separate files:
# * calculations.py - includes basic functions for calculating distances, etc.
# * hydrogen.py - includes bonding and protonation functions
# * energy.py - includes energy functions (dependent on distance functions)


# TODO - the next set of functions form a distinct "module" for distance calculation


# Maximum distance used to bound calculations of smallest distance
MAX_DISTANCE = 1e6


def squared_distance(atom1, atom2):
    """Calculate the squared distance between two atoms.

    Args:
        atom1:  first atom for distance calculation
        atom2:  second atom for distance calculation
    Returns:
        distance squared
    """
    dx = atom2.x - atom1.x
    dy = atom2.y - atom1.y
    dz = atom2.z - atom1.z
    res = dx*dx+dy*dy+dz*dz
    return res


def distance(atom1, atom2):
    """Calculate the distance between two atoms.

    Args:
        atom1:  first atom for distance calculation
        atom2:  second atom for distance calculation
    Returns:
        distance
    """
    return math.sqrt(squared_distance(atom1, atom2))


def get_smallest_distance(atoms1, atoms2):
    """Calculate the smallest distance between two groups of atoms.

    Args:
        atoms1:  atom group 1
        atoms2:  atom group 2
    Returns:
        smallest distance between groups
    """
    res_dist = MAX_DISTANCE
    res_atom1 = None
    res_atom2 = None
    for atom1 in atoms1:
        for atom2 in atoms2:
            dist = squared_distance(atom1, atom2)
            if dist < res_dist:
                res_dist = dist
                res_atom1 = atom1
                res_atom2 = atom2
    return [res_atom1, math.sqrt(res_dist), res_atom2]
