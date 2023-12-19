"""
Working with Determinants
=========================

Functions to manipulate :class:`propka.determinant.Determinant` objects.

.. TODO::

   It is confusing to have both `determinant.py` and `determinants.py`.
   Should these be merged?

.. SeeAlso::
   :mod:`propka.determinant`

"""
import math
from typing import List

import propka.calculations
import propka.iterative
import propka.lib
import propka.vector_algebra
from propka.calculations import squared_distance, get_smallest_distance
from propka.energy import angle_distance_factors, hydrogen_bond_energy
from propka.determinant import Determinant
from propka.group import Group
from propka.iterative import Interaction
from propka.version import Version


# Cutoff for angle factor
# TODO - this constant appears elsewhere in the package.
# It should be moved to a configuration file.
FANGLE_MIN = 0.001


def set_determinants(propka_groups: List[Group], version: Version, options=None):
    """Add side-chain and coulomb determinants/perturbations to all residues.

    NOTE - backbone determinants are set separately

    Args:
        propka_groups:  groups to adjust
        version:  version object
        options:  options object
    """
    iterative_interactions: List[Interaction] = []
    # --- NonIterative section ---#
    for group1 in propka_groups:
        for group2 in propka_groups:
            if group1 == group2:
                break
            # do not calculate interactions for coupled groups
            if group2 in group1.covalently_coupled_groups:
                break
            distance = propka.calculations.distance(group1, group2)
            if distance < version.parameters.coulomb_cutoff2:
                interaction_type = (
                    version.parameters.interaction_matrix.get_value(
                        group1.type, group2.type))
                if interaction_type == 'I':
                    propka.iterative.add_to_determinant_list(
                        group1, group2, distance, iterative_interactions,
                        version=version)
                elif interaction_type == 'N':
                    add_determinants(group1, group2, distance, version)
    # --- Iterative section ---#
    propka.iterative.add_determinants(iterative_interactions, version)


def add_determinants(group1, group2, distance, version):
    """Add determinants and perturbations for distance(R1,R2) < coulomb_cutoff.

    Args:
        group1:  first group to add
        group2:  second group to add
        distance:  distance between groups
        version:  version object
    """
    # side-chain determinant
    add_sidechain_determinants(group1, group2, version)
    # Coulomb determinant
    add_coulomb_determinants(group1, group2, distance, version)


def add_sidechain_determinants(group1: Group, group2: Group, version: Version):
    """Add side-chain determinants and perturbations.

    NOTE - res_num1 > res_num2

    Args:
        group1:  first group to add
        group2:  second group to add
        version:  version object
    """
    hbond_interaction = version.hydrogen_bond_interaction(group1, group2)
    if hbond_interaction:
        if group1.charge == group2.charge:
            # acid pair or base pair
            if group1.model_pka < group2.model_pka:
                new_determinant1 = Determinant(group2, -hbond_interaction)
                new_determinant2 = Determinant(group1, hbond_interaction)
            else:
                new_determinant1 = Determinant(group2, hbond_interaction)
                new_determinant2 = Determinant(group1, -hbond_interaction)
        else:
            new_determinant1 = Determinant(
                group2, hbond_interaction*group1.charge)
            new_determinant2 = Determinant(
                group1, hbond_interaction*group2.charge)
        group1.determinants['sidechain'].append(new_determinant1)
        group2.determinants['sidechain'].append(new_determinant2)


def add_coulomb_determinants(group1, group2, distance, version):
    """Add non-iterative Coulomb determinants and perturbations.

    Args:
        group1:  first group to add
        group2:  second group to add
        distance:  distance between groups
        version:  version object
    """
    coulomb_interaction = version.electrostatic_interaction(
        group1, group2, distance)
    if coulomb_interaction:
        q1 = group1.charge
        q2 = group2.charge
        # assigning the Coulombic interaction
        if   q1 < 0.0 and q2 < 0.0:
            # both are acids
            add_coulomb_acid_pair(group1, group2, coulomb_interaction)
        elif q1 > 0.0 and q2 > 0.0:
            # both are bases
            add_coulomb_base_pair(group1, group2, coulomb_interaction)
        else:
            # one of each
            add_coulomb_ion_pair(group1, group2, coulomb_interaction)


def add_coulomb_acid_pair(object1, object2, value):
    """Add the Coulomb interaction (an acid pair).

    The higher pKa is raised.

    Args:
        object1:  first part of pair
        object2:  second part of pair
        value:  determinant value
    """
    if object1.model_pka > object2.model_pka:
        new_determinant = Determinant(object2, value)
        object1.determinants['coulomb'].append(new_determinant)
    else:
        new_determinant = Determinant(object1, value)
        object2.determinants['coulomb'].append(new_determinant)


def add_coulomb_base_pair(object1, object2, value):
    """Add the Coulomb interaction (a base pair).

    The lower pKa is lowered.

    Args:
        object1:  first part of pair
        object2:  second part of pair
        value:  determinant value
    """
    if object1.model_pka < object2.model_pka:
        new_determinant = Determinant(object2, -value)
        object1.determinants['coulomb'].append(new_determinant)
    else:
        new_determinant = Determinant(object1, -value)
        object2.determinants['coulomb'].append(new_determinant)


def add_coulomb_ion_pair(object1, object2, value):
    """Add the Coulomb interaction (an acid-base pair).

    The pKa of the acid is lowered & the pKa of the base is raised.

    Args:
        object1:  first part of pair
        object2:  second part of pair
        value:  determinant value
    """
    # residue1
    q1 = object1.charge
    new_determinant = Determinant(object2, q1*value)
    object1.determinants['coulomb'].append(new_determinant)
    # residue2
    q2 = object2.charge
    new_determinant = Determinant(object1, q2*value)
    object2.determinants['coulomb'].append(new_determinant)


def set_ion_determinants(conformation_container, version):
    """Add ion determinants and perturbations.

    Args:
        conformation_container:  conformation to set
        version:  version object
    """
    for titratable_group in conformation_container.get_titratable_groups():
        for ion_group in conformation_container.get_ions():
            dist_sq = squared_distance(titratable_group, ion_group)
            if dist_sq < version.parameters.coulomb_cutoff2_squared:
                weight = version.calculate_pair_weight(
                    titratable_group.num_volume, ion_group.num_volume)
                # the pKa of both acids and bases are shifted up by negative
                # ions (and vice versa)
                value = (
                    -ion_group.charge
                    * version.calculate_coulomb_energy(
                        math.sqrt(dist_sq), weight))
                new_det = Determinant(ion_group, value)
                titratable_group.determinants['coulomb'].append(new_det)


def set_backbone_determinants(titratable_groups, backbone_groups, version):
    """Set determinants between titrable and backbone groups.

    Args:
        titratable_groups:  list of titratable groups
        backbone_groups:  list of backbone groups
        version:  version object
    """
    for titratable_group in titratable_groups:
        titratable_group_interaction_atoms = (
            titratable_group.interaction_atoms_for_acids)
        if not titratable_group_interaction_atoms:
            continue
        # find out which backbone groups this titratable is interacting with
        for backbone_group in backbone_groups:
            # find the interacting atoms
            backbone_interaction_atoms = (
                backbone_group.get_interaction_atoms(titratable_group))
            if not backbone_interaction_atoms:
                continue
            # find the smallest distance
            [backbone_atom, distance, titratable_atom] = (
                get_smallest_distance(
                    backbone_interaction_atoms,
                    titratable_group_interaction_atoms))
            assert backbone_atom is not None
            assert titratable_atom is not None
            # get the parameters
            parameters = (
                version.get_backbone_hydrogen_bond_parameters(
                    backbone_atom, titratable_atom))
            if not parameters:
                continue
            [dpka_max, [cutoff1, cutoff2]] = parameters
            if distance < cutoff2:
                # calculate angle factor
                f_angle = 1.0
                # for BBC groups, the hydrogen is on the titratable group
                #
                #        Titra.
                #       /
                #      H
                #     .
                #    O
                #    ||
                #    C
                if backbone_group.type == 'BBC':
                    if (titratable_group.type
                            in version.parameters.angular_dependent_sidechain_interactions):
                        if titratable_atom.element == 'H':
                            heavy_atom = titratable_atom.bonded_atoms[0]
                            hydrogen_atom = titratable_atom
                            [_, f_angle, _] = angle_distance_factors(
                                atom1=heavy_atom, atom2=hydrogen_atom,
                                atom3=backbone_atom)
                        else:
                            # Either the structure is corrupt (no hydrogen),
                            # or the heavy atom is closer to the titratable
                            # atom than the hydrogen. In either case we set the
                            # angle factor to 0
                            f_angle = 0.0
                # for BBN groups, the hydrogen is on the backbone group
                #
                #      Titra.
                #     .
                #    H
                #    |
                #    N
                #  /   \
                if backbone_group.type == 'BBN':
                    if backbone_atom.element == 'H':
                        backbone_n = backbone_atom.bonded_atoms[0]
                        backbone_h = backbone_atom
                        [_, f_angle, _] = (
                            angle_distance_factors(
                                atom1=titratable_atom, atom2=backbone_h,
                                atom3=backbone_n))
                    else:
                        # Either the structure is corrupt (no hydrogen), or the
                        # heavy atom is closer to the titratable atom than the
                        # hydrogen. In either case we set the angle factor to 0
                        f_angle = 0.0
                if f_angle > FANGLE_MIN:
                    value = (
                        titratable_group.charge
                        * hydrogen_bond_energy(
                            distance, dpka_max, [cutoff1, cutoff2], f_angle))
                    new_determinant = Determinant(backbone_group, value)
                    titratable_group.determinants['backbone'].append(
                        new_determinant)
