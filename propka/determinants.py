
from __future__ import division
from __future__ import print_function

import math, time

import propka.iterative, propka.lib, propka.vector_algebra
import propka.calculations
from propka.determinant import Determinant


def setDeterminants(propka_groups, version=None, options=None):
    """
    adding side-chain and coulomb determinants/perturbations to all residues - note, backbone determinants are set separately
    """

    iterative_interactions = []
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
                interaction_type = version.parameters.interaction_matrix.get_value(group1.type,group2.type)
                if interaction_type == 'I':
                    propka.iterative.addtoDeterminantList(group1, group2, distance, iterative_interactions, version=version)
                elif interaction_type == 'N':
                    addDeterminants(group1, group2, distance, version)


    # --- Iterative section ---#
    propka.iterative.addDeterminants(iterative_interactions, version, options=options)


def addDeterminants(group1, group2, distance, version):
    """
    adding determinants/perturbations, distance(R1, R2) < coulomb_cutoff always
    """

    # side-chain determinant
    addSidechainDeterminants(group1, group2, version)

    # Coulomb determinant
    addCoulombDeterminants(group1, group2, distance, version)

    return

def addSidechainDeterminants(group1, group2, version=None):
    """
    adding side-chain determinants/perturbations
    Note, resNumb1 > resNumb2
    """

    hbond_interaction = version.hydrogen_bond_interaction(group1, group2)

    if hbond_interaction:

        if group1.charge == group2.charge:
            # acid pair or base pair
            if group1.model_pka < group2.model_pka:
                newDeterminant1 = Determinant(group2, -hbond_interaction)
                newDeterminant2 = Determinant(group1,  hbond_interaction)
            else:
                newDeterminant1 = Determinant(group2,  hbond_interaction)
                newDeterminant2 = Determinant(group1, -hbond_interaction)
        else:
            newDeterminant1 = Determinant(group2, hbond_interaction*group1.charge)
            newDeterminant2 = Determinant(group1, hbond_interaction*group2.charge)

        group1.determinants['sidechain'].append(newDeterminant1)
        group2.determinants['sidechain'].append(newDeterminant2)

    return

def addCoulombDeterminants(group1, group2, distance, version):
    """
    adding NonIterative Coulomb determinants/perturbations
    """

    coulomb_interaction = version.electrostatic_interaction(group1, group2, distance)

    if coulomb_interaction:
        Q1 = group1.charge
        Q2 = group2.charge

        # assigning the Coulombic interaction
        if   Q1 < 0.0 and Q2 < 0.0:
            """ both are acids """
            addCoulombAcidPair(group1, group2, coulomb_interaction)
        elif Q1 > 0.0 and Q2 > 0.0:
            """ both are bases """
            addCoulombBasePair(group1, group2, coulomb_interaction)
        else:
            """ one of each """
            addCoulombIonPair(group1, group2, coulomb_interaction)

    return


def addCoulombAcidPair(object1, object2, value):
    """
    Adding the Coulomb interaction (an acid pair):
    the higher pKa is raised
    """

    if object1.model_pka > object2.model_pka:
        newDeterminant = Determinant(object2, value)
        object1.determinants['coulomb'].append(newDeterminant)
    else:
        newDeterminant = Determinant(object1, value)
        object2.determinants['coulomb'].append(newDeterminant)


def addCoulombBasePair(object1, object2, value):
    """
    Adding the Coulomb interaction (a base pair):
    the lower pKa is lowered
    """
    if object1.model_pka < object2.model_pka:
        newDeterminant = Determinant(object2, -value)
        object1.determinants['coulomb'].append(newDeterminant)
    else:
        newDeterminant = Determinant(object1, -value)
        object2.determinants['coulomb'].append(newDeterminant)


def addCoulombIonPair(object1, object2, value):
    """
    Adding the Coulomb interaction (an acid-base pair):
    the pKa of the acid is lowered & the pKa of the base is raised
    """

    # residue1
    Q1 = object1.charge
    newDeterminant = Determinant(object2, Q1*value)
    object1.determinants['coulomb'].append(newDeterminant)

    # residue2
    Q2 = object2.charge
    newDeterminant = Determinant(object1, Q2*value)
    object2.determinants['coulomb'].append(newDeterminant)




def setIonDeterminants(conformation_container, version):
    """
    adding ion determinants/perturbations
    """
    for titratable_group in conformation_container.get_titratable_groups():
        for ion_group in conformation_container.get_ions():
            squared_distance = propka.calculations.squared_distance(titratable_group, ion_group)
            if squared_distance < version.parameters.coulomb_cutoff2_squared:
                weight = version.calculatePairWeight(titratable_group.Nmass, ion_group.Nmass)
                # the pKa of both acids and bases are shifted up by negative ions (and vice versa)
                value  =  (-ion_group.charge) * version.calculateCoulombEnergy(math.sqrt(squared_distance), weight)
                newDeterminant = Determinant(ion_group, value)
                titratable_group.determinants['coulomb'].append(newDeterminant)

    return

def setBackBoneDeterminants(titratable_groups, backbone_groups, version):

    for titratable_group in titratable_groups:
        titratable_group_interaction_atoms = titratable_group.interaction_atoms_for_acids
        if not titratable_group_interaction_atoms:
            continue

        # find out which backbone groups this titratable is interacting with
        for backbone_group in backbone_groups:
            # find the interacting atoms
            backbone_interaction_atoms = backbone_group.get_interaction_atoms(titratable_group)
            if not backbone_interaction_atoms:
                continue

            # find the smallest distance
            [backbone_atom, distance, titratable_atom] = propka.calculations.get_smallest_distance(backbone_interaction_atoms,
                                                                                                   titratable_group_interaction_atoms)
            # get the parameters
            parameters = version.get_backbone_hydrogen_bond_parameters(backbone_atom, titratable_atom)
            if not parameters:
                continue
            [dpKa_max, [cutoff1, cutoff2]] = parameters


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
                    if titratable_group.type in version.parameters.angular_dependent_sidechain_interactions:
                        if titratable_atom.element == 'H':
                            heavy_atom    = titratable_atom.bonded_atoms[0]
                            hydrogen_atom = titratable_atom
                            [d1, f_angle, d2] = propka.calculations.AngleFactorX(atom1=heavy_atom,
                                                                                 atom2=hydrogen_atom,
                                                                                 atom3=backbone_atom)
                        else:
                            # Either the structure is corrupt (no hydrogen), or the heavy atom is closer to
                            # the titratable atom than the hydrogen. In either case we set the angle factor
                            # to 0
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
                        backbone_N = backbone_atom.bonded_atoms[0]
                        backbone_H = backbone_atom
                        [d1, f_angle, d2] = propka.calculations.AngleFactorX(atom1=titratable_atom,
                                                                             atom2=backbone_H,
                                                                             atom3=backbone_N)
                    else:
                        # Either the structure is corrupt (no hydrogen), or the heavy atom is closer to
                        # the titratable atom than the hydrogen. In either case we set the angle factor
                        # to 0
                        f_angle = 0.0


                if f_angle > 0.001:
                    value = titratable_group.charge * propka.calculations.HydrogenBondEnergy(distance, dpKa_max, [cutoff1,cutoff2], f_angle)

                    newDeterminant = Determinant(backbone_group, value)
                    titratable_group.determinants['backbone'].append(newDeterminant)


    return
