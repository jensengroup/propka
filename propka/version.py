from __future__ import division
from __future__ import print_function
import math
import sys, os

import propka.lib as lib
from propka.lib import info, warning
import propka.calculations as calculations
import propka.parameters


class version:
    def __init__(self,parameters):
        self.parameters = parameters
        return

    # desolvation
    def calculate_desolvation(self, group):
        return self.desolvation_model(self.parameters, group)

    def calculatePairWeight(self, Nmass1, Nmass2):
        return self.weight_pair_method(self.parameters, Nmass1, Nmass2)

    # side chains
    def hydrogen_bond_interaction(self, group1, group2):
        return self.hydrogen_bond_interaction_model(group1, group2, self)

    def calculateSideChainEnergy(self, distance, dpka_max, cutoff, weight, f_angle):
        return self.sidechain_interaction_model(distance, dpka_max, cutoff, f_angle) # weight is ignored in 3.0 Sep07

    # coulomb
    def electrostatic_interaction(self, group1, group2, distance):
        return self.electrostatic_interaction_model(group1, group2, distance, self)

    def calculateCoulombEnergy(self, distance, weight):
        return self.coulomb_interaction_model(distance, weight, self.parameters)

    def checkCoulombPair(self, group1, group2, distance):
        return self.check_coulomb_pair_method(self.parameters, group1, group2, distance)

    # backbone re-organisation
    def calculateBackBoneReorganization(self, conformation):
        return self.backbone_reorganisation_method(self.parameters, conformation)

    # exceptions
    def checkExceptions(self, group1, group2):
        return self.exception_check_method(self, group1, group2)

    def setup_bonding_and_protonation(self, molecular_container):
        return self.molecular_preparation_method(self.parameters, molecular_container)

    def setup_bonding(self, molecular_container):
        return self.prepare_bonds(self.parameters, molecular_container)



class version_A(version):
    def __init__(self, parameters):
        # set the calculation rutines used in this version
        version.__init__(self, parameters)

        # atom naming, bonding, and protonation
        self.molecular_preparation_method = propka.calculations.setup_bonding_and_protonation
        self.prepare_bonds = propka.calculations.setup_bonding


        # desolvation related methods
        self.desolvation_model = calculations.radial_volume_desolvation
        self.weight_pair_method = calculations.calculatePairWeight

        # side chain methods
        self.sidechain_interaction_model = propka.calculations.HydrogenBondEnergy
        self.hydrogen_bond_interaction_model = propka.calculations.hydrogen_bond_interaction

        # colomb methods
        self.electrostatic_interaction_model = propka.calculations.electrostatic_interaction
        self.check_coulomb_pair_method = propka.calculations.checkCoulombPair
        self.coulomb_interaction_model = propka.calculations.CoulombEnergy

        #backbone
        self.backbone_interaction_model = propka.calculations.HydrogenBondEnergy
        self.backbone_reorganisation_method = propka.calculations.BackBoneReorganization

        # exception methods
        self.exception_check_method = propka.calculations.checkExceptions
        return

    def get_hydrogen_bond_parameters(self, atom1, atom2):
        dpka_max = self.parameters.sidechain_interaction
        cutoff   = self.parameters.sidechain_cutoffs.get_value(atom1.group_type, atom2.group_type)
        return [dpka_max, cutoff]

    def get_backbone_hydrogen_bond_parameters(self, backbone_atom, atom):
        if backbone_atom.group_type == 'BBC':
            if atom.group_type in self.parameters.backbone_CO_hydrogen_bond.keys():
                [v,c1,c2] = self.parameters.backbone_CO_hydrogen_bond[atom.group_type]
                return [v,[c1,c2]]

        if backbone_atom.group_type == 'BBN':
            if atom.group_type in self.parameters.backbone_NH_hydrogen_bond.keys():
                [v,c1,c2] = self.parameters.backbone_NH_hydrogen_bond[atom.group_type]
                return [v,[c1,c2]]

        return None




class simple_hb(version_A):
    def __init__(self, parameters):
        # set the calculation rutines used in this version
        version_A.__init__(self, parameters)
        info('Using simple hb model')
        return

    def get_hydrogen_bond_parameters(self, atom1, atom2):
        return self.parameters.hydrogen_bonds.get_value(atom1.element, atom2.element)


    def get_backbone_hydrogen_bond_parameters(self, backbone_atom, atom):
        return self.parameters.hydrogen_bonds.get_value(backbone_atom.element, atom.element)




class element_based_ligand_interactions(version_A):
    def __init__(self, parameters):
        # set the calculation rutines used in this version
        version_A.__init__(self, parameters)
        info('Using detailed SC model!')
        return

    def get_hydrogen_bond_parameters(self, atom1, atom2):
        if not 'hetatm' in [atom1.type, atom2.type]:
            # this is a protein-protein interaction
            dpka_max = self.parameters.sidechain_interaction.get_value(atom1.group_type, atom2.group_type)
            cutoff   = self.parameters.sidechain_cutoffs.get_value(atom1.group_type, atom2.group_type)
            return [dpka_max, cutoff]

        # at least one ligand atom is involved in this interaction
        # make sure that we are using the heavy atoms for finding paramters
        elements = []
        for a in [atom1, atom2]:
            if a.element == 'H': elements.append(a.bonded_atoms[0].element)
            else: elements.append(a.element)

        return self.parameters.hydrogen_bonds.get_value(elements[0], elements[1])


    def get_backbone_hydrogen_bond_parameters(self, backbone_atom, atom):
        if atom.type == 'atom':
            # this is a backbone-protein interaction
            if backbone_atom.group_type == 'BBC' and\
                    atom.group_type in self.parameters.backbone_CO_hydrogen_bond.keys():
                [v,c1,c2] = self.parameters.backbone_CO_hydrogen_bond[atom.group_type]
                return [v,[c1,c2]]

            if backbone_atom.group_type == 'BBN' and\
                    atom.group_type in self.parameters.backbone_NH_hydrogen_bond.keys():
                [v,c1,c2] = self.parameters.backbone_NH_hydrogen_bond[atom.group_type]
                return [v,[c1,c2]]
        else:
            # this is a backbone-ligand interaction
            # make sure that we are using the heavy atoms for finding paramters
            elements = []
            for a in [backbone_atom, atom]:
                if a.element == 'H': elements.append(a.bonded_atoms[0].element)
                else: elements.append(a.element)

            res = self.parameters.hydrogen_bonds.get_value(elements[0], elements[1])
            if not res:
                info('Could not determine backbone interaction parameters for:',
                     backbone_atom, atom)

            return

        return None



class propka30(version):
    def __init__(self, parameters):
        # set the calculation rutines used in this version
        version.__init__(self, parameters)

        # atom naming, bonding, and protonation
        self.molecular_preparation_method = propka.calculations.setup_bonding_and_protonation_30_style

        # desolvation related methods
        self.desolvation_model = calculations.radial_volume_desolvation
        self.weight_pair_method = calculations.calculatePairWeight

        # side chain methods
        self.sidechain_interaction_model = propka.calculations.HydrogenBondEnergy

        # colomb methods
        self.check_coulomb_pair_method = propka.calculations.checkCoulombPair
        self.coulomb_interaction_model = propka.calculations.CoulombEnergy

        #backbone
        self.backbone_reorganisation_method = propka.calculations.BackBoneReorganization

        # exception methods
        self.exception_check_method = propka.calculations.checkExceptions


        return

    def get_hydrogen_bond_parameters(self, atom1, atom2):
        dpka_max = self.parameters.sidechain_interaction.get_value(atom1.group_type, atom2.group_type)
        cutoff   = self.parameters.sidechain_cutoffs.get_value(atom1.group_type, atom2.group_type)
        return [dpka_max, cutoff]




