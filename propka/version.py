"""
Version-based configuration
===========================

Contains version-specific methods and parameters.

TODO - this module unnecessarily confuses the code.  Can we eliminate it?
"""
import logging
from propka.atom import Atom
from propka.hydrogens import setup_bonding_and_protonation, setup_bonding
from propka.hydrogens import setup_bonding_and_protonation_30_style
from propka.energy import radial_volume_desolvation, calculate_pair_weight
from propka.energy import hydrogen_bond_energy, hydrogen_bond_interaction
from propka.energy import electrostatic_interaction, check_coulomb_pair
from propka.energy import coulomb_energy, check_exceptions
from propka.energy import backbone_reorganization


_LOGGER = logging.getLogger(__name__)


class Version:
    """Store version-specific methods and parameters."""
    def __init__(self, parameters):
        self.parameters = parameters
        self.desolvation_model = self.empty_function
        self.weight_pair_method = self.empty_function
        self.hydrogen_bond_interaction_model = self.empty_function
        self.sidechain_interaction_model = self.empty_function
        self.electrostatic_interaction_model = self.empty_function
        self.coulomb_interaction_model = self.empty_function
        self.check_coulomb_pair_method = self.empty_function
        self.backbone_reorganisation_method = self.empty_function
        self.exception_check_method = self.empty_function
        self.molecular_preparation_method = self.empty_function
        self.prepare_bonds = self.empty_function

    @staticmethod
    def empty_function(*args):
        """Placeholder function so we don't use uninitialized variables.

        Args:
            args:  whatever arguments would have been passed to the function
        Raises:
            NotImplementedError
        """
        err = "Called an empty Version function with args {0:s}".format(args)
        raise NotImplementedError(err)

    def calculate_desolvation(self, group):
        """Calculate desolvation energy using assigned model."""
        return self.desolvation_model(self.parameters, group)

    def calculate_pair_weight(self, num_volume1, num_volume2):
        """Calculate pair weight using assigned model."""
        return self.weight_pair_method(
            self.parameters, num_volume1, num_volume2)

    def hydrogen_bond_interaction(self, group1, group2):
        """Calculate H-bond energy using assigned model."""
        return self.hydrogen_bond_interaction_model(group1, group2, self)

    def calculate_side_chain_energy(self, distance, dpka_max, cutoff, _,
                                    f_angle):
        """Calculate sidechain energy using assigned model."""
        return self.sidechain_interaction_model(
            distance, dpka_max, cutoff, f_angle)

    def electrostatic_interaction(self, group1, group2, distance):
        """Calculate electrostatic energy using assigned model."""
        return self.electrostatic_interaction_model(
            group1, group2, distance, self)

    def calculate_coulomb_energy(self, distance, weight):
        """Calculate Coulomb energy using assigned model."""
        return self.coulomb_interaction_model(
            distance, weight, self.parameters)

    def check_coulomb_pair(self, group1, group2, distance):
        """Check Coulomb pair using assigned model."""
        return self.check_coulomb_pair_method(
            self.parameters, group1, group2, distance)

    def calculate_backbone_reorganization(self, conformation):
        """Calculate backbone reorganization using assigned model."""
        return self.backbone_reorganisation_method(
            self.parameters, conformation)

    def check_exceptions(self, group1, group2):
        """Calculate exceptions using assigned model."""
        return self.exception_check_method(self, group1, group2)

    def setup_bonding_and_protonation(self, molecular_container):
        """Setup bonding and protonation using assigned model."""
        return self.molecular_preparation_method(molecular_container)

    def setup_bonding(self, molecular_container):
        """Setup bonding using assigned model."""
        return self.prepare_bonds(self.parameters, molecular_container)

    def get_hydrogen_bond_parameters(self, atom1: Atom, atom2: Atom) -> tuple:
        """Get hydrogen bond parameters for two atoms."""
        raise NotImplementedError("abstract method")


class VersionA(Version):
    """TODO - figure out what this is."""

    def __init__(self, parameters):
        """Initialize object with parameters."""
        # set the calculation rutines used in this version
        super().__init__(parameters)
        self.molecular_preparation_method = setup_bonding_and_protonation
        self.prepare_bonds = setup_bonding
        self.desolvation_model = radial_volume_desolvation
        self.weight_pair_method = calculate_pair_weight
        self.sidechain_interaction_model = hydrogen_bond_energy
        self.hydrogen_bond_interaction_model = hydrogen_bond_interaction
        self.electrostatic_interaction_model = electrostatic_interaction
        self.check_coulomb_pair_method = check_coulomb_pair
        self.coulomb_interaction_model = coulomb_energy
        self.backbone_interaction_model = hydrogen_bond_energy
        self.backbone_reorganisation_method = backbone_reorganization
        self.exception_check_method = check_exceptions

    def get_hydrogen_bond_parameters(self, atom1, atom2):
        """Get hydrogen bond parameters for two atoms.

        Args:
            atom1:  first atom
            atom2:  second atom
        Returns:
            [dpka_max, cutoff]
        """
        dpka_max = self.parameters.sidechain_interaction
        cutoff = self.parameters.sidechain_cutoffs.get_value(
            atom1.group_type, atom2.group_type)
        return [dpka_max, cutoff]

    def get_backbone_hydrogen_bond_parameters(self, backbone_atom, atom):
        """Get hydrogen bond parameters between backbone atom and other atom.

        Args:
            backbone_atom:  backbone atom
            atom:  other atom
        Returns
            [v, [c1, c3]]  TODO - figure out what this is
        """
        if backbone_atom.group_type == 'BBC':
            if (
                    atom.group_type in
                    self.parameters.backbone_CO_hydrogen_bond.keys()
                    ):
                [v, c1, c2] = self.parameters.backbone_CO_hydrogen_bond[
                    atom.group_type]
                return [v, [c1, c2]]
        if backbone_atom.group_type == 'BBN':
            if (
                    atom.group_type in
                    self.parameters.backbone_NH_hydrogen_bond.keys()
                    ):
                [v, c1, c2] = self.parameters.backbone_NH_hydrogen_bond[
                    atom.group_type]
                return [v, [c1, c2]]
        return None


class SimpleHB(VersionA):
    """A simple hydrogen bond version."""

    def __init__(self, parameters):
        """Initialize object with parameters."""
        # set the calculation rutines used in this version
        super().__init__(parameters)
        _LOGGER.info('Using simple hb model')

    def get_hydrogen_bond_parameters(self, atom1, atom2):
        """Get hydrogen bond parameters for two atoms.

        Args:
            atom1:  first atom
            atom2:  second atom
        Returns:
            [dpka_max, cutoff]
        """
        return self.parameters.hydrogen_bonds.get_value(
            atom1.element, atom2.element)

    def get_backbone_hydrogen_bond_parameters(self, backbone_atom, atom):
        """Get hydrogen bond parameters between backbone atom and other atom.

        Args:
            backbone_atom:  backbone atom
            atom:  other atom
        Returns
            [v, [c1, c3]]  TODO - figure out what this is
        """
        return self.parameters.hydrogen_bonds.get_value(
            backbone_atom.element, atom.element)


class ElementBasedLigandInteractions(VersionA):
    """TODO - figure out what this is."""

    def __init__(self, parameters):
        """Initialize object with parameters."""
        # set the calculation rutines used in this version
        super().__init__(parameters)
        _LOGGER.info('Using detailed SC model!')
        return

    def get_hydrogen_bond_parameters(self, atom1, atom2):
        """Get hydrogen bond parameters for two atoms.

        Args:
            atom1:  first atom
            atom2:  second atom
        Returns:
            [dpka_max, cutoff]
        """
        if 'hetatm' not in [atom1.type, atom2.type]:
            # this is a protein-protein interaction
            dpka_max = self.parameters.sidechain_interaction.get_value(
                atom1.group_type, atom2.group_type)
            cutoff = self.parameters.sidechain_cutoffs.get_value(
                atom1.group_type, atom2.group_type)
            return [dpka_max, cutoff]
        # at least one ligand atom is involved in this interaction
        # make sure that we are using the heavy atoms for finding paramters
        elements = []
        for atom in [atom1, atom2]:
            if atom.element == 'H':
                elements.append(atom.bonded_atoms[0].element)
            else:
                elements.append(atom.element)
        return self.parameters.hydrogen_bonds.get_value(
            elements[0], elements[1])

    def get_backbone_hydrogen_bond_parameters(self, backbone_atom, atom):
        """Get hydrogen bond parameters between backbone atom and other atom.

        Args:
            backbone_atom:  backbone atom
            atom:  other atom
        Returns
            [v, [c1, c3]]  TODO - figure out what this is
        """
        if atom.type == 'atom':
            # this is a backbone-protein interaction
            if (backbone_atom.group_type == 'BBC'
                    and atom.group_type
                    in self.parameters.backbone_CO_hydrogen_bond.keys()):
                [v, c1, c2] = self.parameters.backbone_CO_hydrogen_bond[
                    atom.group_type]
                return [v, [c1, c2]]

            if (backbone_atom.group_type == 'BBN'
                    and atom.group_type
                    in self.parameters.backbone_NH_hydrogen_bond.keys()):
                [v, c1, c2] = self.parameters.backbone_NH_hydrogen_bond[
                    atom.group_type]
                return [v, [c1, c2]]
        else:
            # this is a backbone-ligand interaction
            # make sure that we are using the heavy atoms for finding paramters
            elements = []
            for atom2 in [backbone_atom, atom]:
                if atom2.element == 'H':
                    elements.append(atom2.bonded_atoms[0].element)
                else:
                    elements.append(atom2.element)
            res = self.parameters.hydrogen_bonds.get_value(
                elements[0], elements[1])
            if not res:
                _LOGGER.info(
                    'Could not determine backbone interaction parameters '
                    'for: %s %s',
                    backbone_atom, atom)
            return None
        return None


class Propka30(Version):
    """Version class for PROPKA 3.0."""

    def __init__(self, parameters):
        """Initialize object with parameters."""
        # set the calculation routines used in this version
        super().__init__(parameters)
        self.molecular_preparation_method = (
            setup_bonding_and_protonation_30_style)
        self.desolvation_model = radial_volume_desolvation
        self.weight_pair_method = calculate_pair_weight
        self.sidechain_interaction_model = hydrogen_bond_energy
        self.check_coulomb_pair_method = check_coulomb_pair
        self.coulomb_interaction_model = coulomb_energy
        self.backbone_reorganisation_method = backbone_reorganization
        self.exception_check_method = check_exceptions

    def get_hydrogen_bond_parameters(self, atom1, atom2):
        """Get hydrogen bond parameters for two atoms.

        Args:
            atom1:  first atom
            atom2:  second atom
        Returns:
            [dpka_max, cutoff]
        """
        dpka_max = self.parameters.sidechain_interaction.get_value(
            atom1.group_type, atom2.group_type)
        cutoff = self.parameters.sidechain_cutoffs.get_value(
            atom1.group_type, atom2.group_type)
        return [dpka_max, cutoff]
