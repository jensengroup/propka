"""
Protonate a structure
=====================

The :class:`Protonate` processes a
:class:`propka.molecular_container.MolecularContainer` and adds
protons.

"""
import logging
import math
from typing import Iterable, TYPE_CHECKING

import propka.bonds
import propka.atom
from propka.atom import Atom
from propka.vector_algebra import rotate_vector_around_an_axis, Vector

if TYPE_CHECKING:
    from propka.molecular_container import MolecularContainer


_LOGGER = logging.getLogger(__name__)


class Protonate:
    """ Protonates atoms using VSEPR theory """

    def __init__(self, verbose=False):
        """Initialize with flag for verbosity

        Args:
            verbose:  True for verbose output
        """
        self.verbose = verbose
        self.valence_electrons = {
            "H": 1, "He": 2, "Li": 1, "Be": 2, "B": 3, "C": 4, "N": 5,
            "O": 6, "F": 7, "Ne": 8, "Na": 1, "Mg": 2, "Al": 3, "Si": 4,
            "P": 5, "S": 6, "Cl": 7, "Ar": 8, "K": 1, "Ca": 2, "Sc": 2,
            "Ti": 2, "V": 2, "Cr": 1, "Mn": 2, "Fe": 2, "Co": 2, "Ni": 2,
            "Cu": 1, "Zn": 2, "Ga": 3, "Ge": 4, "As": 5, "Se": 6, "Br": 7,
            "Kr": 8, "Rb": 1, "Sr": 2, "Y": 2, "Zr": 2, "Nb": 1, "Mo": 1,
            "Tc": 2, "Ru": 1, "Rh": 1, "Pd": 8, "Ag": 1, "Cd": 2, "In": 3,
            "Sn": 4, "Sb": 5, "Te": 6, "I": 7, "Xe": 8, "Cs": 1, "Ba": 2,
            "La": 2, "Ce": 2, "Pr": 2, "Nd": 2, "Pm": 2, "Sm": 2, "Eu": 2,
            "Gd": 2, "Tb": 2, "Dy": 2, "Ho": 2, "Er": 2, "Tm": 2, "Yb": 2,
            "Lu": 2, "Hf": 2, "Ta": 2, "W": 2, "Re": 2, "Os": 2, "Ir": 2,
            "Pt": 1, "Au": 1, "Hg": 2, "Tl": 3, "Pb": 4, "Bi": 5, "Po": 6,
            "At": 7, "Rn": 8, "Fr": 1, "Ra": 2, "Ac": 2, "Th": 2, "Pa": 2,
            "U": 2, "Np": 2, "Pu": 2, "Am": 2, "Cm": 2, "Bk": 2,  "Cf": 2,
            "Es": 2, "Fm": 2, "Md": 2, "No": 2, "Lr": 3, "Rf": 2, "Db": 2,
            "Sg": 2, "Bh": 2, "Hs": 2, "Mt": 2, "Ds": 2, "Rg": 2, "Cn": 2,
            "Nh": 3, "Fl": 4, "Mc": 5, "Lv": 6, "Ts": 7, "Og": 8, "Uue": 1
            }
        # TODO - consider putting charges in a configuration file
        self.standard_charges = {
            'ARG-NH1': 1.0, 'ASP-OD2': -1.0, 'GLU-OE2': -1.0, 'HIS-ND1': 1.0,
            'LYS-NZ': 1.0, 'N+': 1.0, 'C-': -1.0}
        self.sybyl_charges = {
            'N.pl3': 1, 'N.3': 1, 'N.4': 1, 'N.ar': 1, 'O.co2-': 1}
        # TODO - consider putting bond lengths in a configuration file
        self.bond_lengths = {
            'C': 1.09, 'N': 1.01, 'O': 0.96, 'F': 0.92, 'Cl': 1.27,
            'Br': 1.41, 'I': 1.61, 'S': 1.35}
        self.protonation_methods = {4: self.tetrahedral, 3: self.trigonal}

    def protonate(self, molecules: "MolecularContainer"):
        """Protonate all atoms in the molecular container.

        Args:
            molecules:  molecular containers
        """
        _LOGGER.debug('----- Protonation started -----')
        # Remove all currently present hydrogen atoms
        self.remove_all_hydrogen_atoms(molecules)
        # protonate all atoms
        for name in molecules.conformation_names:
            non_h_atoms = (molecules.conformations[name]
                           .get_non_hydrogen_atoms())
            for atom in non_h_atoms:
                self.protonate_atom(atom)

    @staticmethod
    def remove_all_hydrogen_atoms(molecular_container: "MolecularContainer"):
        """Remove all hydrogen atoms from molecule.

        Args:
            molecular_container:  molecule to remove hydrogens from
        """
        for name in molecular_container.conformation_names:
            molecular_container.conformations[name].atoms = (
                molecular_container.conformations[name]
                .get_non_hydrogen_atoms())

    def set_charge(self, atom: Atom):
        """Set charge for atom.

        Args:
            atom:  atom to be charged
        """
        # atom is a protein atom
        if atom.type == 'atom':
            key = '{0:3s}-{1:s}'.format(atom.res_name, atom.name)
            if atom.terminal:
                _LOGGER.debug("%s", atom.terminal)
                key = atom.terminal
            if key in self.standard_charges:
                atom.charge = self.standard_charges[key]
                _LOGGER.debug('Charge %s %s', atom, atom.charge)
                atom.charge_set = True
        # atom is a ligand atom
        elif atom.type == 'hetatm':
            if atom.sybyl_type in self.sybyl_charges:
                atom.charge = self.sybyl_charges[atom.sybyl_type]
                atom.sybyl_type = atom.sybyl_type.replace('-', '')
                atom.charge_set = True

    def protonate_atom(self, atom: Atom):
        """Protonate an atom.

        Args:
            atom:  atom to be protonated
        """
        if atom.is_protonated:
            return
        if atom.element == 'H':
            return
        self.set_charge(atom)
        self.set_number_of_protons_to_add(atom)
        self.set_steric_number_and_lone_pairs(atom)
        self.add_protons(atom)
        atom.is_protonated = True

    @staticmethod
    def set_proton_names(heavy_atoms: Iterable[Atom]):
        """Set names for protons.

        Args:
            heavy_atoms:  list of heavy atoms with protons to be named
        """
        for heavy_atom in heavy_atoms:
            i = 1
            for bonded in heavy_atom.bonded_atoms:
                if bonded.element == 'H':
                    bonded.name += str(i)
                    i += 1

    def set_number_of_protons_to_add(self, atom: Atom):
        """Set the number of protons to add to this atom.

        Args:
            atom:  atom for calculation
        """
        _LOGGER.debug('*'*10)
        _LOGGER.debug('Setting number of protons to add for %s', atom)
        atom.number_of_protons_to_add = 8
        _LOGGER.debug("                     8")
        if atom.element not in self.valence_electrons:
            _LOGGER.warning(
                    f'Unknown valence electron for element {atom.element}')
            self.valence_electrons[atom.element] = 4
        atom.number_of_protons_to_add -= self.valence_electrons[atom.element]
        _LOGGER.debug('Valence electrons: {0:>4d}'.format(
            -self.valence_electrons[atom.element]))
        atom.number_of_protons_to_add -= len(atom.bonded_atoms)
        _LOGGER.debug(
            'Number of bonds:  {0:>4d}'.format(-len(atom.bonded_atoms))
        )
        atom.number_of_protons_to_add -= atom.num_pi_elec_2_3_bonds
        _LOGGER.debug(
            'Pi electrons:     {0:>4d}'.format(-atom.num_pi_elec_2_3_bonds)
        )
        atom.number_of_protons_to_add += int(atom.charge)
        _LOGGER.debug('Charge:           {0:>4.1f}'.format(atom.charge))
        _LOGGER.debug('-'*10)
        _LOGGER.debug(atom.number_of_protons_to_add)

    def set_steric_number_and_lone_pairs(self, atom: Atom):
        """Set steric number and lone pairs for atom.

        Args:
            atom:  atom for calculation
        """
        # If we already did this, there is no reason to do it again
        if atom.steric_num_lone_pairs_set:
            return
        _LOGGER.debug('='*10)
        _LOGGER.debug('Setting steric number and lone pairs for %s', atom)
        atom.steric_number = 0
        if atom.element not in self.valence_electrons:
            self.valence_electrons[atom.element] = 4
            _LOGGER.warning(
                    f"Not found valence for element {atom.element}, use 4")
        _LOGGER.debug('{0:>65s}: {1:>4d}'.format(
            'Valence electrons', self.valence_electrons[atom.element]))
        atom.steric_number += self.valence_electrons[atom.element]
        _LOGGER.debug('{0:>65s}: {1:>4d}'.format(
            'Number of bonds', len(atom.bonded_atoms)))
        atom.steric_number += len(atom.bonded_atoms)
        _LOGGER.debug('{0:>65s}: {1:>4d}'.format(
            'Number of hydrogen atoms to add', atom.number_of_protons_to_add))
        atom.steric_number += atom.number_of_protons_to_add
        _LOGGER.debug('{0:>65s}: {1:>4d}'.format(
            'Number of pi-electrons in double and triple bonds(-)',
            atom.num_pi_elec_2_3_bonds))
        atom.steric_number -= atom.num_pi_elec_2_3_bonds
        _LOGGER.debug('{0:>65s}: {1:>4d}'.format(
            'Number of pi-electrons in conjugated double and triple bonds(-)',
            atom.num_pi_elec_conj_2_3_bonds))
        atom.steric_number -= atom.num_pi_elec_conj_2_3_bonds
        _LOGGER.debug('{0:>65s}: {1:>4d}'.format(
            'Number of donated co-ordinated bonds', 0))
        atom.steric_number += 0
        _LOGGER.debug('{0:>65s}: {1:>4.1f}'.format(
            'Charge(-)', atom.charge))
        atom.steric_number = math.floor((atom.steric_number - atom.charge) / 2)
        atom.number_of_lone_pairs = (
            atom.steric_number - len(atom.bonded_atoms)
            - atom.number_of_protons_to_add
        )
        _LOGGER.debug('-'*70)
        _LOGGER.debug('{0:>65s}: {1:>4d}'.format(
            'Steric number', atom.steric_number))
        _LOGGER.debug('{0:>65s}: {1:>4d}'.format(
            'Number of lone pairs', atom.number_of_lone_pairs))
        atom.steric_num_lone_pairs_set = True

    def add_protons(self, atom: Atom):
        """Add protons to atom.

        Args:
            atom:  atom for calculation
        """
        # decide which method to use
        _LOGGER.debug('PROTONATING %s', atom)
        if atom.steric_number in list(self.protonation_methods.keys()):
            self.protonation_methods[atom.steric_number](atom)
        else:
            _LOGGER.warning(
                'Do not have a method for protonating %s %s', atom,
                '(steric number: {0:d})'.format(atom.steric_number)
            )

    def trigonal(self, atom: Atom):
        """Add hydrogens in trigonal geometry.

        Args:
            atom:  atom to protonate
        """
        _LOGGER.debug(
            'TRIGONAL - {0:d} bonded atoms'.format(len(atom.bonded_atoms))
        )
        rot_angle = math.radians(120.0)
        cvec = Vector(atom1=atom)
        # 0 bonds
        if len(atom.bonded_atoms) == 0:
            pass
        # 1 bond
        if len(atom.bonded_atoms) == 1 and atom.number_of_protons_to_add > 0:
            # Add another atom with the right angle to the first one
            avec = Vector(atom1=atom, atom2=atom.bonded_atoms[0])
            # use plane of bonded trigonal atom - e.g. arg
            self.set_steric_number_and_lone_pairs(atom.bonded_atoms[0])
            if (atom.bonded_atoms[0].steric_number == 3
                    and len(atom.bonded_atoms[0].bonded_atoms) > 1):
                # use other atoms bonded to the neighbour to establish the
                # plane, if possible
                other_atom_indices = []
                for i, bonded_atom in enumerate(
                        atom.bonded_atoms[0].bonded_atoms):
                    if bonded_atom != atom:
                        other_atom_indices.append(i)
                vec1 = Vector(atom1=atom, atom2=atom.bonded_atoms[0])
                vec2 = Vector(atom1=atom.bonded_atoms[0],
                              atom2=atom.bonded_atoms[0]
                              .bonded_atoms[other_atom_indices[0]])
                axis = vec1.cross(vec2)
                # this is a trick to make sure that the order of atoms doesn't
                # influence the final postions of added protons
                if len(other_atom_indices) > 1:
                    vec3 = Vector(atom1=atom.bonded_atoms[0],
                                  atom2=atom.bonded_atoms[0]
                                  .bonded_atoms[other_atom_indices[1]])
                    axis2 = vec1.cross(vec3)
                    if axis.dot(axis2) > 0:
                        axis = axis+axis2
                    else:
                        axis = axis-axis2
            else:
                axis = avec.orthogonal()
            avec = rotate_vector_around_an_axis(rot_angle, axis, avec)
            avec = self.set_bond_distance(avec, atom.element)
            self.add_proton(atom, cvec+avec)
        # 2 bonds
        if len(atom.bonded_atoms) == 2 and atom.number_of_protons_to_add > 0:
            # Add another atom with the right angle to the first two
            avec1 = Vector(atom1=atom, atom2=atom.bonded_atoms[0]).rescale(1.0)
            avec2 = Vector(atom1=atom, atom2=atom.bonded_atoms[1]).rescale(1.0)

            new_a = -avec1 - avec2
            new_a = self.set_bond_distance(new_a, atom.element)
            self.add_proton(atom, cvec+new_a)

    def tetrahedral(self, atom: Atom):
        """Protonate atom in tetrahedral geometry.

        Args:
            atom:  atom to protonate.
        """
        _LOGGER.debug(
            'TETRAHEDRAL - {0:d} bonded atoms'.format(len(atom.bonded_atoms)))
        # TODO - might be good to move tetrahedral angle to constant
        rot_angle = math.radians(109.5)
        cvec = Vector(atom1=atom)
        # 0 bonds
        if len(atom.bonded_atoms) == 0:
            pass
        # 1 bond
        if len(atom.bonded_atoms) == 1 and atom.number_of_protons_to_add > 0:
            # Add another atom with the right angle to the first one
            avec = Vector(atom1=atom, atom2=atom.bonded_atoms[0])
            axis = avec.orthogonal()
            avec = rotate_vector_around_an_axis(rot_angle, axis, avec)
            avec = self.set_bond_distance(avec, atom.element)
            self.add_proton(atom, cvec+avec)
        # 2 bonds
        if len(atom.bonded_atoms) == 2 and atom.number_of_protons_to_add > 0:
            # Add another atom with the right angle to the first two
            avec1 = Vector(atom1=atom, atom2=atom.bonded_atoms[0]).rescale(1.0)
            avec2 = Vector(atom1=atom, atom2=atom.bonded_atoms[1]).rescale(1.0)
            axis = avec1 + avec2
            new_a = rotate_vector_around_an_axis(math.radians(90), axis,
                                                 -avec1)
            new_a = self.set_bond_distance(new_a, atom.element)
            self.add_proton(atom, cvec+new_a)
        # 3 bonds
        if len(atom.bonded_atoms) == 3 and atom.number_of_protons_to_add > 0:
            avec1 = Vector(atom1=atom, atom2=atom.bonded_atoms[0]).rescale(1.0)
            avec2 = Vector(atom1=atom, atom2=atom.bonded_atoms[1]).rescale(1.0)
            avec3 = Vector(atom1=atom, atom2=atom.bonded_atoms[2]).rescale(1.0)
            new_a = -avec1-avec2-avec3
            new_a = self.set_bond_distance(new_a, atom.element)
            self.add_proton(atom, cvec+new_a)

    @staticmethod
    def add_proton(atom: Atom, position: Vector):
        """Add a proton to an atom at a specific position.

        Args:
            atom:  atom to protonate
            position:  position for proton
        """
        assert atom.conformation_container is not None
        # Create the new proton
        new_h = propka.atom.Atom()
        new_h.set_property(
            numb=None,
            name='H{0:s}'.format(atom.name[1:]),
            res_name=atom.res_name,
            chain_id=atom.chain_id,
            res_num=atom.res_num,
            x=round(position.x, 3),  # round of to three decimal points to
                                     # avoid round-off differences in input
                                     # file
            y=round(position.y, 3),
            z=round(position.z, 3),
            occ=None,
            beta=None)
        new_h.element = 'H'
        new_h.type = atom.type
        new_h.bonded_atoms = [atom]
        new_h.charge = 0
        new_h.steric_number = 0
        new_h.number_of_lone_pairs = 0
        new_h.number_of_protons_to_add = 0
        new_h.num_pi_elec_2_3_bonds = 0
        atom.bonded_atoms.append(new_h)
        atom.number_of_protons_to_add -= 1
        atom.conformation_container.add_atom(new_h)
        # update names of all protons on this atom
        new_h.residue_label = "{0:<3s}{1:>4d}{2:>2s}".format(
            new_h.name, new_h.res_num, new_h.chain_id)
        no_protons = atom.count_bonded_elements('H')
        if no_protons > 1:
            i = 1
            for proton in atom.get_bonded_elements('H'):
                proton.name = 'H{0:s}{1:d}'.format(
                    atom.name[1:], i)
                proton.residue_label = "{0:<3s}{1:>4d}{2:>2s}".format(
                    proton.name, proton.res_num, proton.chain_id)
                i += 1
        _LOGGER.debug('added %s %s %s', new_h, 'to', atom)

    def set_bond_distance(self, bvec: Vector, element: str) -> Vector:
        """Set bond distance between atom and element.

        Args:
            bvec:  bond vector
            element:  bonded element
        Returns:
            scaled bond vector
        """
        dist = 1.0
        if element in list(self.bond_lengths.keys()):
            dist = self.bond_lengths[element]
        else:
            str_ = (
                'Bond length for {0:s} not found, using the standard value '
                'of {1:f}'.format(element, dist))
            _LOGGER.warning(str_)
        bvec = bvec.rescale(dist)
        return bvec
