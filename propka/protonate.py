"""Protonate a structure."""
import math
import propka.bonds
import propka.atom
from propka.vector_algebra import rotate_vector_around_an_axis, Vector
from propka.lib import warning, debug


class Protonate:
    """ Protonates atoms using VSEPR theory """

    def __init__(self, verbose=False):
        """Initialize with flag for verbosity

        Args:
            verbose:  True for verbose output
        """
        self.verbose = verbose
        self.valence_electrons = {'H': 1, 'He': 2, 'Li': 1, 'Be': 2, 'B': 3,
                                  'C': 4, 'N': 5, 'O': 6, 'F': 7, 'Ne': 8,
                                  'Na': 1, 'Mg': 2, 'Al': 3, 'Si': 4, 'P': 5,
                                  'S': 6, 'Cl': 7, 'Ar': 8, 'K': 1, 'Ca': 2,
                                  'Sc': 2, 'Ti': 2, 'Va': 2, 'Cr': 1, 'Mn': 2,
                                  'Fe': 2, 'Co': 2, 'Ni': 2, 'Cu': 1, 'Zn': 2,
                                  'Ga': 3, 'Ge': 4, 'As': 5, 'Se': 6, 'Br': 7,
                                  'Kr': 8, 'I': 7}
        # TODO - consider putting charges in a configuration file
        self.standard_charges = {'ARG-NH1': 1.0, 'ASP-OD2': -1.0,
                                 'GLU-OE2': -1.0, 'HIS-ND1': 1.0,
                                 'LYS-NZ': 1.0, 'N+': 1.0, 'C-': -1.0}
        self.sybyl_charges = {'N.pl3': 1, 'N.3': 1, 'N.4': 1, 'N.ar': 1,
                              'O.co2-': 1}
        # TODO - consider putting bond lengths in a configuration file
        self.bond_lengths = {'C': 1.09, 'N': 1.01, 'O': 0.96, 'F': 0.92,
                             'Cl': 1.27, 'Br': 1.41, 'I': 1.61, 'S': 1.35}
        self.protonation_methods = {4: self.tetrahedral, 3: self.trigonal}

    def protonate(self, molecules):
        """Protonate all atoms in the molecular container.

        Args:
            molecules:  molecular containers
        """
        debug('----- Protonation started -----')
        # Remove all currently present hydrogen atoms
        self.remove_all_hydrogen_atoms(molecules)
        # protonate all atoms
        for name in molecules.conformation_names:
            non_h_atoms = (molecules.conformations[name]
                           .get_non_hydrogen_atoms())
            for atom in non_h_atoms:
                self.protonate_atom(atom)

    @staticmethod
    def remove_all_hydrogen_atoms(molecular_container):
        """Remove all hydrogen atoms from molecule.

        Args:
            molecular_container:  molecule to remove hydrogens from
        """
        for name in molecular_container.conformation_names:
            molecular_container.conformations[name].atoms = (
                molecular_container.conformations[name]
                .get_non_hydrogen_atoms())

    def set_charge(self, atom):
        """Set charge for atom.

        Args:
            atom:  atom to be charged
        """
        # atom is a protein atom
        if atom.type == 'atom':
            key = '%3s-%s' % (atom.res_name, atom.name)
            if atom.terminal:
                debug(atom.terminal)
                key = atom.terminal
            if key in self.standard_charges:
                atom.charge = self.standard_charges[key]
                debug('Charge', atom, atom.charge)
                atom.charge_set = True
        # atom is a ligand atom
        elif atom.type == 'hetatm':
            if atom.sybyl_type in self.sybyl_charges:
                atom.charge = self.sybyl_charges[atom.sybyl_type]
                atom.sybyl_type = atom.sybyl_type.replace('-', '')
                atom.charge_set = True

    def protonate_atom(self, atom):
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
    def set_proton_names(heavy_atoms):
        """Set names for protons.

        Args:
            heavy_atoms:  list of heavy atoms with protons to be named
        """
        for heavy_atom in heavy_atoms:
            i = 1
            for bonded in heavy_atom.bonded_atoms:
                if bonded.element == 'H':
                    bonded.name += '%d' % i
                    i += 1

    def set_number_of_protons_to_add(self, atom):
        """Set the number of protons to add to this atom.

        Args:
            atom:  atom for calculation
        """
        debug('*'*10)
        debug('Setting number of protons to add for', atom)
        atom.number_of_protons_to_add = 8
        debug('                  %4d' % 8)
        atom.number_of_protons_to_add -= self.valence_electrons[atom.element]
        debug('Valence eletrons: %4d' % -self.valence_electrons[atom.element])
        atom.number_of_protons_to_add -= len(atom.bonded_atoms)
        debug('Number of bonds:  %4d' % -len(atom.bonded_atoms))
        atom.number_of_protons_to_add -= atom.num_pi_elec_2_3_bonds
        debug('Pi electrons:     %4d' % -atom.num_pi_elec_2_3_bonds)
        atom.number_of_protons_to_add += int(atom.charge)
        debug('Charge:           %4.1f' % atom.charge)
        debug('-'*10)
        debug(atom.number_of_protons_to_add)

    def set_steric_number_and_lone_pairs(self, atom):
        """Set steric number and lone pairs for atom.

        Args:
            atom:  atom for calculation
        """
        # If we already did this, there is no reason to do it again
        if atom.steric_num_lone_pairs_set:
            return
        debug('='*10)
        debug('Setting steric number and lone pairs for', atom)
        atom.steric_number = 0
        debug('%65s: %4d' % ('Valence electrons',
                             self.valence_electrons[atom.element]))
        atom.steric_number += self.valence_electrons[atom.element]
        debug('%65s: %4d' % ('Number of bonds',
                             len(atom.bonded_atoms)))
        atom.steric_number += len(atom.bonded_atoms)
        debug('%65s: %4d' % ('Number of hydrogen atoms to add',
                             atom.number_of_protons_to_add))
        atom.steric_number += atom.number_of_protons_to_add
        debug('%65s: %4d' % ('Number of pi-electrons in double '
                             'and triple bonds(-)',
                             atom.num_pi_elec_2_3_bonds))
        atom.steric_number -= atom.num_pi_elec_2_3_bonds
        debug('%65s: %4d' % ('Number of pi-electrons in conjugated double and '
                             'triple bonds(-)',
                             atom.num_pi_elec_conj_2_3_bonds))
        atom.steric_number -= atom.num_pi_elec_conj_2_3_bonds
        debug('%65s: %4d' % ('Number of donated co-ordinated bonds', 0))
        atom.steric_number += 0
        debug('%65s: %4.1f' % ('Charge(-)', atom.charge))
        atom.steric_number -= atom.charge
        atom.steric_number = math.floor(atom.steric_number/2.0)
        atom.number_of_lone_pairs = (atom.steric_number
                                     - len(atom.bonded_atoms)
                                     - atom.number_of_protons_to_add)
        debug('-'*70)
        debug('%65s: %4d' % ('Steric number', atom.steric_number))
        debug('%65s: %4d' % ('Number of lone pairs',
                             atom.number_of_lone_pairs))
        atom.steric_num_lone_pairs_set = True

    def add_protons(self, atom):
        """Add protons to atom.

        Args:
            atom:  atom for calculation
        """
        # decide which method to use
        debug('PROTONATING', atom)
        if atom.steric_number in list(self.protonation_methods.keys()):
            self.protonation_methods[atom.steric_number](atom)
        else:
            warning('Do not have a method for protonating',
                    atom, '(steric number: %d)' % atom.steric_number)

    def trigonal(self, atom):
        """Add hydrogens in trigonal geometry.

        Args:
            atom:  atom to protonate
        """
        debug('TRIGONAL - %d bonded atoms' % len(atom.bonded_atoms))
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
                for i, bonded_atom in enumerate(atom.bonded_atoms[0].bonded_atoms):
                    if bonded_atom != atom:
                        other_atom_indices.append(i)
                vec1 = Vector(atom1=atom, atom2=atom.bonded_atoms[0])
                vec2 = Vector(atom1=atom.bonded_atoms[0],
                              atom2=atom.bonded_atoms[0]
                              .bonded_atoms[other_atom_indices[0]])
                axis = vec1**vec2
                # this is a trick to make sure that the order of atoms doesn't
                # influence the final postions of added protons
                if len(other_atom_indices) > 1:
                    vec3 = Vector(atom1=atom.bonded_atoms[0],
                                  atom2=atom.bonded_atoms[0]
                                  .bonded_atoms[other_atom_indices[1]])
                    axis2 = vec1**vec3
                    if axis*axis2 > 0:
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

    def tetrahedral(self, atom):
        """Protonate atom in tetrahedral geometry.

        Args:
            atom:  atom to protonate.
        """
        debug('TETRAHEDRAL - %d bonded atoms' % len(atom.bonded_atoms))
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
    def add_proton(atom, position):
        """Add a proton to an atom at a specific position.

        Args:
            atom:  atom to protonate
            position:  position for proton
        """
        # Create the new proton
        new_h = propka.atom.Atom()
        new_h.set_property(
            numb=None,
            name='H%s' % atom.name[1:],
            res_name=atom.res_name,
            chain_id=atom.chain_id,
            res_num=atom.res_num,
            x=round(position.x, 3), # round of to three decimal points to
                                    # avoid round-off differences in input file
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
        new_h.is_protonates = True
        atom.bonded_atoms.append(new_h)
        atom.number_of_protons_to_add -= 1
        atom.conformation_container.add_atom(new_h)
        # update names of all protons on this atom
        new_h.residue_label = "%-3s%4d%2s" % (new_h.name, new_h.res_num,
                                              new_h.chain_id)
        no_protons = atom.count_bonded_elements('H')
        if no_protons > 1:
            i = 1
            for proton in atom.get_bonded_elements('H'):
                proton.name = 'H%s%d' % (atom.name[1:], i)
                proton.residue_label = "%-3s%4d%2s" % (proton.name,
                                                       proton.res_num,
                                                       proton.chain_id)
                i += 1
        debug('added', new_h, 'to', atom)

    def set_bond_distance(self, bvec, element):
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
            str_ = ('Bond length for %s not found, using the standard value '
                    'of %f' % (element, dist))
            warning(str_)
        bvec = bvec.rescale(dist)
        return bvec
