"""Atom class - contains all atom information found in the PDB file"""
import string
import propka.lib
import propka.group
from . import hybrid36


class Atom(object):
    """Atom class - contains all atom information found in the PDB file"""

    def __init__(self, line=None, verbose=False):
        """Initialize Atom object.

        Args:
            line:  Line from a PDB file to set properties of atom.
            verbose:  TODO - this does not appear to be used.  Can we remove it?
        """
        self.occ = None
        self.numb = None
        self.res_name = None
        self.type = None
        self.chain_id = None
        self.beta = None
        self.icode = None
        self.res_num = None
        self.name = None
        self.element = None
        self.x = None
        self.y = None
        self.z = None
        self.group = None
        self.group_type = None
        self.number_of_bonded_elements = {}
        self.cysteine_bridge = False
        self.bonded_atoms = []
        self.residue = None
        self.conformation_container = None
        self.molecular_container = None
        self.is_protonated = False
        self.steric_num_lone_pairs_set = False
        self.terminal = None
        self.charge = 0
        self.charge_set = False
        self.steric_number = 0
        self.number_of_lone_pairs = 0
        self.number_of_protons_to_add = 0
        self.num_pi_elec_2_3_bonds = 0
        self.num_pi_elec_conj_2_3_bonds = 0
        self.groups_extracted = 0
        self.set_properties(line)
        self.residue_label = "%-3s%4d%2s" % (self.name, self.res_num, self.chain_id)

        # ligand atom types
        self.sybyl_type = ''
        self.sybyl_assigned = False
        self.marvin_pka = False

    def set_properties(self, line):
        """Line from PDB file to set properties of atom.

        Args:
            line:  PDB file line
        """
        self.name = ''
        self.numb = 0
        self.x = 0.0
        self.y = 0.0
        self.z = 0.0
        self.res_num = 0
        self.res_name = ''
        self.chain_id = 'A'
        self.type = ''
        self.occ = '1.0'
        self.beta = '0.0'
        self.element = ''
        self.icode = ''

        if line:
            self.name = line[12:16].strip()
            self.numb = int(hybrid36.decode(line[6:11]))
            self.x = float(line[30:38].strip())
            self.y = float(line[38:46].strip())
            self.z = float(line[46:54].strip())
            self.res_num = int(line[22:26].strip())
            self.res_name = "%-3s" % (line[17:20].strip())
            self.chain_id = line[21]
            # Set chain id to "_" if it is just white space.
            if not self.chain_id.strip():
                self.chain_id = '_'
            self.type = line[:6].strip().lower()

            if self.res_name in ['DA ', 'DC ', 'DG ', 'DT ']:
                self.type = 'hetatm'

            self.occ = line[55:60].strip()
            self.beta = line[60:66].strip()
            self.icode = line[26:27]

            # Set the element using the position of the name in the pdb file
            self.element = line[12:14].strip().strip(string.digits)
            if len(self.name) == 4:
                self.element = self.element[0]
            if len(self.element) == 2:
                self.element = '%1s%1s' % (self.element[0], self.element[1].lower())

    def set_group_type(self, type_):
        """Set group type of atom.

        Args:
            type_:  group type of atom
        """
        self.group_type = type_

    def count_bonded_elements(self, element):
        """Count number of bonded atoms with same element.

        Args:
            element:  element type for test.
        Returns:
            number of bonded atoms.
        """
        return len(self.get_bonded_elements(element))

    def get_bonded_elements(self, element):
        """Get bonded atoms with same element.

        Args:
            element:  element type for test.
        Returns:
            array of bonded atoms.
        """
        res = []
        for ba in self.bonded_atoms:
            if ba.element == element:
                res.append(ba)
        return res

    def get_bonded_heavy_atoms(self):
        """Get the atoms bonded to this one that aren't hydrogen.

        Returns:
            list of atoms.
        """
        array = [ba for ba in self.bonded_atoms if ba.element != 'H']
        return array

    def is_atom_within_bond_distance(self, other_atom, max_bonds, cur_bond):
        """Check if <other_atom> is found within <max_bonds> bonds of self.

        Args:
            other_atom:  atom to check
            max_bonds:  number of bonds to check for other atom bonding to self
        Returns:
            Boolean for atom bond distance
        """
        for ba in self.bonded_atoms:
            if ba == other_atom:
                return True
            if max_bonds > cur_bond:
                if ba.is_atom_within_bond_distance(other_atom, max_bonds, cur_bond+1):
                    return True
        return False

    def set_property(self, numb=None, name=None, res_name=None, chain_id=None,
                     res_num=None, x=None, y=None, z=None, occ=None, beta=None):
        """Set properties of the atom object.

        Args:
            numb:  Atom number
            name:  Atom name
            res_name:  residue name
            chain_id:  chain ID
            res_num:  residue number
            x:  atom x-coordinate
            y:  atom y-coordinate
            z:  atom z-coordinate
            occ:  atom occupancy
            beta:  atom temperature factor
        """
        if numb is not None:
            self.numb = numb
        if name is not None:
            self.name = name
        if res_name is not None:
            self.res_name = res_name
        if chain_id is not None:
            self.chain_id = chain_id
        if res_num is not None:
            self.res_num = res_num
        if x is not None:
            self.x = x
        if y is not None:
            self.y = y
        if z is not None:
            self.z = z
        if occ is not None:
            self.occ = occ
        if beta is not None:
            self.beta = beta

    def make_copy(self):
        """Make a copy of this atom.

        TODO - this could be a property method/attribute

        Returns:
            Another atom object copy of this one.
        """
        new_atom = Atom()
        new_atom.type = self.type
        new_atom.numb = self.numb
        new_atom.name = self.name
        new_atom.element = self.element
        new_atom.res_name = self.res_name
        new_atom.res_num = self.res_num
        new_atom.chain_id = self.chain_id
        new_atom.x = self.x
        new_atom.y = self.y
        new_atom.z = self.z
        new_atom.occ = self.occ
        new_atom.beta = self.beta
        new_atom.terminal = self.terminal
        new_atom.residue_label = self.residue_label
        new_atom.icode = self.icode
        return new_atom

    def make_input_line(self):
        """PDB line for this atom.

        TODO - Could be @property method/attribute
        TODO - figure out difference between make_pdb_line, make_input_line, and make_pdb_line2

        Returns:
            String with PDB-format line.
        """
        group = '-'
        model_pka = '-'
        if self.group:
            group = self.group.type
            if self.terminal == 'C-':
                group = 'C-' ## circumventing C-/COO parameter unification

            if self.group.titratable:
                model_pka = '%6.2f'%self.group.model_pka
        str_ = "%-6s%5d %s " % (self.type.upper(), self.numb,
                                propka.lib.makeTidyAtomLabel(self.name, self.element))
        str_ += "%s%2s%4d%12.3lf%8.3lf%8.3lf%6s%6s \n" % (self.res_name, self.chain_id,
                                                          self.res_num, self.x, self.y,
                                                          self.z, group, model_pka)
        return str_

    def make_conect_line(self):
        """PDB line for bonding within this molecule.

        Returns:
            String with PDB line.
        """
        res = 'CONECT%5d' % self.numb

        bonded = []
        for atom in self.bonded_atoms:
            bonded.append(atom.numb)
        bonded.sort()

        for b in bonded:
            res += '%5d'%b
        res += '\n'
        return res

    def get_input_parameters(self):
        """Extract the input parameters stored in the occupancy and b-factor
        fields in input files"""
        # Set the group type
        if self.occ != '-':
            # make sure to set the terminal
            if self.occ in ['N+', 'C-']:
                self.terminal = self.occ
            # save the ligand group charge
            if self.occ == 'BLG':
                self.charge = +1
            elif self.occ == 'ALG':
                self.charge = -1
            # generic ions
            if self.occ in ['1P', '2P', '1N', '2N']:
                self.res_name = self.occ
                self.occ = 'Ion'
            # correct the group type
            self.occ = self.occ.replace('N+', 'Nterm')
            self.occ = self.occ.replace('C-', 'Cterm')
            self.occ = self.occ.replace('ION', 'Ion')
            self.occ = self.occ.replace('ALG', 'titratable_ligand')
            self.occ = self.occ.replace('BLG', 'titratable_ligand')
            self.occ = self.occ.replace('LG', 'non_titratable_ligand')
            # try to initialise the group
            try:
                # TODO - get rid of this exec() statement for security reasons
                exec('self.group = propka.group.%s_group(self)' % self.occ)
            except:
                raise Exception('%s in input_file is not recognized as a group' % self.occ)
        # set the model pKa value
        if self.beta != '-':
            self.group.model_pka = float(self.beta)
            self.group.model_pka_set = True
        # set occ and beta to standard values
        self.occ = '1.00'
        self.beta = '0.00'

    def make_pdb_line(self):
        """Create PDB line.

        TODO - this could/should be a @property method/attribute
        TODO - figure out difference between make_pdb_line, make_input_line, and make_pdb_line2

        Returns:
            String with PDB line.
        """
        str_ = "%-6s%5d " % (self.type.upper(), self.numb)
        str_ += "%s %s" % (propka.lib.makeTidyAtomLabel(self.name, self.element),
                           self.res_name)
        str_ += "%2s%4d%12.3lf%8.3lf%8.3lf%6s%6s\n" % (self.chain_id, self.res_num,
                                                       self.x, self.y, self.z,
                                                       self.occ, self.beta)
        return str_

    def make_mol2_line(self, id_):
        """Create MOL2 line.

        Format: 1      S1     3.6147     2.0531     1.4795     S.3     1       noname  -0.1785

        TODO - this could/should be a @property method/attribute

        Returns:
            String with MOL2 line.
        """
        str_ = "%-4d %-4s " % (id_, propka.lib.makeTidyAtomLabel(self.name,
                                                                 self.element))
        str_ += "%10.4f %10.4f %10.4f " % (self.x, self.y, self.z)
        str_ += "%6s %6d %10s %10.4f\n" % (self.sybyl_type.replace('-', ''),
                                           self.res_num, self.res_name, 0.0)
        return str_

    def make_pdb_line2(self, numb=None, name=None, res_name=None, chain_id=None,
                       res_num=None, x=None, y=None, z=None, occ=None, beta=None):
        """Create a PDB line.

        TODO - this could/should be a @property method/attribute
        TODO - figure out difference between make_pdb_line, make_input_line, and make_pdb_line2

        Returns:
            String with PDB line.
        """
        if numb is None:
            numb = self.numb
        if name is None:
            name = self.name
        if res_name is None:
            res_name = self.res_name
        if chain_id is None:
            chain_id = self.chain_id
        if res_num is None:
            res_num = self.res_num
        if x is None:
            x = self.x
        if y is None:
            y = self.y
        if z is None:
            z = self.z
        if occ is None:
            occ = self.occ
        if beta is None:
            beta = self.beta
        str_ = "ATOM "
        str_ += "%6d" % (numb)
        str_ += " %s" % (propka.lib.makeTidyAtomLabel(name, self.element))
        str_ += " %s" % (res_name)
        str_ += "%2s" % (chain_id)
        str_ += "%4d" % (res_num)
        str_ += "%12.3lf" % (x)
        str_ += "%8.3lf" % (y)
        str_ += "%8.3lf" % (z)
        str_ += "%6.2lf" % (occ)
        str_ += "%6.2lf" % (beta)
        str_ += '\n'
        return str_

    def get_tidy_label(self):
        """Returns a 'tidier' atom label for printing the new pdbfile

        TODO - this could/should be a @property method/attribute

        Returns:
            String with label"""
        return propka.lib.makeTidyAtomLabel(self.name, self.element)

    def __str__(self):
        """Return an undefined-format string version of this atom."""
        return '%5d-%4s %5d-%3s (%1s) [%8.3f %8.3f %8.3f] %s' % (self.numb, self.name,
                                                                 self.res_num, self.res_name,
                                                                 self.chain_id, self.x, self.y,
                                                                 self.z, self.element)

    def set_residue(self, residue):
        """ Makes a reference to the parent residue

        Args:
            residue:  the parent residue
        """
        if self.residue is None:
            self.residue = residue
