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
        self.resName = None
        self.type = None
        self.chainID = None
        self.beta = None
        self.icode = None
        self.resNumb = None
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
        self.steric_number_and_lone_pairs_set = False
        self.terminal = None
        self.charge = 0
        self.charge_set = False
        self.steric_number = 0
        self.number_of_lone_pairs = 0
        self.number_of_protons_to_add = 0
        self.number_of_pi_electrons_in_double_and_triple_bonds = 0
        self.number_of_pi_electrons_in_conjugate_double_and_triple_bonds = 0
        self.groups_extracted = 0
        self.set_properties(line)
        self.residue_label = "%-3s%4d%2s" % (self.name, self.resNumb, self.chainID)

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
        self.resNumb = 0
        self.resName = ''
        self.chainID = 'A'
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
            self.resNumb = int(line[22:26].strip())
            self.resName = "%-3s" % (line[17:20].strip())
            self.chainID = line[21]
            # Set chain id to "_" if it is just white space.
            if not self.chainID.strip():
                self.chainID = '_'
            self.type = line[:6].strip().lower()

            if self.resName in ['DA ', 'DC ', 'DG ', 'DT ']:
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

        TODO - this function is silly.  It should just be the len() of the
        array returned by get_bonded_elements()

        Args:
            element:  element type for test.
        Returns:
            number of bonded atoms.
        """
        res = 0
        for ba in self.bonded_atoms:
            if element == ba.element:
                res += 1
        return res

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

        TODO - this could be a @property attribute/method

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

    def setProperty(self, numb=None, name=None, resName=None, chainID=None,
                    resNumb=None, x=None, y=None, z=None, occ=None, beta=None):
        """Set properties of the atom object.

        Args:
            numb:  Atom number
            name:  Atom name
            resName:  residue name
            chainId:  chain ID
            resNumb:  residue number
            x:  atom x-coordinate
            y:  atom y-coordinate
            z:  atom z-coordinate
            occ:  atom occupancy
            beta:  atom temperature factor
        """
        if numb != None:
            self.numb = numb
        if name != None:
            self.name = name
        if resName != None:
            self.resName = resName
        if chainID != None:
            self.chainID = chainID
        if resNumb != None:
            self.resNumb = resNumb
        if x != None:
            self.x = x
        if y != None:
            self.y = y
        if z != None:
            self.z = z
        if occ != None:
            self.occ = occ
        if beta != None:
            self.beta = beta

    def makeCopy(self):
        """Make a copy of this atom.

        TODO - this could be a property method/attribute

        Returns:
            Another atom object copy of this one.
        """
        newAtom = Atom()
        newAtom.type = self.type
        newAtom.numb = self.numb
        newAtom.name = self.name
        newAtom.element = self.element
        newAtom.resName = self.resName
        newAtom.resNumb = self.resNumb
        newAtom.chainID = self.chainID
        newAtom.x = self.x
        newAtom.y = self.y
        newAtom.z = self.z
        newAtom.occ = self.occ
        newAtom.beta = self.beta
        newAtom.terminal = self.terminal
        newAtom.residue_label = self.residue_label
        newAtom.icode = self.icode
        return newAtom

    def make_input_line(self):
        """PDB line for this atom.

        TODO - Could be @property method/attribute
        TODO - figure out difference between make_pdb_line, make_input_line, and makePDBLine

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
        str_ += "%s%2s%4d%12.3lf%8.3lf%8.3lf%6s%6s \n" % (self.resName, self.chainID,
                                                          self.resNumb, self.x, self.y,
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
                self.resName = self.occ
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
        TODO - figure out difference between make_pdb_line, make_input_line, and makePDBLine

        Returns:
            String with PDB line.
        """
        str_ = "%-6s%5d " % (self.type.upper(), self.numb)
        str_ += "%s %s" % (propka.lib.makeTidyAtomLabel(self.name, self.element),
                           self.resName)
        str_ += "%2s%4d%12.3lf%8.3lf%8.3lf%6s%6s\n" % (self.chainID, self.resNumb,
                                                       self.x, self.y, self.z,
                                                       self.occ, self.beta)
        return str_

    def make_mol2_line(self, id_):
        """Create MOL2 line.

        Format: 1      S1     3.6147     2.0531     1.4795     S.3     1       noname  -0.1785

        TODO - this could/shoudl be a @property method/attribute

        Returns:
            String with MOL2 line.
        """
        str_ = "%-4d %-4s " % (id_, propka.lib.makeTidyAtomLabel(self.name,
                                                                 self.element))
        str_ += "%10.4f %10.4f %10.4f " % (self.x, self.y, self.z)
        str_ += "%6s %6d %10s %10.4f\n" % (self.sybyl_type.replace('-', ''),
                                           self.resNumb, self.resName, 0.0)
        return str_

    def makePDBLine(self, numb=None, name=None, resName=None, chainID=None,
                    resNumb=None, x=None, y=None, z=None, occ=None, beta=None):
        """Create a PDB line.

        TODO - this could/should be a @property method/attribute
        TODO - figure out difference between make_pdb_line, make_input_line, and makePDBLine

        Returns:
            String with PDB line.
        """
        if numb is None:
            numb = self.numb
        if name is None:
            name = self.name
        if resName is None:
            resName = self.resName
        if chainID is None:
            chainID = self.chainID
        if resNumb is None:
            resNumb = self.resNumb
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
        str_ += " %s" % (resName)
        str_ += "%2s" % (chainID)
        str_ += "%4d" % (resNumb)
        str_ += "%12.3lf" % (x)
        str_ += "%8.3lf" % (y)
        str_ += "%8.3lf" % (z)
        str_ += "%6.2lf" % (occ)
        str_ += "%6.2lf" % (beta)
        str_ += '\n'
        return str_

    def getTidyLabel(self):
        """Returns a 'tidier' atom label for printing the new pdbfile

        TODO - this could/should be a @property method/attribute

        Returns:
            String with label"""
        return propka.lib.makeTidyAtomLabel(self.name, self.element)

    def __str__(self):
        """Return an undefined-format string version of this atom."""
        return '%5d-%4s %5d-%3s (%1s) [%8.3f %8.3f %8.3f] %s' % (self.numb, self.name,
                                                                 self.resNumb, self.resName,
                                                                 self.chainID, self.x, self.y,
                                                                 self.z, self.element)

    def set_residue(self, residue):
        """ Makes a reference to the parent residue

        Args:
            residue:  the parent residue
        """
        if self.residue is None:
            self.residue = residue
