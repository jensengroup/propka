"""
Atom
====

The :class:`Atom` class contains all atom information found in the PDB file.

"""

import string
from typing import cast, List, NoReturn, Optional, TYPE_CHECKING
import warnings

from propka.lib import make_tidy_atom_label
from . import hybrid36

if TYPE_CHECKING:
    from propka.group import Group
    from propka.molecular_container import MolecularContainer
    from propka.conformation_container import ConformationContainer

# Format strings that get used in multiple places (or are very complex)
PDB_LINE_FMT1 = (
    "{type:6s}{r.numb:>5d} {atom_label} {r.res_name}{r.chain_id:>2s}"
    "{r.res_num:>4d}{r.x:>12.3f}{r.y:>8.3f}{r.z:>8.3f}{r.occ:>6s}"
    "{r.beta:>6s}\n")
MOL2_LINE_FMT = (
    "{id:<4d} {atom_label:4s} "
    "{r.x:>10.4f} {r.y:>10.4f} {r.z:>10.4f} "
    "{r.sybyl_type:>6s} {r.res_num:>6d} {r.res_name:>10s}     0.0000\n")
PDB_LINE_FMT2 = (
    "ATOM {numb:>6d} {atom_label} {res_name}{chain_id:>2s}{res_num:>4d}"
    "{x:>12.3f}{y:>8.3f}{z:>8.3f}{occ:>6.2f}{beta:>6.2f}\n")
STR_FMT = (
    "{r.numb:>5d}-{r.name:>4s} {r.res_num:>5d}-{r.res_name:>3s} "
    "({r.chain_id:1s}) [{r.x:>8.3f} {r.y:>8.3f} {r.z:>8.3f}] {r.element:s}")


class Atom:
    """Atom class - contains all atom information found in the PDB file


    .. versionchanged:: 3.4.0
       :meth:`make_input_line` and :meth:`get_input_parameters` have been
       removed as reading/writing PROPKA input is no longer supported.
    """
    group: Optional["Group"] = None
    group_type: Optional[str] = None
    cysteine_bridge: bool = False
    residue: NoReturn = None  # type: ignore[assignment]
    conformation_container: Optional["ConformationContainer"] = None
    molecular_container: Optional["MolecularContainer"] = None
    is_protonated: bool = False
    steric_num_lone_pairs_set: bool = False
    terminal: Optional[str] = None
    charge: float = 0.0
    charge_set: bool = False
    steric_number: int = 0
    number_of_lone_pairs: int = 0
    number_of_protons_to_add: int = 0
    num_pi_elec_2_3_bonds: int = 0
    num_pi_elec_conj_2_3_bonds: int = 0
    groups_extracted: bool = False

    # PDB attributes
    name: str = ''
    numb: int = 0
    x: float = 0.0
    y: float = 0.0
    z: float = 0.0
    res_num: int = 0
    res_name: str = ''
    chain_id: str = 'A'
    type: str = ''
    occ: str = '1.0'
    beta: str = '0.0'
    element: str = ''
    icode: str = ''

    # ligand atom types
    sybyl_type = ''
    sybyl_assigned = False
    marvin_pka = False

    def __init__(self, line: Optional[str] = None):
        """Initialize Atom object.

        Args:
            line:  Line from a PDB file to set properties of atom.
        """
        self.number_of_bonded_elements: NoReturn = cast(NoReturn, {})  # FIXME unused?
        self.bonded_atoms: List[Atom] = []
        self.set_properties(line)
        fmt = "{r.name:3s}{r.res_num:>4d}{r.chain_id:>2s}"
        self.residue_label = fmt.format(r=self)

    def set_properties(self, line: Optional[str]):
        """Line from PDB file to set properties of atom.

        Args:
            line:  PDB file line
        """
        if line:
            self.name = line[12:16].strip()
            self.numb = int(hybrid36.decode(line[6:11]))
            self.x = float(line[30:38].strip())
            self.y = float(line[38:46].strip())
            self.z = float(line[46:54].strip())
            self.res_num = int(line[22:26].strip())
            self.res_name = "{0:<3s}".format(line[17:20].strip())
            # Set chain id to "_" if it is just white space.
            self.chain_id = line[21].strip() or '_'
            self.type = line[:6].strip().lower()

            # TODO - define nucleic acid residue names elsewhere
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
                self.element = '{0:1s}{1:1s}'.format(
                    self.element[0], self.element[1].lower())

    def set_group_type(self, type_: str):
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
        for bond_atom in self.bonded_atoms:
            if bond_atom.element == element:
                res.append(bond_atom)
        return res

    def get_bonded_heavy_atoms(self):
        """Get the atoms bonded to this one that aren't hydrogen.

        Returns:
            list of atoms.
        """
        return [ba for ba in self.bonded_atoms if ba.element != 'H']

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
                if ba.is_atom_within_bond_distance(other_atom, max_bonds,
                                                   cur_bond+1):
                    return True
        return False

    def set_property(self,
                     numb: Optional[int] = None,
                     name: Optional[str] = None,
                     res_name: Optional[str] = None,
                     chain_id: Optional[str] = None,
                     res_num: Optional[int] = None,
                     x: Optional[float] = None,
                     y: Optional[float] = None,
                     z: Optional[float] = None,
                     occ: Optional[str] = None,
                     beta: Optional[str] = None):
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

    def make_conect_line(self):
        """PDB line for bonding within this molecule.

        Returns:
            String with PDB line.
        """
        res = 'CONECT{0:5d}'.format(self.numb)

        bonded = []
        for atom in self.bonded_atoms:
            bonded.append(atom.numb)
        bonded.sort()

        for bond in bonded:
            res += '{0:5d}'.format(bond)
        res += '\n'
        return res

    def make_pdb_line(self):
        """Create PDB line.

        TODO - this could/should be a @property method/attribute
        TODO - figure out difference between make_pdb_line, and make_pdb_line2

        Returns:
            String with PDB line.
        """
        str_ = PDB_LINE_FMT1.format(
            type=self.type.upper(), r=self,
            atom_label=make_tidy_atom_label(self.name, self.element))
        return str_

    def make_mol2_line(self, id_):
        """Create MOL2 line.

        Format: 1      S1     3.6147     2.0531     1.4795     S.3     1       noname  -0.1785

        TODO - this could/should be a @property method/attribute

        Returns:
            String with MOL2 line.
        """
        str_ = MOL2_LINE_FMT.format(
            id=id_, r=self,
            atom_label=make_tidy_atom_label(self.name, self.element))
        return str_

    def make_pdb_line2(self, numb=None, name=None, res_name=None, chain_id=None,
                       res_num=None, x=None, y=None, z=None, occ=None,
                       beta=None):
        """Create a PDB line.

        TODO - this could/should be a @property method/attribute
        TODO - figure out difference between make_pdb_line, and make_pdb_line2

        Returns:
            String with PDB line.
        """
        warnings.warn("only used by unused function")
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
        str_ = PDB_LINE_FMT2.format(
            numb=numb, res_name=res_name, chain_id=chain_id, res_num=res_num,
            x=x, y=y, z=z, occ=occ, beta=beta,
            atom_label=make_tidy_atom_label(name, self.element)
        )
        return str_

    def get_tidy_label(self):
        """Returns a 'tidier' atom label for printing the new pdbfile

        TODO - this could/should be a @property method/attribute

        Returns:
            String with label"""
        return make_tidy_atom_label(self.name, self.element)

    def __str__(self):
        """Return an undefined-format string version of this atom."""
        return STR_FMT.format(r=self)

    def set_residue(self, residue: NoReturn):
        """ Makes a reference to the parent residue

        Args:
            residue:  the parent residue
        """
        raise NotImplementedError("unused")
        if self.residue is None:
            self.residue = residue
