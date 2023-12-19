"""
Hydrogens
=========

Calculations related to hydrogen placement.

"""
import math
import logging
from typing import List, Optional, Tuple, TYPE_CHECKING

from propka.protonate import Protonate
from propka.bonds import BondMaker
from propka.atom import Atom

if TYPE_CHECKING:
    from propka.molecular_container import MolecularContainer

_LOGGER = logging.getLogger(__name__)


def setup_bonding_and_protonation(molecular_container: "MolecularContainer") -> None:
    """Set up bonding and protonation for a molecule.

    Args:
        parameters:  not used
        molecular_container:  molecule container.
    """
    # make bonds
    my_bond_maker = setup_bonding(molecular_container)
    # set up ligand atom names
    set_ligand_atom_names(molecular_container)
    # apply information on pi electrons
    my_bond_maker.add_pi_electron_information(molecular_container)
    # Protonate atoms
    if molecular_container.options.protonate_all:
        protonator = Protonate(verbose=False)
        protonator.protonate(molecular_container)


def setup_bonding(molecular_container: "MolecularContainer") -> BondMaker:
    """Set up bonding for a molecular container.

    Args:
        molecular_container:  the molecular container in question
    Returns:
        BondMaker object
    """
    my_bond_maker = BondMaker()
    my_bond_maker.find_bonds_for_molecules_using_boxes(molecular_container)
    return my_bond_maker


def setup_bonding_and_protonation_30_style(molecular_container: "MolecularContainer") -> BondMaker:
    """Set up bonding for a molecular container.

    Args:
        molecular_container:  the molecular container in question
    Returns:
        BondMaker object
    """
    # Protonate atoms
    protonate_30_style(molecular_container)
    # make bonds
    bond_maker = BondMaker()
    bond_maker.find_bonds_for_molecules_using_boxes(molecular_container)
    return bond_maker


def protonate_30_style(molecular_container: "MolecularContainer") -> None:
    """Protonate the molecule.

    Args:
        molecular_container:  molecule
    """
    for name in molecular_container.conformation_names:
        _LOGGER.info('Now protonating %s', name)
        # split atom into residues
        curres = -1000000
        residue: List[Atom] = []
        o_atom: Optional[Atom] = None
        c_atom: Optional[Atom] = None
        for atom in molecular_container.conformations[name].atoms:
            if atom.res_num != curres:
                curres = atom.res_num
                if len(residue) > 0:
                    # backbone
                    [o_atom, c_atom] = add_backbone_hydrogen(
                        residue, o_atom, c_atom)
                    # arginine
                    if residue[0].res_name == 'ARG':
                        add_arg_hydrogen(residue)
                    # histidine
                    if residue[0].res_name == 'HIS':
                        add_his_hydrogen(residue)
                    # tryptophan
                    if residue[0].res_name == 'TRP':
                        add_trp_hydrogen(residue)
                    # amides
                    if residue[0].res_name in ['GLN', 'ASN']:
                        add_amd_hydrogen(residue)
                    residue = []
            if atom.type == 'atom':
                residue.append(atom)


def set_ligand_atom_names(molecular_container: "MolecularContainer") -> None:
    """Set names for ligands in molecular container.

    Args:
        molecular_container:  molecular container for ligand names
    """
    for name in molecular_container.conformation_names:
        molecular_container.conformations[name].set_ligand_atom_names()


def add_arg_hydrogen(residue: List[Atom]) -> List[Atom]:
    """Adds Arg hydrogen atoms to residues according to the 'old way'.

    Args:
        residue:  arginine residue to protonate
    Returns:
        list of hydrogen atoms
    """
    cd_atom: Optional[Atom] = None
    cz_atom: Optional[Atom] = None
    ne_atom: Optional[Atom] = None
    nh1_atom: Optional[Atom] = None
    nh2_atom: Optional[Atom] = None

    for atom in residue:
        if atom.name == "CD":
            cd_atom = atom
        elif atom.name == "CZ":
            cz_atom = atom
        elif atom.name == "NE":
            ne_atom = atom
        elif atom.name == "NH1":
            nh1_atom = atom
        elif atom.name == "NH2":
            nh2_atom = atom

    if (cd_atom is None or cz_atom is None or ne_atom is None or nh1_atom is None
            or nh2_atom is None):
        raise ValueError("Unable to find all atoms")

    h1_atom = protonate_sp2(cd_atom, ne_atom, cz_atom)
    h1_atom.name = "HE"
    h2_atom = protonate_direction(nh1_atom, ne_atom, cz_atom)
    h2_atom.name = "HN1"
    h3_atom = protonate_direction(nh1_atom, ne_atom, cd_atom)
    h3_atom.name = "HN2"
    h4_atom = protonate_direction(nh2_atom, ne_atom, cz_atom)
    h4_atom.name = "HN3"
    h5_atom = protonate_direction(nh2_atom, ne_atom, h1_atom)
    h5_atom.name = "HN4"
    return [h1_atom, h2_atom, h3_atom, h4_atom, h5_atom]


def add_his_hydrogen(residue: List[Atom]) -> None:
    """Adds His hydrogen atoms to residues according to the 'old way'.

    Args:
        residue:  histidine residue to protonate
    """
    cg_atom: Optional[Atom] = None
    nd_atom: Optional[Atom] = None
    cd_atom: Optional[Atom] = None
    ce_atom: Optional[Atom] = None
    ne_atom: Optional[Atom] = None

    for atom in residue:
        if atom.name == "CG":
            cg_atom = atom
        elif atom.name == "ND1":
            nd_atom = atom
        elif atom.name == "CD2":
            cd_atom = atom
        elif atom.name == "CE1":
            ce_atom = atom
        elif atom.name == "NE2":
            ne_atom = atom

    if (cg_atom is None or nd_atom is None or cd_atom is None or ce_atom is None
            or ne_atom is None):
        raise ValueError("Unable to find all atoms")

    hd_atom = protonate_sp2(cg_atom, nd_atom, ce_atom)
    hd_atom.name = "HND"
    he_atom = protonate_sp2(cd_atom, ne_atom, ce_atom)
    he_atom.name = "HNE"


def add_trp_hydrogen(residue: List[Atom]) -> None:
    """Adds Trp hydrogen atoms to residues according to the 'old way'.

    Args:
        residue:  tryptophan residue to protonate
    """
    cd_atom = None
    ne_atom = None
    ce_atom = None
    for atom in residue:
        if atom.name == "CD1":
            cd_atom = atom
        elif atom.name == "NE1":
            ne_atom = atom
        elif atom.name == "CE2":
            ce_atom = atom
    if (cd_atom is None) or (ne_atom is None) or (ce_atom is None):
        str_ = "Unable to find all atoms for {0:s} {1:s}".format(
            residue[0].res_name, residue[0].res_num)
        raise ValueError(str_)
    he_atom = protonate_sp2(cd_atom, ne_atom, ce_atom)
    he_atom.name = "HNE"


def add_amd_hydrogen(residue: List[Atom]) -> None:
    """Adds Gln & Asn hydrogen atoms to residues according to the 'old way'.

    Args:
        residue:  glutamine or asparagine residue to protonate
    """
    c_atom = None
    o_atom = None
    n_atom = None
    for atom in residue:
        if ((atom.res_name == "GLN" and atom.name == "CD")
                or (atom.res_name == "ASN" and atom.name == "CG")):
            c_atom = atom
        elif ((atom.res_name == "GLN" and atom.name == "OE1")
              or (atom.res_name == "ASN" and atom.name == "OD1")):
            o_atom = atom
        elif ((atom.res_name == "GLN" and atom.name == "NE2")
              or (atom.res_name == "ASN" and atom.name == "ND2")):
            n_atom = atom
    if (c_atom is None) or (o_atom is None) or (n_atom is None):
        str_ = "Unable to find all atoms for {0:s} {1:s}".format(
            residue[0].res_name, residue[0].res_num)
        raise ValueError(str_)
    h1_atom = protonate_direction(n_atom, o_atom, c_atom)
    h1_atom.name = "HN1"
    h2_atom = protonate_average_direction(n_atom, c_atom, o_atom)
    h2_atom.name = "HN2"


def add_backbone_hydrogen(residue: List[Atom],
                          o_atom: Optional[Atom],
                          c_atom: Optional[Atom]) -> Tuple[Optional[Atom], Optional[Atom]]:
    """Adds hydrogen backbone atoms to residues according to the old way.

    dR is wrong for the N-terminus (i.e. first residue) but it doesn't affect
    anything at the moment. Could be improved, but works for now.

    Args:
        residue:  residue to protonate
        o_atom:  backbone oxygen atom
        c_atom:  backbone carbon atom
    Returns:
        [new backbone oxygen atom, new backbone carbon atom]
    """
    new_c_atom = None
    new_o_atom = None
    n_atom = None
    for atom in residue:
        if atom.name == "N":
            n_atom = atom
        if atom.name == "C":
            new_c_atom = atom
        if atom.name == "O":
            new_o_atom = atom
    if c_atom is None or o_atom is None or n_atom is None:
        return (new_o_atom, new_c_atom)
    if n_atom.res_name == "PRO":
        # PRO doesn't have an H-atom; do nothing
        pass
    else:
        h_atom = protonate_direction(n_atom, o_atom, c_atom)
        h_atom.name = "H"
    return (new_o_atom, new_c_atom)


def protonate_direction(x1_atom: Atom, x2_atom: Atom, x3_atom: Atom) -> Atom:
    """Protonates an atom, x1_atom, given a direction.

    New direction for x1_atom proton is (x2_atom -> x3_atom).

    Args:
        x1_atom:  atom to be protonated
        x2_atom:  atom for direction
        x3_atom:  other atom for direction
    Returns:
        new hydrogen atom
    """
    dx = (x3_atom.x - x2_atom.x)
    dy = (x3_atom.y - x2_atom.y)
    dz = (x3_atom.z - x2_atom.z)
    length = math.sqrt(dx*dx + dy*dy + dz*dz)
    x = x1_atom.x + dx/length
    y = x1_atom.y + dy/length
    z = x1_atom.z + dz/length
    h_atom = make_new_h(x1_atom, x, y, z)
    h_atom.name = "H"
    return h_atom


def protonate_average_direction(x1_atom: Atom, x2_atom: Atom, x3_atom: Atom) -> Atom:
    """Protonates an atom, x1_atom, given a direction.

    New direction for x1_atom is (x1_atom/x2_atom -> x3_atom).
    Note, this one uses the average of x1_atom & x2_atom (N & O) unlike
    the previous N - C = O

    Args:
        x1_atom:  atom to be protonated
        x2_atom:  atom for direction
        x3_atom:  other atom for direction
    Returns:
        new hydrogen atom
    """
    dx = (x3_atom.x + x1_atom.x)*0.5 - x2_atom.x
    dy = (x3_atom.y + x1_atom.y)*0.5 - x2_atom.y
    dz = (x3_atom.z + x1_atom.z)*0.5 - x2_atom.z
    length = math.sqrt(dx*dx + dy*dy + dz*dz)
    x = x1_atom.x + dx/length
    y = x1_atom.y + dy/length
    z = x1_atom.z + dz/length
    h_atom = make_new_h(x1_atom, x, y, z)
    h_atom.name = "H"
    return h_atom


def protonate_sp2(x1_atom: Atom, x2_atom: Atom, x3_atom: Atom) -> Atom:
    """Protonates a SP2 atom, given a list of atoms

    Args:
        x1_atom:  atom to set direction
        x2_atom:  atom to be protonated
        x3_atom:  other atom to set direction
    Returns:
        new hydrogen atom
    """
    dx = (x1_atom.x + x3_atom.x)*0.5 - x2_atom.x
    dy = (x1_atom.y + x3_atom.y)*0.5 - x2_atom.y
    dz = (x1_atom.z + x3_atom.z)*0.5 - x2_atom.z
    length = math.sqrt(dx*dx + dy*dy + dz*dz)
    x = x2_atom.x - dx/length
    y = x2_atom.y - dy/length
    z = x2_atom.z - dz/length
    h_atom = make_new_h(x2_atom, x, y, z)
    h_atom.name = "H"
    return h_atom


def make_new_h(atom: Atom, x: float, y: float, z: float) -> Atom:
    """Add a new hydrogen to an atom at the specified position.

    Args:
        atom:  atom to protonate
        x:  x position of hydrogen
        y:  y position of hydrogen
        z:  z position of hydrogen
    Returns:
        new hydrogen atom
    """
    new_h = Atom()
    new_h.set_property(
        numb=None, name='H{0:s}'.format(atom.name[1:]),
        res_name=atom.res_name, chain_id=atom.chain_id,
        res_num=atom.res_num, x=x, y=y, z=z, occ=None, beta=None)
    new_h.element = 'H'
    new_h.bonded_atoms = [atom]
    new_h.charge = 0
    new_h.steric_number = 0
    new_h.number_of_lone_pairs = 0
    new_h.number_of_protons_to_add = 0
    new_h.num_pi_elec_2_3_bonds = 0
    atom.bonded_atoms.append(new_h)
    assert atom.conformation_container is not None
    atom.conformation_container.add_atom(new_h)
    return new_h
