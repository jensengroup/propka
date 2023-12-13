"""
Ligand atom typing
==================

This module contains the :func:`assign_sybyl_type` function to analyze
all :class:`propka.atom.Atom` in terms of SYBYL atom types (see
:data:`ALL_SYBYL_TYPES`).

"""

from propka.calculations import squared_distance
from propka.vector_algebra import Vector

#: SYBYL atom types
ALL_SYBYL_TYPES = [
    'C.3',   #  carbon sp3
    'H',     #  hydrogen
    'C.2',   #  carbon sp2
    'H.spc', #  hydrogen in Single Point Charge (SPC) water model
    'C.1',   #  carbon sp
    'H.t3p', #  hydrogen in Transferable intermolecular Potential (TIP3P) water model
    'C.ar',  #  carbon aromatic
    'LP',    #  lone pair
    'C.cat', #  carbocation (C+) used only in a guadinium group
    'Du',    #  dummy atom
    'N.3',   #  nitrogen sp3
    'Du.C',  #  dummy carbon
    'N.2',   #  nitrogen sp2
    'Any',   #  any atom
    'N.1',   #  nitrogen sp
    'Hal',   #  halogen
    'N.ar',  #  nitrogen aromatic
    'Het',   #  heteroatom = N, O, S, P
    'N.am',  #  nitrogen amide
    'Hev',   #  heavy atom (non hydrogen)
    'N.pl3', #  nitrogen trigonal planar
    'Li',    #  lithium
    'N.4',   #  nitrogen sp3 positively charged
    'Na',    #  sodium
    'O.3',   #  oxygen sp3
    'Mg',    #  magnesium
    'O.2',   #  oxygen sp2
    'Al',    #  aluminum
    'O.co2', #  oxygen in carboxylate and phosphate groups
    'Si',    #  silicon
    'O.spc', #  oxygen in Single Point Charge (SPC) water model
    'K',     #  potassium
    'O.t3p', #  oxygen in Transferable Intermolecular Potential (TIP3P) water model
    'Ca',    #  calcium
    'S.3',   #  sulfur sp3
    'Cr.th', #  chromium (tetrahedral)
    'S.2',   #  sulfur sp2
    'Cr.oh', #  chromium (octahedral)
    'S.O',   #  sulfoxide sulfur
    'Mn',    #  manganese
    'S.O2',  #  sulfone sulfur
    'Fe',    #  iron
    'P.3',   #  phosphorous sp3
    'Co.oh', #  cobalt (octahedral)
    'F',     #  fluorine
    'Cu',    #  copper
    'Cl',    #  chlorine
    'Zn',    #  zinc
    'Br',    #  bromine
    'Se',    #  selenium
    'I',     #  iodine
    'Mo',    #  molybdenum
    'Sn']    #  tin

#: PROPKA input types
PROPKA_INPUT_TYPES = ['1P', '1N', '2P', '2N', 'C3', 'H', 'C2', 'Hsp', 'C1',
                      'Ht3', 'Car', 'LP', 'Cca', 'Du', 'N3', 'DuC', 'N2',
                      'Any', 'N1', 'Hal', 'Nar', 'Het', 'Nam', 'Hev', 'Npl',
                      'Li', 'N4', 'Na', 'O3', 'Mg', 'O2', 'Al', 'Oco', 'Si',
                      'Osp', 'K', 'Ot3', 'Ca', 'S3', 'Crt', 'S2', 'Cro', 'SO',
                      'Mn', 'SO2', 'Fe', 'P3', 'Coo', 'F', 'Cu', 'Cl', 'Zn',
                      'Br', 'Se', 'I', 'Mo', 'Sn']


MAX_C_DOUBLE_BOND = 1.3
MAX_C_TRIPLE_BOND = 1.2
MAX_C_DOUBLE_BOND_SQUARED = MAX_C_DOUBLE_BOND*MAX_C_DOUBLE_BOND
MAX_C_TRIPLE_BOND_SQUARED = MAX_C_TRIPLE_BOND*MAX_C_TRIPLE_BOND
PLANARITY_MARGIN = 0.20


def assign_sybyl_type(atom):
    """Assign Sybyl type to atom.

    Args:
        atom:  atom to assign
    """
    # check if we already have assigned a name to this atom
    if atom.sybyl_assigned:
        return
    # find some properties of the atom
    ring_atoms = is_ring_member(atom)
    aromatic = is_aromatic_ring(ring_atoms)
    planar = is_planar(atom)
    bonded_elements = {}
    for i, bonded_atom in enumerate(atom.bonded_atoms):
        bonded_elements[i] = bonded_atom.element
    # Aromatic carbon/nitrogen
    if aromatic:
        for ring_atom in ring_atoms:
            if ring_atom.element in ['C', 'N']:
                set_type(ring_atom, ring_atom.element+'.ar')
        return
    # check for amide
    if atom.element in ['O', 'N', 'C']:
        o_atom = None
        n_atom = None
        c_atom = None
        # oxygen, nitrogen
        if atom.element in ['O', 'N']:
            for bonded_elem in atom.get_bonded_elements('C'):
                for bonded_atom in bonded_elem.bonded_atoms:
                    if (bonded_atom.element == 'N' and atom.element == 'O'):
                        o_atom = atom
                        c_atom = bonded_elem
                        n_atom = bonded_atom
                    elif (bonded_atom.element == 'O' and atom.element == 'N'):
                        n_atom = atom
                        c_atom = bonded_elem
                        o_atom = bonded_atom
        # carbon
        if atom.element == 'C':
            nitrogens = atom.get_bonded_elements('N')
            oxygens = atom.get_bonded_elements('O')
            if len(nitrogens) == 1 and len(oxygens) == 1:
                c_atom = atom
                n_atom = nitrogens[0]
                o_atom = oxygens[0]
        if c_atom and n_atom and o_atom:
            # make sure that the Nitrogen is not aromatic and that it has two
            # heavy atom bonds
            if (not is_aromatic_ring(is_ring_member(n_atom))
                    and len(n_atom.get_bonded_heavy_atoms()) == 2):
                set_type(n_atom, 'N.am')
                set_type(c_atom, 'C.2')
                set_type(o_atom, 'O.2')
                return
    if atom.element == 'C':
        # check for carboxyl
        if (len(atom.bonded_atoms) == 3
                and list(bonded_elements.values()).count('O') == 2):
            index1 = list(bonded_elements.values()).index('O')
            index2 = list(bonded_elements.values()).index('O', index1+1)
            if (len(atom.bonded_atoms[index1].bonded_atoms) == 1
                    and len(atom.bonded_atoms[index2].bonded_atoms) == 1):
                set_type(atom.bonded_atoms[index1], 'O.co2-')
                set_type(atom.bonded_atoms[index2], 'O.co2')
                set_type(atom, 'C.2')
                return
        # sp carbon
        if len(atom.bonded_atoms) <= 2:
            for bonded_atom in atom.bonded_atoms:
                if (squared_distance(atom, bonded_atom)
                        < MAX_C_TRIPLE_BOND_SQUARED):
                    set_type(atom, 'C.1')
                    set_type(bonded_atom, bonded_atom.element + '.1')
            if atom.sybyl_assigned:
                return
        # sp2 carbon
        if planar:
            set_type(atom, 'C.2')
            # check for N.pl3
            for bonded_atom in atom.bonded_atoms:
                if bonded_atom.element == 'N':
                    if (len(bonded_atom.bonded_atoms) < 3
                            or is_planar(bonded_atom)):
                        set_type(bonded_atom, 'N.pl3')
            return
        # sp3 carbon
        set_type(atom, 'C.3')
        return
    # Nitrogen
    if atom.element == 'N':
        # check for planar N
        if len(atom.bonded_atoms) == 1:
            if is_planar(atom.bonded_atoms[0]):
                set_type(atom, 'N.pl3')
                return
        if planar:
            set_type(atom, 'N.pl3')
            return
        set_type(atom, 'N.3')
        return
    # Oxygen
    if atom.element == 'O':
        set_type(atom, 'O.3')
        if len(atom.bonded_atoms) == 1:
            # check for carboxyl
            if atom.bonded_atoms[0].element == 'C':
                the_carbon = atom.bonded_atoms[0]
                if (len(the_carbon.bonded_atoms) == 3
                        and the_carbon.count_bonded_elements('O') == 2):
                    [oxy1, oxy2] = the_carbon.get_bonded_elements('O')
                    if (len(oxy1.bonded_atoms) == 1
                            and len(oxy2.bonded_atoms) == 1):
                        set_type(oxy1, 'O.co2-')
                        set_type(oxy2, 'O.co2')
                        set_type(the_carbon, 'C.2')
                        return
            # check for X=O
            if (squared_distance(atom, atom.bonded_atoms[0])
                    < MAX_C_DOUBLE_BOND_SQUARED):
                set_type(atom, 'O.2')
                if atom.bonded_atoms[0].element == 'C':
                    set_type(atom.bonded_atoms[0], 'C.2')
        return
    # Sulphur
    if atom.element == 'S':
        # check for SO2
        if list(bonded_elements.values()).count('O') == 2:
            index1 = list(bonded_elements.values()).index('O')
            index2 = list(bonded_elements.values()).index('O', index1+1)
            set_type(atom.bonded_atoms[index1], 'O.2')
            set_type(atom.bonded_atoms[index2], 'O.2')
            set_type(atom, 'S.o2')
            return
        # check for SO4
        if list(bonded_elements.values()).count('O') == 4:
            no_o2 = 0
            for i in range(len(atom.bonded_atoms)):
                if len(atom.bonded_atoms[i].bonded_atoms) == 1 and no_o2 < 2:
                    set_type(atom.bonded_atoms[i], 'O.2')
                    no_o2 += 1
                else:
                    set_type(atom.bonded_atoms[i], 'O.3')
        set_type(atom, 'S.3')
        return
    # Phosphorus
    if atom.element == 'P':
        set_type(atom, 'P.3')
        # check for phosphate group
        bonded_oxygens = atom.get_bonded_elements('O')
        for o_atom in bonded_oxygens:
            set_type(o_atom, 'O.3')
        return
    element = atom.element.capitalize()
    set_type(atom, element)


def is_ring_member(atom):
    """Determine if atom is a member of a ring.

    Args:
        atom:  atom to test
    Returns:
        list of atoms
    """
    return identify_ring(atom, atom, 0, [])


def identify_ring(this_atom, original_atom, number, past_atoms):
    """Identify the atoms in a ring

    Args:
        this_atom:  atom to test
        original_atom:  some other atom
        number:  number of atoms
        past_atoms:  atoms that have already been found
    Returns:
        list of atoms
    """
    number += 1
    past_atoms = past_atoms + [this_atom]
    return_atoms = []
    if number > 10:
        return return_atoms
    for atom in this_atom.get_bonded_heavy_atoms():
        if atom == original_atom and number > 2:
            return past_atoms
        if atom not in past_atoms:
            these_return_atoms = identify_ring(atom, original_atom, number,
                                               past_atoms)
            if len(these_return_atoms) > 0:
                if (len(return_atoms) > len(these_return_atoms)
                        or len(return_atoms) == 0):
                    return_atoms = these_return_atoms
    return return_atoms


def set_type(atom, type_):
    """Set atom type..

    Args:
        atom:  atom to set
        type_:  type value to set
    """
    atom.sybyl_type = type_
    atom.sybyl_assigned = True


def is_planar(atom):
    """Finds out if atom forms a plane together with its bonded atoms.

    Args:
        atom:  atom to test
    Returns:
        Boolean
    """
    atoms = [atom] + atom.bonded_atoms
    return are_atoms_planar(atoms)


def are_atoms_planar(atoms):
    """Test whether a group of atoms are planar.

    Args:
        atoms:  list of atoms
    Returns:
        Boolean
    """
    if len(atoms) == 0:
        return False
    if len(atoms) < 4:
        return False
    vec1 = Vector(atom1=atoms[0], atom2=atoms[1])
    vec2 = Vector(atom1=atoms[0], atom2=atoms[2])
    norm = vec1.cross(vec2).rescale(1.0)
    margin = PLANARITY_MARGIN
    for atom in atoms[3:]:
        vec = Vector(atom1=atoms[0], atom2=atom).rescale(1.0)
        if abs(vec.dot(norm)) > margin:
            return False
    return True


def is_aromatic_ring(atoms):
    """Determine whether group of atoms form aromatic ring.

    Args:
        atoms:  list of atoms to test
    Returns:
        Boolean
    """
    if len(atoms) < 5:
        return False
    for i in range(len(atoms)):
        if not are_atoms_planar(atoms[i:]+atoms[:i]):
            return False
    return True
