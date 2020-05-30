"""PDB parsing functionality."""
import propka.lib
from propka.lib import warning
from propka.atom import Atom
from propka.conformation_container import ConformationContainer


EXPECTED_ATOM_NUMBERS = {'ALA': 5, 'ARG': 11, 'ASN': 8, 'ASP': 8, 'CYS': 6,
                         'GLY': 4, 'GLN': 9, 'GLU': 9, 'HIS': 10, 'ILE': 8,
                         'LEU': 8, 'LYS': 9, 'MET': 8, 'PHE': 11, 'PRO': 7,
                         'SER': 6, 'THR': 7, 'TRP': 14, 'TYR': 12, 'VAL': 7}


def read_pdb(pdb_file, parameters, molecule):
    """Parse a PDB file.

    Args:
        pdb_file:  file to read
        parameters:  parameters to guide parsing
        molecule:  molecular container
    Returns:
        list with elements:
            1. list of conformations
            2. list of names
    """
    conformations = {}
    # read in all atoms in the file
    lines = get_atom_lines_from_pdb(
        pdb_file, ignore_residues=parameters.ignore_residues,
        keep_protons=molecule.options.keep_protons,
        chains=molecule.options.chains)
    for (name, atom) in lines:
        if not name in conformations.keys():
            conformations[name] = ConformationContainer(
                name=name, parameters=parameters, molecular_container=molecule)
        conformations[name].add_atom(atom)
    # make a sorted list of conformation names
    names = sorted(conformations.keys(), key=propka.lib.conformation_sorter)
    return [conformations, names]


def protein_precheck(conformations, names):
    """Check protein for correct number of atoms, etc.

    Args:
        names:  conformation names to check
    """
    for name in names:
        atoms = conformations[name].atoms
        # Group the atoms by their residue:
        atoms_by_residue = {}
        for atom in atoms:
            if atom.element != 'H':
                res_id = resid_from_atom(atom)
                try:
                    atoms_by_residue[res_id].append(atom)
                except KeyError:
                    atoms_by_residue[res_id] = [atom]
        for res_id, res_atoms in atoms_by_residue.items():
            res_name = res_atoms[0].res_name
            residue_label = '{0:>3s}{1:>5s}'.format(res_name, res_id)
            # ignore ligand residues
            if res_name not in EXPECTED_ATOM_NUMBERS:
                continue
            # check for c-terminal
            if 'C-' in [a.terminal for a in res_atoms]:
                if len(res_atoms) != EXPECTED_ATOM_NUMBERS[res_name]+1:
                    str_ = ("Unexpected number ({num:d}) of atoms in residue "
                            "{res:s} in conformation {conf:s}".format(
                                num=len(res_atoms), res=residue_label,
                                conf=name))
                    warning(str_)
                continue
            # check number of atoms in residue
            if len(res_atoms) != EXPECTED_ATOM_NUMBERS[res_name]:
                str_ = ("Unexpected number ({num:d}) of atoms in residue "
                        "{res:s} in conformation {conf:s}".format(
                            num=len(res_atoms), res=residue_label,
                            conf=name))
                warning(str_)


def resid_from_atom(atom):
    """Return string with atom residue information.

    Args:
        atom:  atom to generate string for
    Returns
        string
    """
    return '{0:>4d} {1:s} {2:s}'.format(
        atom.res_num, atom.chain_id, atom.icode)


def get_atom_lines_from_pdb(pdb_file, ignore_residues=[], keep_protons=False,
                            tags=['ATOM  ', 'HETATM'], chains=None):
    """Get atom lines from PDB file.

    Args:
        pdb_file:  PDB file to parse
        ignore_residues:  list of residues to ignore
        keep_protons:  bool to keep/ignore protons
        tags:  tags of lines that include atoms
        chains:  list of chains
    """
    lines = propka.lib.open_file_for_reading(pdb_file).readlines()
    nterm_residue = 'next_residue'
    old_residue = None
    terminal = None
    model = 1
    for line in lines:
        tag = line[0:6]
        # set the model number
        if tag == 'MODEL ':
            model = int(line[6:])
            nterm_residue = 'next_residue'
        if tag == 'TER   ':
            nterm_residue = 'next_residue'
        if tag in tags:
            alt_conf_tag = line[16]
            residue_name = line[12: 16]
            residue_number = line[22: 26]
            # check if we want this residue
            if line[17: 20] in ignore_residues:
                continue
            if chains and line[21] not in chains:
                continue
            # set the Nterm residue number - nessecary because we may need to
            # identify more than one N+ group for structures with alt_conf tags
            if nterm_residue == 'next_residue' and tag == 'ATOM  ':
                # make sure that we reached a new residue - nessecary if OXT is
                # not the last atom inthe previous residue
                if old_residue != residue_number:
                    nterm_residue = residue_number
                    old_residue = None
            # Identify the configuration
            # convert digits to letters
            if alt_conf_tag in '123456789':
                alt_conf_tag = chr(ord(alt_conf_tag)+16)
            if alt_conf_tag == ' ':
                alt_conf_tag = 'A'
            conformation = '{0:d}{1:s}'.format(model, alt_conf_tag)
            # set the terminal
            if  tag == 'ATOM  ':
                if (residue_name.strip() == 'N'
                        and nterm_residue == residue_number):
                    terminal = 'N+'
                if  residue_name.strip() in ['OXT', 'O\'\'']:
                    terminal = 'C-'
                    nterm_residue = 'next_residue'
                    old_residue = residue_number
            # and yield the atom
            atom = Atom(line=line)
            atom.terminal = terminal
            #ignore hydrogen
            if not (atom.element == 'H' and not keep_protons):
                yield (conformation, atom)
            terminal = None


def write_pdb(conformation, filename):
    """Write PDB conformation to a file.

    Args:
        conformation:  conformation container
        filename:  filename for output
    """
    write_pdb_for_atoms(conformation.atoms, filename)


def write_pdb_for_atoms(atoms, filename, make_conect_section=False):
    """Write out PDB file for atoms.

    Args:
        atoms:  list of atoms
        filename:  name of file
        make_conect_section:  generate a CONECT PDB section
    """
    out = propka.lib.open_file_for_writing(filename)
    for atom in atoms:
        out.write(atom.make_pdb_line())
    if make_conect_section:
        for atom in atoms:
            out.write(atom.make_conect_line())
    out.close()


def write_mol2_for_atoms(atoms, filename):
    """Write out MOL2 file for atoms.

    Args:
        atoms:  list of atoms
        filename:  name of file
    """
    # TODO - header needs to be converted to format string
    header = '@<TRIPOS>MOLECULE\n\n{natom:d} {id:d}\nSMALL\nUSER_CHARGES\n'
    atoms_section = '@<TRIPOS>ATOM\n'
    for i, atom in enumerate(atoms):
        atoms_section += atom.make_mol2_line(i+1)
    bonds_section = '@<TRIPOS>BOND\n'
    id_ = 1
    for i, atom1 in enumerate(atoms):
        for j, atom2 in enumerate(atoms, i+1):
            if atom1 in atom2.bonded_atoms:
                type_ = get_bond_order(atom1, atom2)
                bonds_section += '{0:>7d} {1:>7d} {2:>7d} {3:>7s}\n'.format(
                    id_, i+1, j+1, type_)
                id_ += 1
    substructure_section = '@<TRIPOS>SUBSTRUCTURE\n\n'
    if len(atoms) > 0:
        substructure_section = (
            '@<TRIPOS>SUBSTRUCTURE\n{0:<7d} {1:>10s} {2:>7d}\n'.format(
                atoms[0].res_num, atoms[0].res_name, atoms[0].numb))
    out = propka.lib.open_file_for_writing(filename)
    out.write(header.format(natom=len(atoms), id=id_-1))
    out.write(atoms_section)
    out.write(bonds_section)
    out.write(substructure_section)
    out.close()


def get_bond_order(atom1, atom2):
    """Get the order of a bond between two atoms.

    Args:
        atom1:  first atom in bond
        atom2:  second atom in bond
    Returns:
        string with bond type
    """
    type_ = '1'
    pi_electrons1 = atom1.num_pi_elec_2_3_bonds
    pi_electrons2 = atom2.num_pi_elec_2_3_bonds
    if '.ar' in atom1.sybyl_type:
        pi_electrons1 -= 1
    if '.ar' in atom2.sybyl_type:
        pi_electrons2 -= 1
    if pi_electrons1 > 0 and pi_electrons2 > 0:
        type_ = '{0:d}'.format(min(pi_electrons1, pi_electrons2)+1)
    if '.ar' in atom1.sybyl_type and '.ar' in atom2.sybyl_type:
        type_ = 'ar'
    return type_


def write_input(molecular_container, filename):
    """Write PROPKA input file for molecular container.

    Args:
        molecular_container:  molecular container
        filename:  output file name
    """
    out = propka.lib.open_file_for_writing(filename)
    for conformation_name in molecular_container.conformation_names:
        out.write('MODEL {0:s}\n'.format(conformation_name))
        # write atoms
        for atom in molecular_container.conformations[conformation_name].atoms:
            out.write(atom.make_input_line())
        # write bonds
        for atom in molecular_container.conformations[conformation_name].atoms:
            out.write(atom.make_conect_line())
        # write covalently coupled groups
        for group in (
                molecular_container.conformations[conformation_name].groups):
            out.write(group.make_covalently_coupled_line())
        # write non-covalently coupled groups
        for group in (
                molecular_container.conformations[conformation_name].groups):
            out.write(group.make_non_covalently_coupled_line())
        out.write('ENDMDL\n')
    out.close()


def read_input(input_file, parameters, molecule):
    """Read PROPKA input file for molecular container.

    Args:
        input_file:  input file
        parameters:  parameters for parsing/setup
        molecule:  molecular container
    Returns:
        list with [conformations, names of conformations]
    """
    conformations = {}
    # read in all atoms in the input file
    lines = get_atom_lines_from_input(input_file)
    for (name, atom) in lines:
        if not name in conformations.keys():
            conformations[name] = ConformationContainer(
                name=name, parameters=parameters,
                molecular_container=molecule)
        conformations[name].add_atom(atom)
    # make a sorted list of conformation names
    names = sorted(conformations.keys(), key=propka.lib.conformation_sorter)
    return [conformations, names]


def get_atom_lines_from_input(input_file, tags=['ATOM  ', 'HETATM']):
    """Get atom lines from a PROPKA input file.

    Args:
        input_file:  input file
        tags:  tags defining atom lines
    Yields:
        conformation container, list of atoms
    """
    lines = propka.lib.open_file_for_reading(input_file).readlines()
    conformation = ''
    atoms = {}
    numbers = []
    for line in lines:
        tag = line[0:6]
        # set the conformation
        if tag == 'MODEL ':
            conformation = line[6:].strip()
        # found an atom - save it
        if tag in tags:
            atom = Atom(line=line)
            atom.get_input_parameters()
            atom.groups_extracted = 1
            atom.is_protonated = True
            atoms[atom.numb] = atom
            numbers.append(atom.numb)
        # found bonding information - apply it
        if tag == 'CONECT' and len(line) > 14:
            conect_numbers = [line[i:i+5] for i in range(6, len(line)-1, 5)]
            center_atom = atoms[int(conect_numbers[0])]
            for num in conect_numbers[1:]:
                bond_atom = atoms[int(num)]
                # remember to check for cysteine bridges
                if center_atom.element == 'S' and bond_atom.element == 'S':
                    center_atom.cysteine_bridge = True
                    bond_atom.cysteine_bridge = True
                # set up bonding
                if not bond_atom in center_atom.bonded_atoms:
                    center_atom.bonded_atoms.append(bond_atom)
                if not center_atom in bond_atom.bonded_atoms:
                    bond_atom.bonded_atoms.append(center_atom)
        # found info on covalent coupling
        if tag == 'CCOUPL' and len(line) > 14:
            conect_numbers = [line[i:i+5] for i in range(6, len(line)-1, 5)]
            center_atom = atoms[int(conect_numbers[0])]
            for num in conect_numbers[1:]:
                cov_atom = atoms[int(num)]
                center_atom.group.couple_covalently(cov_atom.group)
        # found info on non-covalent coupling
        if tag == 'NCOUPL' and len(line) > 14:
            conect_numbers = [line[i:i+5] for i in range(6, len(line)-1, 5)]
            center_atom = atoms[int(conect_numbers[0])]
            for num in conect_numbers[1:]:
                cov_atom = atoms[int(num)]
                center_atom.group.couple_non_covalently(cov_atom.group)
        # this conformation is done - yield the atoms
        if tag == 'ENDMDL':
            for num in numbers:
                yield (conformation, atoms[num])
            # prepare for next conformation
            atoms = {}
            numbers = []
