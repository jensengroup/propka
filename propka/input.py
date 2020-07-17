"""
Input handling
==============

Input routines.
"""
from pathlib import Path
from pkg_resources import resource_filename
from propka.lib import protein_precheck
from propka.atom import Atom
from propka.conformation_container import ConformationContainer
from propka.group import initialize_atom_group


def open_file_for_reading(input_file):
    """Open file or file-like stream for reading.

    TODO - convert this to a context manager

    Args:
        input_file: path to file or file-like object. If file-like object,
        then will attempt seek(0).
    """
    try:
        input_file.seek(0)
        return input_file
    except AttributeError:
        pass

    try:
        file_ = open(input_file, 'rt')
    except:
        raise IOError('Cannot find file {0:s}'.format(input_file))
    return file_


def read_molecule_file(filename: str, mol_container, stream=None):
    """Read input file or stream (PDB or PROPKA) for a molecular container

    Args:
        filename(str):  name of input file. If not using a filestream via the
            ``stream`` argument, should be a path to the file to be read.
        mol_container:  :class:`~propka.molecular_container.MolecularContainer`
            object.
        stream: optional filestream handle. If ``None``, then open
            ``filename`` as a local file for reading.

    Returns:
        updated :class:`~propka.molecular_container.MolecularContainer` object.

    Raises:
        ValuError: if invalid input given

    Examples:
        There are two main cases for using ``read_molecule_file``. The first
        (and most common) is to pass the input file (``filename``) as a
        string which gives the path of the molecule file to be read (here we
        also pass a :class:`~propka.molecular_container.MolecularContainer`
        object named ``mol_container``).

        >>> read_molecule_file('test.pdb', mol_container)
        <propka.molecular_container.MolecularContainer at 0x7f6e0c8f2310>

        The other use case is when passing a file-like object, e.g. a
        :class:`io.StringIO` class, instance. This is done by passing the
        object via the ``stream`` argument. Since file-like objects do not
        usually have an associated file name, an appropirate file name should
        be passed to the ``filename`` argument. In this case, ``filename`` is
        not opened for reading, but instead is used to help recognise the file
        type (based on the extension being either `.pdb` or `.propka_input`)
        and also uses that given ``filename`` to assign a name to the input
        :class:`~propka.molecular_container.MolecularContainer` object.

        >>> read_molecule_file('test.pdb', mol_container,
                               stream=string_io_object)
        <propka.molecular_container.MolecularContainer at 0x7f6e0c8f2310>

    """
    input_path = Path(filename)
    mol_container.name = input_path.stem
    input_file_extension = input_path.suffix

    if stream is not None:
        input_file = stream
    else:
        input_file = filename

    if input_file_extension.lower() == '.pdb':
        # input is a pdb file. read in atoms and top up containers to make
        # sure that all atoms are present in all conformations
        conformations, conformation_names = read_pdb(
            input_file, mol_container.version.parameters, mol_container)
        if len(conformations) == 0:
            str_ = ('Error: The pdb file does not seem to contain any '
                    'molecular conformations')
            raise ValueError(str_)
        mol_container.conformations = conformations
        mol_container.conformation_names = conformation_names
        mol_container.top_up_conformations()
        # make a structure precheck
        protein_precheck(
            mol_container.conformations, mol_container.conformation_names)
        # set up atom bonding and protonation
        mol_container.version.setup_bonding_and_protonation(mol_container)
        # Extract groups
        mol_container.extract_groups()
        # sort atoms
        for name in mol_container.conformation_names:
            mol_container.conformations[name].sort_atoms()
        # find coupled groups
        mol_container.find_covalently_coupled_groups()
    elif input_file_extension.lower() == '.propka_input':
        # input is a propka_input file
        conformations, conformation_names = read_propka(
            input_file, mol_container.version.parameters, mol_container)
        mol_container.conformations = conformations
        mol_container.conformation_names = conformation_names
        # Extract groups - this merely sets up the groups found in the
        # input file
        mol_container.extract_groups()
        # do some additional set up
        mol_container.additional_setup_when_reading_input_file()
    else:
        str_ = "Unknown input file type {0!s} for file {1!s}".format(
            input_file_extension, input_path)
        raise ValueError(str_)
    return mol_container


def read_parameter_file(input_file, parameters):
    """Read a parameter file.

    Args:
        input_file:  input file to read
        parameters:  Parameters object
    Returns:
        updated Parameters object
    """
    # try to locate the parameter file
    try:
        ifile = resource_filename(__name__, input_file)
        input_ = open_file_for_reading(ifile)
    except (IOError, FileNotFoundError, ValueError, KeyError):
        input_ = open_file_for_reading(input_file)
    for line in input_:
        parameters.parse_line(line)
    return parameters


def conformation_sorter(conf):
    """TODO - figure out what this function does."""
    model = int(conf[:-1])
    altloc = conf[-1:]
    return model*100+ord(altloc)


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
    lines = open_file_for_reading(pdb_file).readlines()
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


def read_propka(input_file, parameters, molecule):
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
    names = sorted(conformations.keys(), key=conformation_sorter)
    return [conformations, names]


def get_atom_lines_from_input(input_file, tags=['ATOM  ', 'HETATM']):
    """Get atom lines from a PROPKA input file.

    Args:
        input_file:  input file
        tags:  tags defining atom lines
    Yields:
        conformation container, list of atoms
    """
    lines = open_file_for_reading(input_file).readlines()
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
            initialize_atom_group(atom)
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
    names = sorted(conformations.keys(), key=conformation_sorter)
    return [conformations, names]
