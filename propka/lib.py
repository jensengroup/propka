"""Implements many of the main functions used to call PROPKA."""
import sys
import logging
import argparse
import pkg_resources


_LOGGER = logging.getLogger("propka")
_STDOUT_HANDLER = logging.StreamHandler(sys.stdout)
_STDOUT_HANDLER.setFormatter(logging.Formatter("%(message)s"))
_LOGGER.addHandler(_STDOUT_HANDLER)


def open_file_for_reading(input_file):
    """Open file or file-like stream for reading.

    TODO - convert this to a context manager

    Args:
        input_file: path to file or file-like object. If file-like object,
        then will attempt fseek(0).
    """
    try:
        input_file.fseek(0)
        return input_file
    except AttributeError:
        pass

    try:
        file_ = open(input_file, 'rt')
    except:
        raise IOError('Cannot find file {0:s}'.format(input_file))
    return file_


def open_file_for_writing(input_file):
    """Open file or file-like stream for writing.

    TODO - convert this to a context manager.

    Args:
        input_file: path to file or file-like object. If file-like object,
        then will attempt to get file mode.
    """
    try:
        mode = input_file.mode
        if not ("w" in mode or "a" in mode or "+" in mode):
            raise IOError("File/stream not open for writing")
        return input_file
    except AttributeError:
        pass
    try:
        file_ = open(input_file, 'wt')
    except FileNotFoundError:
        raise Exception('Could not open {0:s}'.format(input_file))
    return file_


def conformation_sorter(conf):
    """TODO - figure out what this function does."""
    model = int(conf[:-1])
    altloc = conf[-1:]
    return model*100+ord(altloc)


def split_atoms_into_molecules(atoms):
    """Maps atoms into molecules.

    Args:
        atoms:  list of atoms
    Returns:
        list of molecules
    """
    molecules = []
    while len(atoms) > 0:
        initial_atom = atoms.pop()
        molecules.append(make_molecule(initial_atom, atoms))
    return molecules


def make_molecule(atom, atoms):
    """Make a molecule from atoms.

    Args:
        atom:  one of the atoms
        atoms:  a list of the remaining atoms
    Return:
        list of atoms
    """
    bonded_atoms = [a for a in atoms if atom in a.bonded_atoms]
    res_atoms = [atom,]
    for bond_atom in bonded_atoms:
        if bond_atom in atoms:
            atoms.remove(bond_atom)
            res_atoms.extend(make_molecule(bond_atom, atoms))
    return res_atoms


def make_grid(min_, max_, step):
    """Make a grid across the specified tange.

    TODO - figure out if this duplicates existing generators like `range` or
    numpy function.

    Args:
        min_:  minimum value of grid
        max_:  maximum value of grid
        step:  grid step size
    """
    x = min_
    while x <= max_:
        yield x
        x += step


def generate_combinations(interactions):
    """Generate combinations of interactions.

    Args:
        interactions:  list of interactions
    Returns:
        list of combinations
    """
    res = [[]]
    for interaction in interactions:
        res = make_combination(res, interaction)
    res.remove([])
    return res


def make_combination(combis, interaction):
    """Make a specific set of combinations.

    Args:
        combis:  list of combinations
        interaction:  interaction to add to combinations
    Returns:
        list of combinations
    """
    res = []
    for combi in combis:
        res.append(combi+[interaction])
        res.append(combi)
    return res


def parse_res_string(res_str):
    """Parse a residue string.

    Args:
        res_string:  residue string in format "chain:resnum[inscode]"
    Returns:
        a tuple of (chain, resnum, inscode).
    Raises:
        ValueError if the input string is invalid.
    """
    try:
        chain, resnum_str = res_str.split(":")
    except ValueError:
        raise ValueError("Invalid residue string (must contain 2 "
                         "colon-separated values)")
    try:
        resnum = int(resnum_str)
    except ValueError:
        try:
            resnum = int(resnum_str[:-1])
        except ValueError:
            raise ValueError("Invalid residue number (not an int)")
        else:
            inscode = resnum_str[-1]
    else:
        inscode = " "
    return chain, resnum, inscode


def build_parser(parser=None):
    """Build an argument parser for PROPKA.

    Args:
        parser:  existing parser. If this is not None, then the PROPKA parser will
                 be created as a subparser to this existing parser.  Otherwise, a
                 new parser will be created.
    Returns:
        ArgumentParser object.
    """
    if parser is not None:
        group = parser.add_argument_group(title="PROPKA invoation options")
    else:
        parser = argparse.ArgumentParser(
            description=("PROPKA predicts the pKa values of ionizable "
                         "groups in proteins and protein-ligand "
                         "complexes based in the 3D structure"),
            formatter_class=argparse.ArgumentDefaultsHelpFormatter)
        # This is duck-typing at its finest
        group = parser
        group.add_argument("input_pdb", help="read data from <filename>")
    group.add_argument(
        "-f", "--file", action="append", dest="filenames", default=[],
        help="read data from <filename>, i.e. <filename> is added to arguments")
    group.add_argument(
        "-r", "--reference", dest="reference", default="neutral",
        help=("setting which reference to use for stability calculations "
              "[neutral/low-pH]"))
    group.add_argument(
        "-c", "--chain", action="append", dest="chains",
        help=('creating the protein with only a specified chain. Specify '
              '" " for chains without ID [all]'))
    group.add_argument(
        "-i", "--titrate_only", dest="titrate_only",
        help=('Treat only the specified residues as titratable. Value should '
              'be a comma-separated list of "chain:resnum" values; for example: '
              '-i "A:10,A:11"'))
    group.add_argument(
        "-t", "--thermophile", action="append", dest="thermophiles",
        help=("defining a thermophile filename; usually used in "
              "'alignment-mutations'"))
    group.add_argument(
        "-a", "--alignment", action="append", dest="alignment",
        help=("alignment file connecting <filename> and <thermophile> "
              "[<thermophile>.pir]"))
    group.add_argument(
        "-m", "--mutation", action="append", dest="mutations",
        help=("specifying mutation labels which is used to modify "
              "<filename> according to, e.g. N25R/N181D"))
    group.add_argument(
        "-v", "--version", dest="version_label", default="Jan15",
        help="specifying the sub-version of propka [Jan15/Dec19]")
    group.add_argument(
        "-p", "--parameters", dest="parameters",
        default=pkg_resources.resource_filename(__name__, "propka.cfg"),
        help="set the parameter file [{default:s}]")
    try:
        group.add_argument(
            "--log-level",
            choices=["DEBUG", "INFO", "WARNING", "ERROR", "CRITICAL"],
            help="logging level verbosity", default="INFO")
    except argparse.ArgumentError:
        # It is possible that --log-level has already been set by APBS
        pass
    group.add_argument(
        "-o", "--pH", dest="pH", type=float, default=7.0,
        help="setting pH-value used in e.g. stability calculations [7.0]")
    group.add_argument(
        "-w", "--window", dest="window", nargs=3, type=float,
        default=(0.0, 14.0, 1.0),
        help=("setting the pH-window to show e.g. stability profiles "
              "[0.0, 14.0, 1.0]"))
    group.add_argument(
        "-g", "--grid", dest="grid", nargs=3, type=float,
        default=(0.0, 14.0, 0.1),
        help=("setting the pH-grid to calculate e.g. stability "
              "related properties [0.0, 14.0, 0.1]"))
    group.add_argument(
        "--mutator", dest="mutator",
        help="setting approach for mutating <filename> [alignment/scwrl/jackal]")
    group.add_argument(
        "--mutator-option", dest="mutator_options", action="append",
        help="setting property for mutator [e.g. type=\"side-chain\"]")
    group.add_argument(
        "-d", "--display-coupled-residues", dest="display_coupled_residues",
        action="store_true",
        help=("Displays alternative pKa values due "
              "to coupling of titratable groups"))
    group.add_argument(
        "-l", "--reuse-ligand-mol2-files", dest="reuse_ligand_mol2_file",
        action="store_true", default=False,
        help=("Reuses the ligand mol2 files allowing the user to alter "
              "ligand bond orders"))
    group.add_argument(
        "-k", "--keep-protons", dest="keep_protons", action="store_true",
        help="Keep protons in input file", default=False)
    group.add_argument(
        "-q", "--quiet", action="store_const", const="WARNING",
        dest="log_level", help="suppress non-warning messages")
    group.add_argument(
        "--protonate-all", dest="protonate_all", action="store_true",
        help="Protonate all atoms (will not influence pKa calculation)",
        default=False)
    return parser


def loadOptions(args):
    """Load the arguments parser with options.

    NOTE - verbosity is set as soon as this function is invoked.

    Arguments:
        args:  list of arguments
    Returns:
        argparse namespace
    """
    # loading the parser
    parser = build_parser()
    # parsing and returning options and arguments
    if len(args) == 0:
        # command line
        options = parser.parse_args()
    else:
        options = parser.parse_args(args)
    # adding specified filenames to arguments
    options.filenames.append(options.input_pdb)
    # Convert titrate_only string to a list of (chain, resnum) items:
    if options.titrate_only is not None:
        res_list = []
        for res_str in options.titrate_only.split(','):
            try:
                chain, resnum, inscode = parse_res_string(res_str)
            except ValueError:
                _LOGGER.critical(
                    'Invalid residue string: "{0:s}"'.format(res_str))
                sys.exit(1)
            res_list.append((chain, resnum, inscode))
        options.titrate_only = res_list
    # Set the no-print variable
    level = getattr(logging, options.log_level)
    _LOGGER.setLevel(level)
    # done!
    return options


def make_tidy_atom_label(name, element):
    """Returns a 'tidier' atom label for printing to the new PDB file.

    Args:
        name:  atom name
        element:  atom element
    Returns:
        string
    """
    if len(name) > 4: # if longer than 4, just truncate the name
        label = name[0:4]
    elif len(name) == 4: # if length is 4, otherwise use the name as it is
        label = name
    else: # if less than 4 characters long, insert white space as needed
        if len(element) == 1:
            label = ' {0:<3s}'.format(name)
        else: # The element should occupy the two first chars
            label = '{0:<4s}'.format(name)
    return label


def get_sorted_configurations(configuration_keys):
    """Extract and sort configurations.

    Args:
        configuration_keys:  list of configuration keys
    Returns:
        list of configurations
    """
    configurations = list(configuration_keys)
    configurations.sort(key=configuration_compare)
    return configurations


def configuration_compare(conf):
    """TODO - figure out what this function does."""
    return 100*int(conf[1:-2]) + ord(conf[-1])


def write_file(filename, lines):
    """Writes a new file.

    Args:
        filename:  name of file
        lines:  lines to write to file
    """
    file_ = open_file_for_writing(filename)
    for line in lines:
        file_.write("{0:s}\n".format(line))
    file_.close()


def _args_to_str(arg_list):
    """Summarize list of arguments in string.

    Args:
        arg_list:  list of arguments
    Returns:
        string
    """
    return " ".join(map(str, arg_list))


def info(*args):
    """Log a message to info.

    Level defaults to INFO unless overridden.

    Args:
        args:  argument list
    """
    _LOGGER.info(_args_to_str(args))


def debug(*args):
    """Log a message to debug.

    Level defaults to DEBUG unless overridden.

    Args:
        args:  argument list
    """
    _LOGGER.debug(_args_to_str(args))


def warning(*args):
    """Log a message to warning.

    Level defaults to WARNING unless overridden.

    Args:
        args:  argument list
    """
    _LOGGER.warning(_args_to_str(args))
