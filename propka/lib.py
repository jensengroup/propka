from __future__ import division
from __future__ import print_function

import sys
import pkg_resources
import logging
import argparse


logger = logging.getLogger("propka")
stdout_handler = logging.StreamHandler(sys.stdout)
stdout_handler.setFormatter(logging.Formatter("%(message)s"))
logger.addHandler(stdout_handler)

#
# file I/O
#
def open_file_for_reading(filename):
    """Open file or file-like stream  *filename* for reading.

    *filename* may be a string and then it is opened but if it is a
    file-like object (such as an open :class:`file` or
    :class:`StringIO.StringIO` --- really anything with ``next()``,
    ``read()``, ``readlines()``, ``readline``, ``close`` methods) then
    the object is just passed through (the stream is attempted to be
    reset to the beginning with ``fseek(0)``).
    """
    if (hasattr(filename, 'next') or hasattr(filename, '__next__')) \
        and hasattr(filename, 'read') \
        and hasattr(filename, 'readline') and hasattr(filename, 'readlines') \
        and hasattr(filename, 'close'):
        # already a stream
        try:
            filename.fseek(0)
        except AttributeError:
            pass
        return filename

    try:
        f = open(filename,'r')
    except:
        raise IOError('Cannot find file %s' %filename)
    return f

def open_file_for_writing(filename):
    """Open file or file-like stream for writing"""
    if hasattr(filename, 'write') and hasattr(filename, 'writeline') and hasattr(filename, 'writelines') \
            and hasattr(filename, 'close'):
        # already a stream
        try:
            mode = filename.mode
        except AttributeError:
            mode = "w"
        if not ("w" in mode or "a" in mode or "+" in mode):
            raise IOError("File/stream not open for writing")
        return filename

    try:
        f = open(filename,'w')
    except:
        raise Exception('Could not open %s'%filename)
    return f

#
# bookkeeping etc.
#
def conformation_sorter(conf):
    model = int(conf[:-1])
    altloc = conf[-1:]
    return model*100+ord(altloc)

def split_atoms_into_molecules(atoms):
    molecules = []

    while len(atoms)>0:
        initial_atom = atoms.pop()
        molecules.append( make_molecule(initial_atom,atoms))

    return molecules

def make_molecule(atom, atoms):
    bonded_atoms = [a for a in atoms if atom in a.bonded_atoms]
    res_atoms = [atom,]

    for ba in bonded_atoms:
        if ba in atoms:
            atoms.remove(ba)
            res_atoms.extend(make_molecule(ba, atoms))

    return res_atoms


def make_grid(min,max,step):
    x = min
    while x <= max:
        yield x
        x += step
    return

def generate_combinations(interactions):
    res = [[]]
    for interaction in interactions:
        res = make_combination(res, interaction)
    res.remove([])

    return res


def make_combination(combis, interaction):
    res = []
    for combi in combis:
        res.append(combi+[interaction])
        res.append(combi)
    return res


def parse_res_string(res_str):
    """
    Parse the residue string, in format "chain:resnum[inscode]", and return
    a tuple of (chain, resnum, inscode). Raises ValueError if the input
    string is invalid.
    """
    try:
        chain, resnum_str = res_str.split(":")
    except ValueError:
        raise ValueError("Invalid residue string (must contain 2 colon-separated values)")
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
        parser_:  existing parser.  If this is not None, then the PROPKA parser will
                 be created as a subparser to this existing parser.  Otherwise, a 
                 new parser will be created.
    Returns:
        ArgumentParser object.
    """
    if parser is not None:
        group = parser.add_argument_group(title="PROPKA invoation options")
    else:
        parser = argparse.ArgumentParser(description=("PROPKA predicts the pKa values of ionizable "
                                                      "groups in proteins and protein-ligand "
                                                      "complexes based in the 3D structure"),
                                         formatter_class=argparse.ArgumentDefaultsHelpFormatter)
        # This is duck-typing at its finest
        group = parser
        group.add_argument("input_pdb", help="read data from <filename>")

    group.add_argument("-f", "--file", action="append", dest="filenames", default=[],
                       help="read data from <filename>, i.e. <filename> is added to arguments")
    group.add_argument("-r", "--reference", dest="reference", default="neutral",
                       help=("setting which reference to use for stability calculations "
                             "[neutral/low-pH]"))
    group.add_argument("-c", "--chain", action="append", dest="chains",
                       help=('creating the protein with only a specified chain. Specify '
                             '" " for chains without ID [all]'))
    group.add_argument("-i", "--titrate_only", dest="titrate_only",
                       help=('Treat only the specified residues as titratable. Value should '
                             'be a comma-separated list of "chain:resnum" values; for example: '
                             '-i "A:10,A:11"'))
    group.add_argument("-t", "--thermophile", action="append", dest="thermophiles",
                       help=("defining a thermophile filename; usually used in "
                             "'alignment-mutations'"))
    group.add_argument("-a", "--alignment", action="append", dest="alignment",
                       help=("alignment file connecting <filename> and <thermophile> "
                             "[<thermophile>.pir]"))
    group.add_argument("-m", "--mutation", action="append", dest="mutations",
                       help=("specifying mutation labels which is used to modify "
                             "<filename> according to, e.g. N25R/N181D"))
    group.add_argument("-v", "--version", dest="version_label", default="Jan15",
                       help="specifying the sub-version of propka [Jan15/Dec19]")
    group.add_argument("-p", "--parameters", dest="parameters",
                       default=pkg_resources.resource_filename(__name__, "propka.cfg"),
                       help="set the parameter file [%(default)s]")
    try:
        group.add_argument("--log-level", choices=["DEBUG", "INFO", "WARNING", "ERROR", "CRITICAL"],
                           help="logging level verbosity", default="INFO")
    except argparse.ArgumentError:
        # It is possible that --log-level has already been set by APBS
        pass
    group.add_argument("-o", "--pH", dest="pH", type=float, default=7.0,
                       help="setting pH-value used in e.g. stability calculations [7.0]")
    group.add_argument("-w", "--window", dest="window", nargs=3, type=float,
                       default=(0.0, 14.0, 1.0),
                       help=("setting the pH-window to show e.g. stability profiles "
                             "[0.0, 14.0, 1.0]"))
    group.add_argument("-g", "--grid", dest="grid", nargs=3, type=float,
                       default=(0.0, 14.0, 0.1),
                       help=("setting the pH-grid to calculate e.g. stability "
                             "related properties [0.0, 14.0, 0.1]"))
    group.add_argument("--mutator", dest="mutator",
                       help="setting approach for mutating <filename> [alignment/scwrl/jackal]")
    group.add_argument("--mutator-option", dest="mutator_options", action="append",
                       help="setting property for mutator [e.g. type=\"side-chain\"]")
    group.add_argument("-d", "--display-coupled-residues", dest="display_coupled_residues",
                       action="store_true", help=("Displays alternative pKa values due "
                                                  "to coupling of titratable groups"))
    group.add_argument("-l", "--reuse-ligand-mol2-files", dest="reuse_ligand_mol2_file",
                       action="store_true", default=False,
                       help=("Reuses the ligand mol2 files allowing the user to alter "
                             "ligand bond orders"))
    group.add_argument("-k", "--keep-protons", dest="keep_protons", action="store_true",
                       help="Keep protons in input file", default=False)
    group.add_argument("-q", "--quiet", action="store_const", const="WARNING",
                       dest="log_level", help="suppress non-warning messages")
    group.add_argument("--protonate-all", dest="protonate_all", action="store_true",
                       help="Protonate all atoms (will not influence pKa calculation)",
                       default=False)
    return parser


def loadOptions(args):
    """
    Load the arguments parser with options. Note that verbosity is set as soon
    as this function is invoked.

    Arguments:
        args:  list of arguments
    Returns:
        argparse namespace
    """
    # defining a 'usage' message
    usage = "usage: %prog [options] filename"

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
                logger.critical('Invalid residue string: "%s"' % res_str)
                sys.exit(1)
            res_list.append((chain, resnum, inscode))
        options.titrate_only = res_list


    # Set the no-print variable
    level = getattr(logging, options.log_level)
    logger.setLevel(level)

    # done!
    return options


def makeTidyAtomLabel(name,element):
    """
    Returns a 'tidier' atom label for printing the new pdbfile
    """

    if len(name)>4:# if longer than 4, just truncate the name
        label=name[0:4]
    elif len(name)==4:# if lenght is 4, otherwise use the name as it is
        label = name
    else: # if less than 4 characters long, insert white space as needed
        if len(element)==1:
            label = ' %-3s'%name
        else: # The element shoul occupy the two first chars
            label = '%-4s'%name

    return label



def get_sorted_configurations(configuration_keys):
    """
    extract and sort configurations
    """
    configurations = list(configuration_keys)
    configurations.sort(key=configuration_compare)
    return configurations

def configuration_compare(conf):
    return 100*int(conf[1:-2]) + ord(conf[-1])




def writeFile(filename, lines):
    """
    Writes a new file
    """
    f = open_file_for_writing(filename)

    for line in lines:
        f.write( "%s\n" % (line) )
    f.close()



def _args_to_str(arg_list):
    return " ".join(map(str, arg_list))

def info(*args):
    """Log a message. Level defaults to INFO unless overridden."""
    logger.info(_args_to_str(args))

def debug(*args):
    """Log a message on the DEBUG level."""
    logger.debug(_args_to_str(args))

def warning(*args):
    """Log a WARN message"""
    logger.warning(_args_to_str(args))
