"""
Set-up of a PROPKA calculation
==============================

Implements many of the main functions used to call PROPKA.
"""

import logging
import argparse
from pathlib import Path
from typing import Iterable, Iterator, List, TYPE_CHECKING, NoReturn, Optional, Tuple, TypeVar

if TYPE_CHECKING:
    from propka.atom import Atom


T = TypeVar("T")
Number = TypeVar("Number", int, float)

_T_RESIDUE_TUPLE = Tuple[str, int, str]

_LOGGER = logging.getLogger(__name__)


EXPECTED_ATOM_NUMBERS = {'ALA': 5, 'ARG': 11, 'ASN': 8, 'ASP': 8, 'CYS': 6,
                         'GLY': 4, 'GLN': 9, 'GLU': 9, 'HIS': 10, 'ILE': 8,
                         'LEU': 8, 'LYS': 9, 'MET': 8, 'PHE': 11, 'PRO': 7,
                         'SER': 6, 'THR': 7, 'TRP': 14, 'TYR': 12, 'VAL': 7}


class Options:
    # Note: All the "NoReturn" members appear to be unused
    alignment: NoReturn  # Optional[List[str]]
    chains: Optional[List[str]]
    display_coupled_residues: bool = False
    filenames: List[str]  # List[Path]?
    grid: Tuple[float, float, float] = (0.0, 14.0, 0.1)
    input_pdb: str  # Path?
    keep_protons: bool = False
    log_level: str = 'INFO'
    mutations: NoReturn  # Optional[List[str]]
    mutator: NoReturn  # Optional[str]  # alignment/scwrl/jackal
    mutator_options: NoReturn  # Optional[List[str]]
    pH: NoReturn  # float = 7.0
    parameters: Path
    protonate_all: bool = False
    reference: NoReturn  # str = 'neutral'
    reuse_ligand_mol2_file: bool = False  # only used by unused function
    thermophiles: NoReturn  # Optional[List[str]]
    titrate_only: Optional[List[_T_RESIDUE_TUPLE]]
    window: Tuple[float, float, float] = (0.0, 14.0, 1.0)


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
                    _LOGGER.warning(str_)
                continue
            # check number of atoms in residue
            if len(res_atoms) != EXPECTED_ATOM_NUMBERS[res_name]:
                str_ = ("Unexpected number ({num:d}) of atoms in residue "
                        "{res:s} in conformation {conf:s}".format(
                            num=len(res_atoms), res=residue_label,
                            conf=name))
                _LOGGER.warning(str_)


def resid_from_atom(atom):
    """Return string with atom residue information.

    Args:
        atom:  atom to generate string for
    Returns
        string
    """
    return '{0:>4d} {1:s} {2:s}'.format(
        atom.res_num, atom.chain_id, atom.icode)


def split_atoms_into_molecules(atoms: List["Atom"]):
    """Maps atoms into molecules.

    Args:
        atoms:  list of atoms
    Returns:
        list of molecules
    """
    molecules: List[List["Atom"]] = []
    while len(atoms) > 0:
        initial_atom = atoms.pop()
        molecules.append(make_molecule(initial_atom, atoms))
    return molecules


def make_molecule(atom: "Atom", atoms: List["Atom"]):
    """Make a molecule from atoms.

    Args:
        atom:  one of the atoms
        atoms:  a list of the remaining atoms
    Return:
        list of atoms
    """
    bonded_atoms = [a for a in atoms if atom in a.bonded_atoms]
    res_atoms = [atom]
    for bond_atom in bonded_atoms:
        if bond_atom in atoms:
            atoms.remove(bond_atom)
            res_atoms.extend(make_molecule(bond_atom, atoms))
    return res_atoms


def make_grid(min_: Number, max_: Number, step: Number) -> Iterator[Number]:
    """Make a grid across the specified tange.

    Like range() for integers or numpy.arange() for floats, except that `max_`
    is not excluded from the range.

    Args:
        min_:  minimum value of grid
        max_:  maximum value of grid
        step:  grid step size
    """
    x = min_
    while x <= max_:
        yield x
        x += step


def generate_combinations(interactions: Iterable[T]) -> List[List[T]]:
    """Generate combinations of interactions.

    Args:
        interactions:  list of interactions
    Returns:
        list of combinations
    """
    res: List[List[T]] = [[]]
    for interaction in interactions:
        res = make_combination(res, interaction)
    res.remove([])
    return res


def make_combination(combis: List[List[T]], interaction: T) -> List[List[T]]:
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


def parse_res_string(res_str: str) -> _T_RESIDUE_TUPLE:
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


def parse_res_list(titrate_only: str):
    res_list: List[_T_RESIDUE_TUPLE] = []
    for res_str in titrate_only.split(','):
        try:
            res_list.append(parse_res_string(res_str))
        except ValueError as ex:
            raise argparse.ArgumentTypeError(f'{ex}: "{res_str:s}"')
    return res_list


def build_parser(parser=None):
    """Build an argument parser for PROPKA.

    Args:
        parser:  existing parser. If this is not None, then the PROPKA parser
                 will be created as a subparser to this existing parser.
                 Otherwise, a new parser will be created.
    Returns:
        ArgumentParser object.


    .. versionchanged:: 3.4.0
       Argument `--generate-propka-input` has been removed as writing PROPKA
       input files is no longer supported.
    """
    import propka

    if parser is not None:
        group = parser.add_argument_group(title="PROPKA invocation options")
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
        help=(
            "read data from <filename>, i.e. <filename> is added to arguments"
        )
    )
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
        type=parse_res_list,
        help=('Treat only the specified residues as titratable. Value should '
              'be a comma-separated list of "chain:resnum" values; for '
              'example: -i "A:10,A:11"'))
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
        "--version", action="version", version=f"%(prog)s {propka.__version__}")
    group.add_argument(
        "-p", "--parameters", dest="parameters",
        type=Path, default=Path(__file__).parent / "propka.cfg",
        help="set the parameter file")
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
        help=(
            "setting approach for mutating <filename> "
            "[alignment/scwrl/jackal]"
        )
    )
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


def loadOptions(args=None) -> Options:
    """
    Load the arguments parser with options. Note that verbosity is set as soon
    as this function is invoked.

    Arguments:
        args:  list of arguments
    Returns:
        argparse namespace
    """
    # loading the parser
    parser = build_parser()
    # parsing and returning options and arguments
    options = parser.parse_args(args, namespace=Options())

    # adding specified filenames to arguments
    options.filenames.append(options.input_pdb)
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
    if len(name) > 4:  # if longer than 4, just truncate the name
        label = name[0:4]
    elif len(name) == 4:  # if length is 4, otherwise use the name as it is
        label = name
    else:  # if less than 4 characters long, insert white space as needed
        if len(element) == 1:
            label = ' {0:<3s}'.format(name)
        else:  # The element should occupy the two first chars
            label = '{0:<4s}'.format(name)
    return label
