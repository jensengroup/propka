"""
Configuration file parameters
=============================

Holds parameters and settings that can be set in :file:`propka.cfg`. The file
format consists of lines of  ``keyword value [value ...]``, blank lines, and
comment lines (introduced with ``#``).

The module attributes below list the names and types of all key words
in configuration file.

"""
import logging
from dataclasses import dataclass, field
from typing import Dict, List

try:
    # New in version 3.10, deprecated since version 3.12
    from typing import TypeAlias
except ImportError:
    TypeAlias = "TypeAlias"  # type: ignore


_LOGGER = logging.getLogger(__name__)


class squared_property:

    def __set_name__(self, owner, name: str):
        assert name.endswith("_squared")
        self._name_not_squared = name[:-len("_squared")]  # removesuffix()

    def __get__(self, instance, owner=None) -> float:
        if instance is None:
            return self  # type: ignore[return-value]
        return getattr(instance, self._name_not_squared)**2

    def __set__(self, instance, value: float):
        setattr(instance, self._name_not_squared, value**0.5)


_T_MATRIX: TypeAlias = "InteractionMatrix"
_T_PAIR_WISE_MATRIX: TypeAlias = "PairwiseMatrix"
_T_NUMBER_DICTIONARY = Dict[str, float]
_T_LIST_DICTIONARY = Dict[str, list]
_T_STRING_DICTIONARY = Dict[str, str]
_T_STRING_LIST = List[str]
_T_STRING = str
_T_BOOL = bool


@dataclass
class Parameters:
    """PROPKA parameter class."""

    # MATRICES
    interaction_matrix: _T_MATRIX = field(
        default_factory=lambda: InteractionMatrix("interaction_matrix"))

    # PAIR_WISE_MATRICES
    sidechain_cutoffs: _T_PAIR_WISE_MATRIX = field(
        default_factory=lambda: PairwiseMatrix("sidechain_cutoffs"))

    # NUMBER_DICTIONARIES
    VanDerWaalsVolume: _T_NUMBER_DICTIONARY = field(default_factory=dict)
    charge: _T_NUMBER_DICTIONARY = field(default_factory=dict)
    model_pkas: _T_NUMBER_DICTIONARY = field(default_factory=dict)
    ions: _T_NUMBER_DICTIONARY = field(default_factory=dict)
    valence_electrons: _T_NUMBER_DICTIONARY = field(default_factory=dict)
    custom_model_pkas: _T_NUMBER_DICTIONARY = field(default_factory=dict)

    # LIST_DICTIONARIES
    backbone_NH_hydrogen_bond: _T_LIST_DICTIONARY = field(default_factory=dict)
    backbone_CO_hydrogen_bond: _T_LIST_DICTIONARY = field(default_factory=dict)

    # STRING_DICTIONARIES
    protein_group_mapping: _T_STRING_DICTIONARY = field(default_factory=dict)

    # STRING_LISTS
    ignore_residues: _T_STRING_LIST = field(default_factory=list)
    angular_dependent_sidechain_interactions: _T_STRING_LIST = field(default_factory=list)
    acid_list: _T_STRING_LIST = field(default_factory=list)
    base_list: _T_STRING_LIST = field(default_factory=list)
    exclude_sidechain_interactions: _T_STRING_LIST = field(default_factory=list)
    backbone_reorganisation_list: _T_STRING_LIST = field(default_factory=list)
    write_out_order: _T_STRING_LIST = field(default_factory=list)

    # DISTANCES
    desolv_cutoff: float = 20.0
    buried_cutoff: float = 15.0
    coulomb_cutoff1: float = 4.0
    coulomb_cutoff2: float = 10.0

    # DISTANCES SQUARED
    desolv_cutoff_squared = squared_property()
    buried_cutoff_squared = squared_property()
    coulomb_cutoff1_squared = squared_property()
    coulomb_cutoff2_squared = squared_property()

    # STRINGS
    version: _T_STRING = "VersionA"
    output_file_tag: _T_STRING = ""
    ligand_typing: _T_STRING = "groups"
    pH: _T_STRING = "variable"
    reference: _T_STRING = "neutral"

    # PARAMETERS
    Nmin: int = 280
    Nmax: int = 560
    desolvationSurfaceScalingFactor: float = 0.25
    desolvationPrefactor: float = -13.0
    desolvationAllowance: float = 0.0
    coulomb_diel: float = 80.0
    # TODO - it would be nice to rename these; they're defined everywhere
    COO_HIS_exception: float = 1.60
    OCO_HIS_exception: float = 1.60
    CYS_HIS_exception: float = 1.60
    CYS_CYS_exception: float = 3.60
    min_ligand_model_pka: float = -10.0
    max_ligand_model_pka: float = 20.0
    # include_H_in_interactions: NoReturn = None
    coupling_max_number_of_bonds: int = 3
    min_bond_distance_for_hydrogen_bonds: int = 4
    # coupling_penalty: NoReturn = None
    shared_determinants: _T_BOOL = False
    common_charge_centre: _T_BOOL = False
    # hide_penalised_group: NoReturn = None
    remove_penalised_group: _T_BOOL = True
    max_intrinsic_pka_diff: float = 2.0
    min_interaction_energy: float = 0.5
    max_free_energy_diff: float = 1.0
    min_swap_pka_shift: float = 1.0
    min_pka: float = 0.0
    max_pka: float = 10.0
    sidechain_interaction: float = 0.85

    def parse_line(self, line):
        """Parse parameter file line."""
        # first, remove comments
        comment_pos = line.find('#')
        if comment_pos != -1:
            line = line[:comment_pos]
        # split the line into words
        words = line.split()
        if len(words) == 0:
            return
        # parse the words
        typeannotation = self.__annotations__.get(words[0])
        if typeannotation is _T_NUMBER_DICTIONARY:
            self.parse_to_number_dictionary(words)
        elif typeannotation is _T_STRING_LIST:
            self.parse_to_string_list(words)
        elif typeannotation is _T_STRING:
            self.parse_string(words)
        elif typeannotation is _T_LIST_DICTIONARY:
            self.parse_to_list_dictionary(words)
        elif typeannotation is _T_MATRIX or typeannotation is _T_PAIR_WISE_MATRIX:
            self.parse_to_matrix(words)
        elif typeannotation is _T_STRING_DICTIONARY:
            self.parse_to_string_dictionary(words)
        else:
            self.parse_parameter(words)

    def parse_to_number_dictionary(self, words):
        """Parse field to number dictionary.

        Args:
            words:  strings to parse.
        """
        assert len(words) == 3, words
        dict_ = getattr(self, words[0])
        key = words[1]
        value = words[2]
        dict_[key] = float(value)

    def parse_to_string_dictionary(self, words):
        """Parse field to string dictionary.

        Args:
            words:  strings to parse
        """
        assert len(words) == 3, words
        dict_ = getattr(self, words[0])
        key = words[1]
        value = words[2]
        dict_[key] = value

    def parse_to_list_dictionary(self, words: List[str]):
        """Parse field to list dictionary.

        Args:
            words:  strings to parse.
        """
        assert len(words) > 2, words
        dict_ = getattr(self, words[0])
        key = words[1]
        if key not in dict_:
            dict_[key] = []
        for value in words[2:]:
            if isinstance(value, list):
                dict_[key].append([float(x) for x in value])
            dict_[key].append(float(value))

    def parse_to_string_list(self, words):
        """Parse field to string list.

        Args:
            words:  strings to parse
        """
        assert len(words) == 2, words
        list_ = getattr(self, words[0])
        value = words[1]
        list_.append(value)

    def parse_to_matrix(self, words):
        """Parse field to matrix.

        Args:
            words:  strings to parse
        """
        matrix = getattr(self, words[0])
        value = tuple(words[1:])
        matrix.add(value)

    def parse_parameter(self, words):
        """Parse field to parameters.

        Args:
            words:  strings to parse
        """
        assert len(words) == 2, words
        value = float(words[1])
        setattr(self, words[0], value)

    def parse_string(self, words):
        """Parse field to strings.

        Args:
            words:  strings to parse
        """
        assert len(words) == 2, words
        setattr(self, words[0], words[1])

    def print_interaction_parameters(self):
        """Print interaction parameters."""
        _LOGGER.info('--------------- Model pKa values ----------------------')
        for k in self.model_pkas:
            _LOGGER.info('{0:>3s} {1:8.2f}'.format(k, self.model_pkas[k]))

        _LOGGER.info('')
        _LOGGER.info('--------------- Interactions --------------------------')
        agroups = [
            'COO', 'HIS', 'CYS', 'TYR', 'SER', 'N+', 'LYS', 'AMD', 'ARG',
            'TRP', 'ROH', 'CG', 'C2N', 'N30', 'N31', 'N32', 'N33', 'NAR',
            'OCO', 'NP1', 'OH', 'O3', 'CL', 'F', 'NAM', 'N1', 'O2', 'OP', 'SH']
        lgroups = [
            'CG', 'C2N', 'N30', 'N31', 'N32', 'N33', 'NAR', 'OCO', 'NP1',
            'OH', 'O3', 'CL', 'F', 'NAM', 'N1', 'O2', 'OP', 'SH']
        map_ = {
            'CG': ['ARG'], 'C2N': ['ARG'], 'N30': ['N+', 'LYS'],
            'N31': ['N+', 'LYS'], 'N32': ['N+', 'LYS'], 'N33': ['N+', 'LYS'],
            'NAR': ['HIS'], 'OCO': ['COO'], 'OP': [], 'SH': ['CYS'],
            'NP1': [], 'OH': ['ROH'], 'O3': [], 'CL': [], 'F': [],
            'NAM': ['AMD'], 'N1': [], 'O2': []}
        for group1 in agroups:
            for group2 in lgroups:
                fmt = "{grp1:>3s} {grp2:>3s} {mat:1s} {val1:4} {val2:4}"
                interaction = fmt.format(
                    grp1=group1, grp2=group2,
                    mat=self.interaction_matrix[group1][group2],
                    val1=self.sidechain_cutoffs.get_value(group1, group2)[0],
                    val2=self.sidechain_cutoffs.get_value(group1, group2)[1])
                map_interaction = ''
                if group2 in map_:
                    for val in map_[group2]:
                        fmt = (
                            "|{grp1:>3s} {grp2:>3s} {mat:1s} {val1:4} {val2:4}"
                        )
                        map_interaction += fmt.format(
                            group1, val, self.interaction_matrix[group1][val],
                            self.sidechain_cutoffs.get_value(group1, val)[0],
                            self.sidechain_cutoffs.get_value(group1, val)[1])
                        if (self.interaction_matrix[group1][val]
                                != self.interaction_matrix[group1][group2]):
                            map_interaction += '* '
                        if (self.sidechain_cutoffs.get_value(group1, val)[0]
                                != self.sidechain_cutoffs.get_value(
                                    group1, group2)[0]
                                or self.sidechain_cutoffs.get_value(
                                    group1, val)[1]
                                != self.sidechain_cutoffs.get_value(
                                    group1, group2)[1]):
                            map_interaction += '! '
                        else:
                            map_interaction += '  '
                    if (len(map_[group2]) == 0
                            and (self.sidechain_cutoffs.get_value(
                                group1, group2)[0]
                                 != 3
                                 or self.sidechain_cutoffs.get_value(
                                     group1, group2)[1]
                                 != 4)):
                        map_interaction += '?  '
                _LOGGER.info("%s %s", interaction, map_interaction)
                if group1 == group2:
                    break
            _LOGGER.info('-')
        _LOGGER.info('--------------- Exceptions ----------------------------')
        _LOGGER.info('COO-HIS %s', self.COO_HIS_exception)
        _LOGGER.info('OCO-HIS %s', self.OCO_HIS_exception)
        _LOGGER.info('CYS-HIS %s', self.CYS_HIS_exception)
        _LOGGER.info('CYS-CYS %s', self.CYS_CYS_exception)

        _LOGGER.info('--------------- Mapping -------------------------------')
        _LOGGER.info("""
Titratable:
CG  ARG
C2N ARG
N30 N+/LYS
N31 N+/LYS
N32 N+/LYS
N33 N+/LYS
NAR HIS
OCO COO
OP  TYR/SER?
SH  CYS

Non-titratable:
NP1 AMD?
OH  ROH
O3  ?
CL
F
NAM
N1
O2
""")

    def print_interaction_parameters_latex(self):
        """Print interaction parameters in LaTeX format."""
        # TODO - if these lists and dictionaries are the same as above, then
        # should be constants at the level of the module
        agroups = ['COO', 'HIS', 'CYS', 'TYR', 'SER', 'N+', 'LYS', 'AMD',
                   'ARG', 'TRP', 'ROH', 'CG', 'C2N', 'N30', 'N31', 'N32',
                   'N33', 'NAR', 'OCO', 'NP1', 'OH', 'O3', 'CL', 'F', 'NAM',
                   'N1', 'O2', 'OP', 'SH']
        lgroups = ['CG', 'C2N', 'N30', 'N31', 'N32', 'N33', 'NAR', 'OCO',
                   'NP1', 'OH', 'O3', 'CL', 'F', 'NAM', 'N1', 'O2', 'OP',
                   'SH']
        lines = [
            "",
            "\\begin{{longtable}}{{lllll}}",
            ("\\caption{{Ligand interaction parameters. For interactions not "
             "listed, the default value of {0:s} is applied.}}").format(
                 self.sidechain_cutoffs.default),
            "\\label{{tab:ligand_interaction_parameters}}\\\\",
            "\\toprule",
            "Group1 & Group2 & Interaction & c1 &c2 \\\\",
            "\\midrule",
            "\\endfirsthead",
            "",
            "\\multicolumn{{5}}{{l}}{\\emph{{continued from the previous "
            "page}}}\\\\",
            "\\toprule",
            "Group1 & Group2 & Interaction & c1 &c2 \\\\",
            "\\midrule",
            "\\endhead",
            "",
            "\\midrule",
            "\\multicolumn{{5}}{{r}}{\\emph{{continued on the next "
            "page}}}\\\\",
            "\\endfoot",
            "",
            "\\bottomrule",
            "\\endlastfoot",
            ""]
        str_ = "\n".join(lines)
        for group1 in agroups:
            for group2 in lgroups:
                if self.interaction_matrix[group1][group2] == '-':
                    continue
                if (self.sidechain_cutoffs.get_value(group1, group2)
                        == self.sidechain_cutoffs.default):
                    continue
                fmt = (
                    "{grp1:>3s} & {grp2:>3s} & {mat:1s} & {val1:4} & "
                    "{val2:4}\\\\ \n")
                str_ += fmt.format(
                    group1, group2,
                    self.interaction_matrix[group1][group2],
                    self.sidechain_cutoffs.get_value(group1, group2)[0],
                    self.sidechain_cutoffs.get_value(group1, group2)[1])
                if group1 == group2:
                    break
        str_ += '  \\end{{longtable}}\n'
        _LOGGER.info(str_)

    def print_interactions_latex(self):
        """Print interactions in LaTeX."""
        # TODO - are these the same lists as above? Convert to module
        # constants.
        agroups = ['COO', 'HIS', 'CYS', 'TYR', 'SER', 'N+', 'LYS', 'AMD',
                   'ARG', 'TRP', 'ROH', 'CG', 'C2N', 'N30', 'N31', 'N32',
                   'N33', 'NAR', 'OCO', 'NP1', 'OH', 'O3', 'CL', 'F', 'NAM',
                   'N1', 'O2', 'OP', 'SH']
        lines = [
            "",
            "\\begin{{longtable}}{{{0:s}}}".format('l'*len(agroups)),
            ("\\caption{{Ligand interaction parameters. For interactions not "
             "listed, the default value of {0:s} is applied.}}").format(
                 str(self.sidechain_cutoffs.default)),
            "\\label{{tab:ligand_interaction_parameters}}\\\\",
            "\\toprule",
            "Group1 & Group2 & Interaction & c1 &c2 \\\\",
            "\\midrule",
            "\\endfirsthead",
            "",
            "\\multicolumn{{5}}{{l}}{\\emph{{continued from the previous "
            "page}}}\\\\",
            "\\toprule",
            "Group1 & Group2 & Interaction & c1 &c2 \\\\",
            "\\midrule",
            "\\endhead",
            "",
            "\\midrule",
            "\\multicolumn{{5}}{{r}}{\\emph{{continued on the next "
            "page}}}\\\\",
            "\\endfoot",
            "",
            "\\bottomrule",
            "\\endlastfoot",
            ""
        ]
        str_ = "\n".join(lines)
        for group1 in agroups:
            for group2 in agroups:
                fmt = (
                    '{g1:>3s} & {g2:>3s} & {mat:1s} & {val1:>4s} & '
                    '{val2:>4s}\\\\ \n'
                )
                str_ += fmt.format(
                    group1, group2, self.interaction_matrix[group1][group2],
                    str(self.sidechain_cutoffs.get_value(group1, group2)[0]),
                    str(self.sidechain_cutoffs.get_value(group1, group2)[1]))
                if group1 == group2:
                    break
        str_ += '  \\end{{longtable}}\n'
        _LOGGER.info(str_)


class InteractionMatrix:
    """Interaction matrix class."""

    def __init__(self, name):
        """Initialize with name of matrix.

        Args:
            name:  name of interaction matrix
        """
        self.name = name
        self.value = None
        self.ordered_keys = []
        self.dictionary = {}

    def add(self, words):
        """Add values to matrix.

        Args:
            words:  values to add
        """
        new_group = words[0]
        self.ordered_keys.append(new_group)
        if new_group not in self.dictionary.keys():
            self.dictionary[new_group] = {}
        for i, group in enumerate(self.ordered_keys):
            if len(words) > i+1:
                try:
                    self.value = float(words[i+1])
                except ValueError:
                    self.value = words[i+1]
                self.dictionary[group][new_group] = self.value
                self.dictionary[new_group][group] = self.value

    def get_value(self, item1, item2):
        """Get specific matrix value.

        Args:
            item1:  matrix row index
            item2:  matrix column index
        Returns:
            matrix value or None
        """
        try:
            return self.dictionary[item1][item2]
        except KeyError:
            return None

    def __getitem__(self, group):
        """Get specific group from matrix.

        Args:
            group:  group to get
        """
        if group not in self.dictionary.keys():
            str_ = '{0:s} not found in interaction matrix {1:s}'.format(
                group, self.name)
            raise KeyError(str_)
        return self.dictionary[group]

    def keys(self):
        """Get keys from matrix.

        Returns:
            dictionary key list
        """
        return self.dictionary.keys()

    def __str__(self):
        str_ = '      '
        for key in self.ordered_keys:
            str_ += '{0:>3s} '.format(key)
        str_ += '\n'
        for key1 in self.ordered_keys:
            str_ += '{0:>3s} '.format(key1)
            for key2 in self.ordered_keys:
                str_ += '{0:>3s} '.format(self[key1][key2])
            str_ += '\n'
        return str_


class PairwiseMatrix:
    """Pairwise interaction matrix class."""

    def __init__(self, name):
        """Initialize pairwise matrix.

        Args:
            name:  name of pairwise interaction
        """
        self.name = name
        self.dictionary = {}
        self.default = [0.0, 0.0]

    def add(self, words):
        """Add information to the matrix.

        TODO - this function unnecessarily bundles arguments into a tuple

        Args:
            words:  tuple with assignment information and value
        """
        # assign the default value
        if len(words) == 3 and words[0] == 'default':
            self.default = [float(words[1]), float(words[2])]
            return
        # assign non-default values
        group1 = words[0]
        group2 = words[1]
        value = [float(words[2]), float(words[3])]
        self.insert(group1, group2, value)
        self.insert(group2, group1, value)

    def insert(self, key1, key2, value):
        """Insert value into matrix.

        Args:
            key1:  first matrix key (row)
            key2:  second matrix key (column)
            value:  value to insert
        """
        if key1 in self.dictionary and key2 in self.dictionary[key1]:
            if key1 != key2:
                str_ = (
                    'Parameter value for {0:s}, {1:s} defined more '
                    'than once'.format(key1, key2))
                _LOGGER.warning(str_)
        if key1 not in self.dictionary:
            self.dictionary[key1] = {}
        self.dictionary[key1][key2] = value

    def get_value(self, item1, item2):
        """Get specified value from matrix.

        Args:
            item1:  row index
            item2:  column index
        Returns:
            matrix value (or default)
        """
        try:
            return self.dictionary[item1][item2]
        except KeyError:
            return self.default

    def __getitem__(self, group):
        """Get item from matrix corresponding to specific group.

        Args:
            group:  group to retrieve
        Returns:
            matrix information
        """
        if group not in self.dictionary.keys():
            str_ = '{0:s} not found in interaction matrix {1:s}'.format(
                group, self.name)
            raise KeyError(str_)
        return self.dictionary[group]

    def keys(self):
        """Get keys from matrix.

        Returns:
            dictionary key list
        """
        return self.dictionary.keys()

    def __str__(self):
        str_ = ''
        for key1 in self.keys():
            for key2 in self[key1].keys():
                str_ += '{0:s} {1:s} {2:s}\n'.format(
                    key1, key2, self[key1][key2])
        return str_
