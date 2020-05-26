"""Holds parameters and settings."""
import pkg_resources
import propka.lib as lib
from propka.lib import info, warning


# names and types of all key words in configuration file
MATRICES = ['interaction_matrix']
PAIR_WISE_MATRICES = ['sidechain_cutoffs']
NUMBER_DICTIONARIES = [
    'VanDerWaalsVolume', 'charge', 'model_pkas', 'ions', 'valence_electrons',
    'custom_model_pkas']
LIST_DICTIONARIES = ['backbone_NH_hydrogen_bond', 'backbone_CO_hydrogen_bond']
STRING_DICTIONARIES = ['protein_group_mapping']
STRING_LISTS = [
    'ignore_residues', 'angular_dependent_sidechain_interactions',
    'acid_list', 'base_list', 'exclude_sidechain_interactions',
    'backbone_reorganisation_list', 'write_out_order']
DISTANCES = ['desolv_cutoff', 'buried_cutoff', 'coulomb_cutoff1',
             'coulomb_cutoff2']
PARAMETERS = [
    'Nmin', 'Nmax', 'desolvationSurfaceScalingFactor', 'desolvationPrefactor',
    'desolvationAllowance', 'coulomb_diel', 'COO_HIS_exception',
    'OCO_HIS_exception', 'CYS_HIS_exception', 'CYS_CYS_exception',
    'min_ligand_model_pka', 'max_ligand_model_pka',
    'include_H_in_interactions', 'coupling_max_number_of_bonds',
    'min_bond_distance_for_hydrogen_bonds', 'coupling_penalty',
    'shared_determinants', 'common_charge_centre', 'hide_penalised_group',
    'remove_penalised_group', 'max_intrinsic_pka_diff',
    'min_interaction_energy', 'max_free_energy_diff', 'min_swap_pka_shift',
    'min_pka', 'max_pka', 'sidechain_interaction']
STRINGS = ['version', 'output_file_tag', 'ligand_typing', 'pH', 'reference']


class Parameters:
    """PROPKA parameter class."""

    def __init__(self, parameter_file):
        """Initialize parameter class.

        Args:
            parameter_file:  file with parameters
        """
        # TODO - need to define all members explicitly
        self.model_pkas = {}
        self.interaction_matrix = InteractionMatrix("interaction_matrix")
        self.sidechain_cutoffs = None
        # TODO - it would be nice to rename these; they're defined everywhere
        self.COO_HIS_exception = None
        self.OCO_HIS_exception = None
        self.CYS_HIS_exception = None
        self.CYS_CYS_exception = None
        # These functions set up remaining data structures implicitly
        self.set_up_data_structures()
        self.read_parameters(parameter_file)

    def read_parameters(self, file_):
        """Read parameters from file.

        Args:
            file_:  file to read
        """
        # try to locate the parameters file
        try:
            ifile = pkg_resources.resource_filename(__name__, file_)
            input_ = lib.open_file_for_reading(ifile)
        except (IOError, FileNotFoundError, ValueError):
            input_ = lib.open_file_for_reading(file_)
        for line in input_:
            self.parse_line(line)

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
        if len(words) == 3 and words[0] in NUMBER_DICTIONARIES:
            self.parse_to_number_dictionary(words)
        elif len(words) == 2 and words[0] in STRING_LISTS:
            self.parse_to_string_list(words)
        elif len(words) == 2 and words[0] in DISTANCES:
            self.parse_distance(words)
        elif len(words) == 2 and words[0] in PARAMETERS:
            self.parse_parameter(words)
        elif len(words) == 2 and words[0] in STRINGS:
            self.parse_string(words)
        elif len(words) > 2 and words[0] in LIST_DICTIONARIES:
            self.parse_to_list_dictionary(words)
        elif words[0] in MATRICES+PAIR_WISE_MATRICES:
            self.parse_to_matrix(words)
        elif len(words) == 3 and words[0] in STRING_DICTIONARIES:
            self.parse_to_string_dictionary(words)

    def parse_to_number_dictionary(self, words):
        """Parse field to number dictionary.

        Args:
            words:  strings to parse.
        """
        dict_ = getattr(self, words[0])
        key = words[1]
        value = words[2]
        dict_[key] = float(value)

    def parse_to_string_dictionary(self, words):
        """Parse field to string dictionary.

        Args:
            words:  strings to parse
        """
        dict_ = getattr(self, words[0])
        key = words[1]
        value = words[2]
        dict_[key] = value

    def parse_to_list_dictionary(self, words):
        """Parse field to list dictionary.

        Args:
            words:  strings to parse.
        """
        dict_ = getattr(self, words[0])
        key = words[1]
        if not key in dict_:
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

    def parse_distance(self, words):
        """Parse field to distance.

        Args:
            words:  strings to parse
        """
        value = float(words[1])
        setattr(self, words[0], value)
        value_sq = value*value
        setattr(self, "%s_squared" % words[0], value_sq)

    def parse_parameter(self, words):
        """Parse field to parameters.

        Args:
            words:  strings to parse
        """
        value = float(words[1])
        setattr(self, words[0], value)

    def parse_string(self, words):
        """Parse field to strings.

        Args:
            words:  strings to parse
        """
        setattr(self, words[0], words[1])

    def set_up_data_structures(self):
        """Set up internal data structures.

        TODO - it would be better to make these assignments explicit in
        __init__.
        """
        for key_word in (NUMBER_DICTIONARIES + LIST_DICTIONARIES
                         + STRING_DICTIONARIES):
            setattr(self, key_word, {})
        for key_word in STRING_LISTS:
            setattr(self, key_word, [])
        for key_word in STRINGS:
            setattr(self, key_word, "")
        for key_word in MATRICES:
            matrix = InteractionMatrix(key_word)
            setattr(self, key_word, matrix)
        for key_word in PAIR_WISE_MATRICES:
            matrix = PairwiseMatrix(key_word)
            setattr(self, key_word, matrix)

    def print_interaction_parameters(self):
        """Print interaction parameters."""
        info('--------------- Model pKa values ----------------------')
        for k in self.model_pkas:
            info('%3s %8.2f' % (k, self.model_pkas[k]))

        info('')
        info('--------------- Interactions --------------------------')
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
                interaction = '%3s %3s %1s %4s %4s' % (
                    group1, group2, self.interaction_matrix[group1][group2],
                    self.sidechain_cutoffs.get_value(group1, group2)[0],
                    self.sidechain_cutoffs.get_value(group1, group2)[1])
                map_interaction = ''
                if group2 in map_:
                    for val in map_[group2]:
                        map_interaction += '|%3s %3s %1s %4s %4s' % (
                            group1, val, self.interaction_matrix[group1][val],
                            self.sidechain_cutoffs.get_value(group1, val)[0],
                            self.sidechain_cutoffs.get_value(group1, val)[1])
                        if (self.interaction_matrix[group1][val]
                                != self.interaction_matrix[group1][group2]):
                            map_interaction += '* '
                        if (self.sidechain_cutoffs.get_value(group1, val)[0]
                                != self.sidechain_cutoffs.get_value(group1,
                                                                    group2)[0]
                                or self.sidechain_cutoffs.get_value(group1,
                                                                    val)[1]
                                != self.sidechain_cutoffs.get_value(group1,
                                                                    group2)[1]):
                            map_interaction += '! '
                        else:
                            map_interaction += '  '
                    if (len(map_[group2]) == 0
                            and (self.sidechain_cutoffs.get_value(group1,
                                                                  group2)[0]
                                 != 3
                                 or self.sidechain_cutoffs.get_value(group1,
                                                                     group2)[1]
                                 != 4)):
                        map_interaction += '?  '
                info(interaction, map_interaction)
                if group1 == group2:
                    break
            info('-')
        info('--------------- Exceptions ----------------------------')
        info('COO-HIS', self.COO_HIS_exception)
        info('OCO-HIS', self.OCO_HIS_exception)
        info('CYS-HIS', self.CYS_HIS_exception)
        info('CYS-CYS', self.CYS_CYS_exception)

        info('--------------- Mapping -------------------------------')
        info("""
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
        str_ = """
\\begin{longtable}{lllll}
\\caption{Ligand interaction parameters. For interactions not listed, the default value of %s is applied.}
\\label{tab:ligand_interaction_parameters}\\\\

\\toprule
Group1 & Group2 & Interaction & c1 &c2 \\\\
\\midrule
\\endfirsthead

\\multicolumn{5}{l}{\\emph{continued from the previous page}}\\\\
\\toprule
Group1 & Group2 & Interaction & c1 &c2 \\\\
\\midrule
\\endhead

\\midrule
\\multicolumn{5}{r}{\\emph{continued on the next page}}\\\\
\\endfoot

\\bottomrule
\\endlastfoot

""" % (self.sidechain_cutoffs.default)
        for group1 in agroups:
            for group2 in lgroups:
                if self.interaction_matrix[group1][group2] == '-':
                    continue
                if (self.sidechain_cutoffs.get_value(group1, group2)
                        == self.sidechain_cutoffs.default):
                    continue
                str_ += ('%3s & %3s & %1s & %4s & %4s\\\\ \n'
                         % (group1, group2,
                            self.interaction_matrix[group1][group2],
                            self.sidechain_cutoffs.get_value(group1,
                                                             group2)[0],
                            self.sidechain_cutoffs.get_value(group1,
                                                             group2)[1]))
                if group1 == group2:
                    break
        str_ += '  \\end{longtable}\n'
        info(str_)

    def print_interactions_latex(self):
        """Print interactions in LaTeX."""
        # TODO - are these the same lists as above? Convert to module constants.
        agroups = ['COO', 'HIS', 'CYS', 'TYR', 'SER', 'N+', 'LYS', 'AMD',
                   'ARG', 'TRP', 'ROH', 'CG', 'C2N', 'N30', 'N31', 'N32',
                   'N33', 'NAR', 'OCO', 'NP1', 'OH', 'O3', 'CL', 'F', 'NAM',
                   'N1', 'O2', 'OP', 'SH']
        str_ = """
\\begin{longtable}{%s}
\\caption{Ligand interaction parameters. For interactions not listed, the default value of %s is applied.}
\\label{tab:ligand_interaction_parameters}\\\\

\\toprule
Group1 & Group2 & Interaction & c1 &c2 \\\\
\\midrule
\\endfirsthead

\\multicolumn{5}{l}{\\emph{continued from the previous page}}\\\\
\\toprule
Group1 & Group2 & Interaction & c1 &c2 \\\\
\\midrule
\\endhead

\\midrule
\\multicolumn{5}{r}{\\emph{continued on the next page}}\\\\
\\endfoot

\\bottomrule
\\endlastfoot

""" % ('l'*len(agroups), self.sidechain_cutoffs.default)
        for group1 in agroups:
            for group2 in agroups:
                str_ += ('%3s & %3s & %1s & %4s & %4s\\\\ \n'
                         % (group1, group2,
                            self.interaction_matrix[group1][group2],
                            self.sidechain_cutoffs.get_value(
                                group1, group2)[0],
                            self.sidechain_cutoffs.get_value(
                                group1, group2)[1]))
                if group1 == group2:
                    break
        str_ += '  \\end{longtable}\n'
        info(str_)


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
        if not new_group in self.dictionary.keys():
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
            str_ = '%s not found in interaction matrix %s' % (group, self.name)
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
            str_ += '%3s ' % key
        str_ += '\n'
        for key1 in self.ordered_keys:
            str_ += '%3s ' % key1
            for key2 in self.ordered_keys:
                str_ += '%3s ' % self[key1][key2]
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
                str_ = ('Parameter value for %s, %s defined more than once'
                        % (key1, key2))
                warning(str_)
        if not key1 in self.dictionary:
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
            str_ = '%s not found in interaction matrix %s' % (group, self.name)
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
                str_ += '%s %s %s\n' % (key1, key2, self[key1][key2])
        return str_
