"""
Data structures for groups
==========================

Routines and classes for storing groups important to PROPKA calculations.


.. versionchanged:: 3.4.0
   Removed :func:`initialize_atom_group` as reading PROPKA inputs is no longer
   supported.
"""
import logging
import math
from typing import cast, Dict, Iterable, List, NoReturn, Optional

import propka.ligand
from propka.parameters import Parameters
import propka.protonate
from propka.atom import Atom
from propka.ligand_pka_values import LigandPkaValues
from propka.determinant import Determinant


_LOGGER = logging.getLogger(__name__)


# Constants that start with "UNK_" are a mystery to me
UNK_PKA_SCALING = -1.36
PROTONATOR = propka.protonate.Protonate(verbose=False)
#: acids
EXPECTED_ATOMS_ACID_INTERACTIONS = {
    'COO': {'O': 2}, 'HIS': {'H': 2, 'N': 2}, 'CYS': {'S': 1}, 'TYR': {'O': 1},
    'LYS': {'N': 1}, 'ARG': {'H': 5, 'N': 3}, 'ROH': {'O': 1},
    'AMD': {'H': 2, 'N': 1}, 'TRP': {'H': 1, 'N': 1}, 'N+': {'N': 1},
    'C-': {'O': 2}, 'BBN': {'H': 1, 'N': 1}, 'BBC': {'O': 1},
    'NAR': {'H': 1, 'N': 1}, 'NAM': {'H': 1, 'N': 1}, 'F': {'F': 1},
    'Cl': {'Cl': 1}, 'OH': {'H': 1, 'O': 1}, 'OP': {'O': 1}, 'O3': {'O': 1},
    'O2': {'O': 1}, 'SH': {'S': 1}, 'CG': {'H': 5, 'N': 3},
    'C2N': {'H': 4, 'N': 2}, 'OCO': {'O': 2}, 'N30': {'H': 4, 'N': 1},
    'N31': {'H': 3, 'N': 1}, 'N32': {'H': 2, 'N': 1}, 'N33': {'H': 1, 'N': 1},
    'NP1': {'H': 2, 'N': 1}, 'N1': {'N': 1}}
#: bases
EXPECTED_ATOMS_BASE_INTERACTIONS = {
    'COO': {'O': 2}, 'HIS': {'N': 2}, 'CYS': {'S': 1}, 'TYR': {'O': 1},
    'LYS': {'N': 1}, 'ARG': {'N': 3}, 'ROH': {'O': 1}, 'AMD': {'O': 1},
    'TRP': {'N': 1}, 'N+': {'N': 1}, 'C-': {'O': 2}, 'BBN': {'H': 1, 'N': 1},
    'BBC': {'O': 1}, 'NAR': {'H': 1, 'N': 1}, 'NAM': {'H': 1, 'N': 1},
    'F': {'F': 1}, 'Cl': {'Cl': 1}, 'OH': {'H': 1, 'O': 1}, 'OP': {'O': 1},
    'O3': {'O': 1}, 'O2': {'O': 1}, 'SH': {'S': 1}, 'CG': {'N': 3},
    'C2N': {'N': 2}, 'OCO': {'O': 2}, 'N30': {'N': 1}, 'N31': {'N': 1},
    'N32': {'N': 1}, 'N33': {'N': 1}, 'NP1': {'N': 1}, 'N1': {'N': 1}}


class Group:
    """Class for storing groups important to pKa calculations.


    .. versionchanged:: 3.4.0
       Removed :meth:`make_covalently_coupled_line` and
       :meth:`make_non_covalently_coupled_line` as writing PROPKA inputs is no
       longer supported.
    """

    def __init__(self, atom: Atom):
        """Initialize with an atom.

        Args:
            atom:  atom object
        """
        self.atom = atom
        self.type = ''
        atom.group = self
        # set up data structures
        self.determinants: Dict[str, List[Determinant]] = {
            'sidechain': [],
            'backbone': [],
            'coulomb': [],
        }
        self.pka_value = 0.0
        self.model_pka = 0.0
        # Energy associated with volume interactions
        self.energy_volume = 0.0
        # Number of atoms associated with volume interactions
        self.num_volume = 0.0
        # Energy associated with local interactions
        self.energy_local = 0.0
        # Number of atoms associated with local interactions
        self.num_local = 0.0
        self.buried = 0.0
        self.x = 0.0
        self.y = 0.0
        self.z = 0.0
        self.charge = 0
        self.parameters: Optional[Parameters] = None
        self.exclude_cys_from_results = False
        self.interaction_atoms_for_acids: List[Atom] = []
        self.interaction_atoms_for_bases: List[Atom] = []
        self.model_pka_set = False
        self.intrinsic_pka = None
        self.titratable = False
        # information on covalent and non-covalent coupling
        self.non_covalently_coupled_groups: List["Group"] = []
        self.covalently_coupled_groups: List["Group"] = []
        self.coupled_titrating_group: Optional["Group"] = None
        self.common_charge_centre = False
        self.residue_type = self.atom.res_name
        if self.atom.terminal:
            self.residue_type = self.atom.terminal
        if self.atom.type == 'atom':
            fmt = "{g.residue_type:<3s}{a.res_num:>4d}{a.chain_id:>2s}"
            self.label = fmt.format(g=self, a=atom)
        elif self.atom.res_name in ['DA ', 'DC ', 'DG ', 'DT ']:
            fmt = "{type:1s}{elem:1s}{name:1s}{res_num:>4d}{chain:>2s}"
            self.label = fmt.format(
                type=self.residue_type[1], elem=atom.element,
                name=atom.name.replace('\'', '')[-1], res_num=atom.res_num,
                chain=atom.chain_id)
        else:
            fmt = "{type:<3s}{name:>4s}{chain:>2s}"
            self.label = fmt.format(
                type=self.residue_type, name=atom.name, chain=atom.chain_id)
        # container for squared distances
        self.squared_distances: NoReturn = cast(NoReturn, {})  # FIXME unused?

    def couple_covalently(self, other: "Group") -> None:
        """Couple this group with another group.

        Args:
            other:  other group for coupling
        """
        # do the coupling
        if other not in self.covalently_coupled_groups:
            self.covalently_coupled_groups.append(other)
        if self not in other.covalently_coupled_groups:
            other.covalently_coupled_groups.append(self)

    def couple_non_covalently(self, other: "Group") -> None:
        """Non-covalenthly couple this group with another group.

        Args:
            other:  other group for coupling
        """
        # do the coupling
        if other not in self.non_covalently_coupled_groups:
            self.non_covalently_coupled_groups.append(other)
        if self not in other.non_covalently_coupled_groups:
            other.non_covalently_coupled_groups.append(self)

    def get_covalently_coupled_groups(self):
        """Get covalently coupled groups.

        Returns:
            list of covalently coupled groups.
        """
        return self.covalently_coupled_groups

    def get_non_covalently_coupled_groups(self):
        """Get non-covalently coupled groups.

        Returns:
            list of covalently coupled groups.
        """
        return self.non_covalently_coupled_groups

    def share_determinants(self, others: Iterable["Group"]) -> None:
        """Share determinants between this group and others.

        Args:
            others:  list of other groups
        """
        raise NotImplementedError("unused")
        # for each determinant type
        for other in others:
            if other == self:
                the_other = other
                continue
            for type_ in ['sidechain', 'backbone', 'coulomb']:
                for det in other.determinants[type_]:
                    self.share_determinant(det, type_)
        # recalculate pka values
        self.calculate_total_pka()
        the_other.calculate_total_pka()

    def share_determinant(self, new_determinant: Determinant, type_: str) -> None:
        """Add determinant to this group's list of determinants.

        Args:
            new_determinant:  determinant to add
            type_:  type of determinant
        """
        added = False
        # first check if we already have a determinant with this label
        for own_determinant in self.determinants[type_]:
            if own_determinant.group == new_determinant.group:
                # if so, find the average value
                avr = 0.5*(own_determinant.value + new_determinant.value)
                own_determinant.value = avr
                new_determinant.value = avr
                added = True
        # otherwise we just add the determinant to our list
        if not added:
            self.determinants[type_].append(
                Determinant(new_determinant.group, new_determinant.value))

    def __eq__(self, other):
        """Needed for creating sets of groups."""
        if self.atom.type == 'atom':
            # In case of protein atoms we trust the labels
            return self.label == other.label
        else:
            # For heterogene atoms we also need to check the residue number
            return (
                (self.label == other.label)
                and (self.atom.res_num == other.atom.res_num))

    def __hash__(self):
        """Needed for creating sets of groups."""
        return id(self)

    def __iadd__(self, other):
        if self.type != other.type:
            str_ = (
                'Cannot add groups of different types '
                '({0:s} and {1:s})'.format(self.type, other.type))
            raise Exception(str_)
        # add all values
        self.pka_value += other.pka_value
        self.num_volume += other.num_volume
        self.energy_volume += other.energy_volume
        self.num_local += other.num_local
        self.energy_local += other.energy_local
        self.buried += other.buried
        # and add all determinants
        # TODO - list ['sidechain', 'backbone', 'coulomb'] should be constant
        # This list appears all over the code and should be moved to a constant
        # higher in the package
        for type_ in ['sidechain', 'backbone', 'coulomb']:
            for determinant in other.determinants[type_]:
                self.add_determinant(determinant, type_)
        return self

    def add_determinant(self, new_determinant: Determinant, type_: str) -> None:
        """Add to current and creates non-present determinants.

        Args:
            new_determinant:  new determinant to add
            type_:  determinant type
        """
        # first check if we already have a determinant with this label
        for own_determinant in self.determinants[type_]:
            if own_determinant.group == new_determinant.group:
                # if so, add the value
                own_determinant.value += new_determinant.value
                return
        # otherwise we just add the determinant to our list
        self.determinants[type_].append(Determinant(new_determinant.group,
                                                    new_determinant.value))

    def set_determinant(self, new_determinant: Determinant, type_: str) -> None:
        """Overwrite current and create non-present determinants.

        Args:
            new_determinant:  new determinant to add
            type_:  determinant type
        """
        # first check if we already have a determinant with this label
        for own_determinant in self.determinants[type_]:
            if own_determinant.group == new_determinant.group:
                # if so, overwrite the value
                own_determinant.value = new_determinant.value
                return
        # otherwise we just add the determinant to our list
        self.determinants[type_].append(Determinant(new_determinant.group,
                                                    new_determinant.value))

    def remove_determinants(self, labels):
        """Remove all determinants with specified labels.

        Args:
            labels:  list of labels to remove
        """
        for type_ in ['sidechain', 'backbone', 'coulomb']:
            matches = list(
                filter(lambda d: d.label
                       in labels, [d for d in self.determinants[type_]]))
            for match in matches:
                self.determinants[type_].remove(match)

    def __truediv__(self, value):
        value = float(value)
        # divide all values
        self.pka_value /= value
        self.num_volume /= value
        self.energy_volume /= value
        self.num_local /= value
        self.energy_local /= value
        self.buried /= value
        # and all determinants
        for type_ in ['sidechain', 'backbone', 'coulomb']:
            for determinant in self.determinants[type_]:
                determinant.value /= value
        return self

    def clone(self):
        """Create a copy of this group.

        Returns:
            Copy of this group.
        """
        res = Group(self.atom)
        res.type = self.type
        res.residue_type = self.residue_type
        res.model_pka = self.model_pka
        res.coupled_titrating_group = self.coupled_titrating_group
        res.covalently_coupled_groups = self.covalently_coupled_groups
        res.non_covalently_coupled_groups = self.non_covalently_coupled_groups
        res.titratable = self.titratable
        res.exclude_cys_from_results = self.exclude_cys_from_results
        res.charge = self.charge
        return res

    def setup(self):
        """Set up a group."""
        assert self.parameters is not None
        # set the charges
        if self.type in self.parameters.charge.keys():
            self.charge = self.parameters.charge[self.type]
        if self.residue_type in self.parameters.ions.keys():
            self.charge = self.parameters.ions[self.residue_type]
        # find the center and the interaction atoms
        self.setup_atoms()
        # set the model pka value
        self.titratable = False
        if self.residue_type in self.parameters.model_pkas.keys():
            if not self.model_pka_set:
                self.model_pka = self.parameters.model_pkas[self.residue_type]
                # check if we should apply a custom model pka
                key = '{0:s}-{1:s}'.format(
                    self.atom.res_name.strip(),
                    self.atom.name.strip())
                if key in self.parameters.custom_model_pkas.keys():
                    self.model_pka = self.parameters.custom_model_pkas[key]
                self.model_pka_set = True
        if self.model_pka_set and not self.atom.cysteine_bridge:
            self.titratable = True
        self.exclude_cys_from_results = False

    def setup_atoms(self):
        """Set up atoms in group.

        This method is overwritten for some types of groups
        """
        # set the center at the position of the main atom
        self.set_center([self.atom])
        # set the main atom as interaction atom
        self.set_interaction_atoms([self.atom], [self.atom])

    def set_interaction_atoms(self, interaction_atoms_for_acids: List[Atom],
                              interaction_atoms_for_bases: List[Atom]):
        """Set interacting atoms and group types.

        Args:
            interaction_atoms_for_acids:  list of atoms for acid interactions
            interaction_atoms_for_base:  list of atoms for base interactions
        """
        for atom in interaction_atoms_for_acids + interaction_atoms_for_bases:
            atom.set_group_type(self.type)
        self.interaction_atoms_for_acids = interaction_atoms_for_acids
        self.interaction_atoms_for_bases = interaction_atoms_for_bases
        # check if all atoms have been identified
        ok = True
        for (expect, found) in [
            (EXPECTED_ATOMS_ACID_INTERACTIONS, self.interaction_atoms_for_acids),
            (EXPECTED_ATOMS_BASE_INTERACTIONS, self.interaction_atoms_for_bases),
        ]:
            if self.type in expect.keys():
                for elem in expect[self.type].keys():
                    if (len([a for a in found if a.element == elem])
                            != expect[self.type][elem]):
                        ok = False
        if not ok:
            str_ = 'Missing atoms or failed protonation for '
            str_ += '{0:s} ({1:s}) -- please check the structure'.format(
                self.label, self.type)
            _LOGGER.warning(str_)
            _LOGGER.warning('{0:s}'.format(str(self)))
            num_acid = sum(
                [EXPECTED_ATOMS_ACID_INTERACTIONS[self.type][e]
                 for e in EXPECTED_ATOMS_ACID_INTERACTIONS[self.type].keys()])
            num_base = sum(
                [EXPECTED_ATOMS_BASE_INTERACTIONS[self.type][e]
                 for e in EXPECTED_ATOMS_BASE_INTERACTIONS[self.type].keys()])
            _LOGGER.warning(
                'Expected {0:d} interaction atoms for acids, found:'.format(
                    num_acid))
            for i in range(len(self.interaction_atoms_for_acids)):
                _LOGGER.warning(
                    '             {0:s}'.format(
                        str(self.interaction_atoms_for_acids[i])))
            _LOGGER.warning(
                'Expected {0:d} interaction atoms for bases, found:'.format(
                    num_base))
            for i in range(len(self.interaction_atoms_for_bases)):
                _LOGGER.warning(
                    '             {0:s}'.format(
                        str(self.interaction_atoms_for_bases[i])))

    def get_interaction_atoms(self, interacting_group: "Group") -> List[Atom]:
        """Get atoms involved in interaction with other group.

        Args:
            interacting_group:  other group
        Returns:
            list of atoms
        """
        assert self.parameters is not None
        if interacting_group.residue_type in self.parameters.base_list:
            return self.interaction_atoms_for_bases
        else:
            # default is acid interaction atoms - cf. 3.0
            return self.interaction_atoms_for_acids

    def set_center(self, atoms):
        """Set center of group based on atoms.

        Args:
            atoms:  list of atoms
        """
        if not atoms:
            raise ValueError("At least one atom must be specified")
        # reset center
        self.x = 0.0
        self.y = 0.0
        self.z = 0.0
        # find the average position of atoms
        for atom in atoms:
            self.x += atom.x
            self.y += atom.y
            self.z += atom.z
        self.x /= float(len(atoms))
        self.y /= float(len(atoms))
        self.z /= float(len(atoms))

    def get_determinant_string(self, remove_penalised_group=False):
        """Create a string to identify this determinant.

        Args:
            remove_penalised_group:  Boolean flag to remove penalized groups
        Returns:
            string
        """
        if self.coupled_titrating_group and remove_penalised_group:
            return ''
        number_of_sidechain = len(self.determinants['sidechain'])
        number_of_backbone = len(self.determinants['backbone'])
        number_of_coulomb = len(self.determinants['coulomb'])
        number_of_lines = max(1, number_of_sidechain, number_of_backbone,
                              number_of_coulomb)
        str_ = ""
        for line_number in range(number_of_lines):
            str_ += "{0:s}".format(self.label)
            if line_number == 0:
                str_ += " {0:6.2f}".format(self.pka_value)
                if len(self.non_covalently_coupled_groups) > 0:
                    str_ += '*'
                else:
                    str_ += ' '
                str_ += " {0:4d}{1:>2s} ".format(int(100.0*self.buried), "%")
                str_ += " {0:6.2f} {1:4d}".format(
                    self.energy_volume, int(self.num_volume))
                str_ += " {0:6.2f} {1:4d}".format(
                    self.energy_local, int(self.num_local))
            else:
                str_ += "{0:>40s}".format(" ")
            # add the determinants
            for type_ in ['sidechain', 'backbone', 'coulomb']:
                str_ += self.get_determinant_for_string(type_, line_number)
            # adding end-of-line
            str_ += "\n"
        str_ += "\n"
        return str_

    def get_determinant_for_string(self, type_, number):
        """Return a string describing determinant.

        Args:
            type_:  determinant type
            number:  determinant index number
        Returns:
            string
        """
        if number >= len(self.determinants[type_]):
            return "    0.00 XXX   0 X"
        else:
            determinant = self.determinants[type_][number]
            return "{0:8.2f} {1:s}".format(
                determinant.value, determinant.label)

    def calculate_total_pka(self):
        """Calculate total pKa based on determinants associated with this
        group."""
        # if this is a cysteine involved in a di-sulphide bond
        if self.atom.cysteine_bridge:
            self.pka_value = 99.99
            return
        self.pka_value = (
            self.model_pka + self.energy_volume + self.energy_local)
        for determinant_type in ['sidechain', 'backbone', 'coulomb']:
            for determinant in self.determinants[determinant_type]:
                self.pka_value += determinant.value

    def calculate_intrinsic_pka(self):
        """Calculate the intrinsic pKa values from the desolvation
        determinants, back-bone hydrogen bonds, and side-chain hydrogen bonds
        to non-titratable residues.
        """
        back_bone = 0.0
        for determinant in self.determinants['backbone']:
            value = determinant.value
            back_bone += value
        side_chain = 0.0
        for determinant in self.determinants['sidechain']:
            if determinant.label[0:3] not in [
                    'ASP', 'GLU', 'LYS', 'ARG', 'HIS', 'CYS', 'TYR', 'C- ',
                    'N+ ']:
                value = determinant.value
                side_chain += value
        self.intrinsic_pka = (
            self.model_pka + self.energy_volume + self.energy_local
            + back_bone + side_chain)

    def get_summary_string(self, remove_penalised_group: bool = False) -> str:
        """Create summary string for this group.

        Args:
            remove_penalised_group:  Boolean to ignore penalized groups
        Returns:
            string
        """
        if self.coupled_titrating_group and remove_penalised_group:
            return ''
        ligand_type = ''
        if self.atom.type == 'hetatm':
            ligand_type = self.type
        penalty = ''
        if self.coupled_titrating_group:
            penalty = (
                ' NB: Discarded due to coupling with {0:s}'.format(
                    self.coupled_titrating_group.label))
        fmt = (
            "   {g.label:>9s} {g.pka_value:8.2f} {g.model_pka:10.2f} "
            "{type:>18s}   {penalty:s}\n")
        return fmt.format(g=self, type=ligand_type, penalty=penalty)

    def __str__(self):
        str_ = 'Group ({0:s}) for {1:s}'.format(self.type, str(self.atom))
        return str_

    def calculate_folding_energy(self, parameters, ph=None, reference=None):
        """Return the electrostatic energy of this residue at specified pH.

        Args:
            parameters:  parameters for energy calculation
            ph:  pH value for calculation
            reference:  reference state for calculation
        Returns:
            float describing energy
        """
        if ph is None:
            ph = parameters.pH
        if reference is None:
            reference = parameters.reference
        # If not titratable, the contribution is zero
        if not self.titratable:
            return 0.00
        # calculating the ddg(neutral --> low-pH) contribution
        ddg_neutral = 0.00
        if reference == 'neutral' and self.charge > 0.00:
            pka_prime = self.pka_value
            for determinant in self.determinants['coulomb']:
                if determinant.value > 0.00:
                    pka_prime -= determinant.value
            ddg_neutral = UNK_PKA_SCALING*(pka_prime - self.model_pka)
        # calculating the ddg(low-pH --> pH) contribution
        # folded
        dpka = ph - self.pka_value
        conc_ratio = 10**dpka
        q_pro = math.log10(1+conc_ratio)
        # unfolded
        dpka = ph - self.model_pka
        conc_ratio = 10**dpka
        q_mod = math.log10(1+conc_ratio)
        ddg_low = UNK_PKA_SCALING*(q_pro - q_mod)
        ddg = ddg_neutral + ddg_low
        return ddg

    def calculate_charge(self, _, ph: float = 7.0, state: str = 'folded') -> float:
        """Calculate the charge of the specified state at the specified pH.

        Args:
            _:  parameters for calculation
            ph:  pH value
            state:  "folded" or "unfolded"
        Returns:
            float with charge
        """
        if state == "unfolded":
            q_dpka = self.charge * (self.model_pka - ph)
        else:
            q_dpka = self.charge * (self.pka_value - ph)
        conc_ratio = 10**q_dpka
        charge = self.charge*(conc_ratio/(1.0+conc_ratio))
        return charge

    def use_in_calculations(self) -> bool:
        """Indicate whether group should be included in results report.

        If --titrate_only option is specified, only residues that are
        titratable and are in that list are included; otherwise all titratable
        residues and CYS residues are included.
        """
        return self.titratable or (
            self.residue_type == 'CYS' and not self.exclude_cys_from_results)


class COOGroup(Group):
    """Carboxyl group."""

    def __init__(self, atom):
        Group.__init__(self, atom)
        self.type = 'COO'

    def setup_atoms(self):
        """Set up group."""
        # Identify the two caroxyl oxygen atoms
        the_oxygens = self.atom.get_bonded_elements('O')
        # set the center using the two oxygen carboxyl atoms (if present)
        if the_oxygens:
            self.set_center(the_oxygens)
        else:
            self.set_center([self.atom])
            # TODO - perhaps it would be better to ignore this group completely
            # if the oxygen is missing from this residue?
        self.set_interaction_atoms(the_oxygens, the_oxygens)


class HISGroup(Group):
    """Histidine group."""

    def __init__(self, atom):
        Group.__init__(self, atom)
        self.type = 'HIS'

    def setup_atoms(self):
        """Set up atoms in group."""
        # Find the atoms in the histidine ring
        ring_atoms = propka.ligand.is_ring_member(self.atom)
        if len(ring_atoms) != 5:
            _LOGGER.warning(
                'His group does not seem to contain a ring %s', self
            )
        # protonate ring
        for ring_atom in ring_atoms:
            PROTONATOR.protonate_atom(ring_atom)
        # set the center using the ring atoms
        if ring_atoms:
            self.set_center(ring_atoms)
        else:
            # Missing side-chain atoms
            self.set_center([self.atom])
            # TODO - perhaps it would be better to ignore this group
            # completely?
        # find the hydrogens on the ring-nitrogens
        hydrogens = []
        nitrogens = [ra for ra in ring_atoms if ra.element == 'N']
        for nitrogen in nitrogens:
            hydrogens.extend(nitrogen.get_bonded_elements('H'))
        self.set_interaction_atoms(hydrogens+nitrogens, nitrogens)


class CYSGroup(Group):
    """Cysteine group."""

    def __init__(self, atom):
        Group.__init__(self, atom)
        self.type = 'CYS'


class TYRGroup(Group):
    """Tyrosine group."""

    def __init__(self, atom):
        Group.__init__(self, atom)
        self.type = 'TYR'


class LYSGroup(Group):
    """Lysine group."""

    def __init__(self, atom):
        Group.__init__(self, atom)
        self.type = 'LYS'


class ARGGroup(Group):
    """Arginine group."""

    def __init__(self, atom):
        Group.__init__(self, atom)
        self.type = 'ARG'

    def setup_atoms(self):
        """Set up group."""
        # set the center at the position of the main atom
        self.set_center([self.atom])
        # find all the hydrogens on the nitrogen atoms
        nitrogens = self.atom.get_bonded_elements('N')
        for nitrogen in nitrogens:
            PROTONATOR.protonate_atom(nitrogen)
        hydrogens = []
        for nitrogen in nitrogens:
            hydrogens.extend(nitrogen.get_bonded_elements('H'))
        self.set_interaction_atoms(nitrogens+hydrogens, nitrogens)


class ROHGroup(Group):
    """Alcohol group."""

    def __init__(self, atom):
        Group.__init__(self, atom)
        self.type = 'ROH'


class SERGroup(Group):
    """Serine group."""

    def __init__(self, atom):
        Group.__init__(self, atom)
        self.type = 'SER'


class AMDGroup(Group):
    """Amide group."""

    def __init__(self, atom):
        Group.__init__(self, atom)
        self.type = 'AMD'

    def setup_atoms(self):
        """Setup group."""
        # Identify the oxygen and nitrogen amide atoms
        the_oxygen = self.atom.get_bonded_elements('O')
        the_nitrogen = self.atom.get_bonded_elements('N')
        if not (the_oxygen and the_nitrogen):
            _LOGGER.warning(f"Missing N or O atom: {self}")
            self.set_center([self.atom])
            return
        # add protons to the nitrogen
        PROTONATOR.protonate_atom(the_nitrogen[0])
        the_hydrogens = the_nitrogen[0].get_bonded_elements('H')
        # set the center using the oxygen and nitrogen amide atoms
        self.set_center(the_oxygen+the_nitrogen)
        self.set_interaction_atoms(the_nitrogen + the_hydrogens, the_oxygen)


class TRPGroup(Group):
    """Tryptophan group."""

    def __init__(self, atom):
        Group.__init__(self, atom)
        self.type = 'TRP'

    def setup_atoms(self):
        """Set up atoms in group."""
        # set the center at the position of the main atom
        self.set_center([self.atom])
        # find the hydrogen on the nitrogen atom
        PROTONATOR.protonate_atom(self.atom)
        the_hydrogen = self.atom.get_bonded_elements('H')
        self.set_interaction_atoms(the_hydrogen+[self.atom], [self.atom])


class NtermGroup(Group):
    """N-terminus group."""
    def __init__(self, atom):
        Group.__init__(self, atom)
        self.type = 'N+'


class CtermGroup(Group):
    """C-terminus group."""
    def __init__(self, atom):
        Group.__init__(self, atom)
        # this is to deal with the COO-C- parameter unification.
        self.type = 'COO'

    def setup_atoms(self):
        """Set up atoms in group."""
        # Identify the carbon and other oxygen carboxyl atoms
        the_carbons = self.atom.get_bonded_elements('C')
        if not the_carbons:
            self.set_center([self.atom])
            # TODO - perhaps it would be better to ignore this group completely
            # if the carbon is missing from this residue?
        else:
            the_other_oxygen = the_carbons[0].get_bonded_elements('O')
            the_other_oxygen.remove(self.atom)
            # set the center and interaction atoms
            the_oxygens = [self.atom] + the_other_oxygen
            self.set_center(the_oxygens)
            self.set_interaction_atoms(the_oxygens, the_oxygens)


class BBNGroup(Group):
    """Backbone nitrogen group."""

    def __init__(self, atom):
        Group.__init__(self, atom)
        self.type = 'BBN'
        self.residue_type = 'BBN'

    def setup_atoms(self):
        """Set up atoms in group."""
        # Identify the hydrogen
        PROTONATOR.protonate_atom(self.atom)
        the_hydrogen = self.atom.get_bonded_elements('H')
        # set the center using the nitrogen
        self.set_center([self.atom])
        self.set_interaction_atoms(
            the_hydrogen+[self.atom], the_hydrogen+[self.atom])


class BBCGroup(Group):
    """Backbone carbon group."""

    def __init__(self, atom):
        Group.__init__(self, atom)
        self.type = 'BBC'
        self.residue_type = 'BBC'

    def setup_atoms(self):
        """Set up atoms in group."""
        # Identify the oxygen
        the_oxygen = self.atom.get_bonded_elements('O')
        # set the center using the nitrogen
        self.set_center([self.atom])
        self.set_interaction_atoms(the_oxygen, the_oxygen)


class NARGroup(Group):
    """Unknown group.

    TODO - identify this group.
    """

    def __init__(self, atom):
        Group.__init__(self, atom)
        self.type = 'NAR'
        self.residue_type = 'NAR'
        _LOGGER.info('Found NAR group: %s', atom)

    def setup_atoms(self):
        """Set up atoms in group."""
        # Identify the hydrogen
        PROTONATOR.protonate_atom(self.atom)
        the_hydrogen = self.atom.get_bonded_elements('H')
        # set the center using the nitrogen
        self.set_center([self.atom])
        self.set_interaction_atoms(
            the_hydrogen+[self.atom], the_hydrogen+[self.atom])


class NAMGroup(Group):
    """Unknown group.

    TODO - identify this group.
    """

    def __init__(self, atom):
        Group.__init__(self, atom)
        self.type = 'NAM'
        self.residue_type = 'NAM'
        _LOGGER.info('Found NAM group: %s', atom)

    def setup_atoms(self):
        """Set up atoms in this group."""
        # Identify the hydrogen
        PROTONATOR.protonate_atom(self.atom)
        the_hydrogen = self.atom.get_bonded_elements('H')
        # set the center using the nitrogen
        self.set_center([self.atom])
        self.set_interaction_atoms(
            the_hydrogen+[self.atom], the_hydrogen+[self.atom])


class FGroup(Group):
    """Fluoride group."""

    def __init__(self, atom):
        Group.__init__(self, atom)
        self.type = 'F'
        self.residue_type = 'F'
        _LOGGER.info('Found F   group: %s', atom)


class ClGroup(Group):
    """Chloride group."""

    def __init__(self, atom):
        Group.__init__(self, atom)
        self.type = 'Cl'
        self.residue_type = 'Cl'
        _LOGGER.info('Found Cl   group: %s', atom)


class OHGroup(Group):
    """Hydroxide group."""

    def __init__(self, atom):
        Group.__init__(self, atom)
        self.type = 'OH'
        self.residue_type = 'OH'
        _LOGGER.info('Found OH  group: %s', atom)

    def setup_atoms(self):
        """Set up atoms in this group."""
        # Identify the hydrogen
        PROTONATOR.protonate_atom(self.atom)
        the_hydrogen = self.atom.get_bonded_elements('H')
        # set the center using the nitrogen
        self.set_center([self.atom])
        self.set_interaction_atoms(
            the_hydrogen+[self.atom], the_hydrogen+[self.atom])


class OPGroup(Group):
    """Phosphate group."""

    def __init__(self, atom):
        Group.__init__(self, atom)
        self.type = 'OP'
        self.residue_type = 'OP'
        _LOGGER.info('Found OP  group: %s', atom)

    def setup_atoms(self):
        """Set up atoms in this group."""
        # Identify the hydrogen
        PROTONATOR.protonate_atom(self.atom)
        # set the center using the oxygen
        self.set_center([self.atom])
        self.set_interaction_atoms([self.atom], [self.atom])


class O3Group(Group):
    """Unknown group.

    TODO - identify this group.
    """

    def __init__(self, atom):
        Group.__init__(self, atom)
        self.type = 'O3'
        self.residue_type = 'O3'
        _LOGGER.info('Found O3  group: %s', atom)


class O2Group(Group):
    """Unknown group.

    TODO - identify this group.
    """

    def __init__(self, atom):
        Group.__init__(self, atom)
        self.type = 'O2'
        self.residue_type = 'O2'
        _LOGGER.info('Found O2  group: %s', atom)


class SHGroup(Group):
    """Sulfhydryl group."""

    def __init__(self, atom):
        Group.__init__(self, atom)
        self.type = 'SH'
        self.residue_type = 'SH'
        _LOGGER.info('Found SH  group: %s', atom)


class CGGroup(Group):
    """Guadinium group."""

    def __init__(self, atom):
        Group.__init__(self, atom)
        self.type = 'CG'
        self.residue_type = 'CG'
        _LOGGER.info('Found CG  group: %s', atom)

    def setup_atoms(self):
        """Set up atoms in this group."""
        # Identify the nitrogens
        the_nitrogens = self.atom.get_bonded_elements('N')
        # set the center using the nitrogen
        self.set_center([self.atom])
        the_hydrogens = []
        for nitrogen in the_nitrogens:
            PROTONATOR.protonate_atom(nitrogen)
            the_hydrogens += nitrogen.get_bonded_elements('H')
        self.set_interaction_atoms(
            the_hydrogens+the_nitrogens, the_nitrogens)


class C2NGroup(Group):
    """Amidinium group."""

    def __init__(self, atom):
        Group.__init__(self, atom)
        self.type = 'C2N'
        self.residue_type = 'C2N'
        _LOGGER.info('Found C2N group: %s', atom)

    def setup_atoms(self):
        """Set up atoms in this group."""
        # Identify the nitrogens
        the_nitrogens = self.atom.get_bonded_elements('N')
        the_nitrogens = [
            n for n in the_nitrogens if len(n.get_bonded_heavy_atoms()) == 1]
        # set the center using the nitrogen
        self.set_center([self.atom])
        the_hydrogens = []
        for nitrogen in the_nitrogens:
            PROTONATOR.protonate_atom(nitrogen)
            the_hydrogens += nitrogen.get_bonded_elements('H')
        self.set_interaction_atoms(the_hydrogens+the_nitrogens, the_nitrogens)


class OCOGroup(Group):
    """Carboxyl group."""

    def __init__(self, atom):
        Group.__init__(self, atom)
        self.type = 'OCO'
        self.residue_type = 'OCO'
        _LOGGER.info('Found OCO group: %s', atom)

    def setup_atoms(self):
        """Set up atoms in group."""
        # Identify the two carboxyl oxygen atoms
        the_oxygens = self.atom.get_bonded_elements('O')
        # set the center using the two oxygen carboxyl atoms
        self.set_center(the_oxygens)
        self.set_interaction_atoms(the_oxygens, the_oxygens)


class N30Group(Group):
    """Unknown group.

    TODO - identify this group.
    """

    def __init__(self, atom):
        Group.__init__(self, atom)
        self.type = 'N30'
        self.residue_type = 'N30'
        _LOGGER.info('Found N30 group: %s', atom)

    def setup_atoms(self):
        """Set up atoms in this group."""
        # Identify the nitrogens
        PROTONATOR.protonate_atom(self.atom)
        the_hydrogens = self.atom.get_bonded_elements('H')
        # set the center using the nitrogen
        self.set_center([self.atom])
        self.set_interaction_atoms(the_hydrogens+[self.atom], [self.atom])


class N31Group(Group):
    """Unknown group.

    TODO - identify this group.
    """

    def __init__(self, atom):
        Group.__init__(self, atom)
        self.type = 'N31'
        self.residue_type = 'N31'
        _LOGGER.info('Found N31 group: %s', atom)

    def setup_atoms(self):
        """Set up atoms in this group."""
        # Identify the nitrogens
        PROTONATOR.protonate_atom(self.atom)
        the_hydrogens = self.atom.get_bonded_elements('H')
        # set the center using the nitrogen
        self.set_center([self.atom])
        self.set_interaction_atoms(the_hydrogens+[self.atom], [self.atom])


class N32Group(Group):
    """Unknown group.

    TODO - identify this group.
    """

    def __init__(self, atom):
        Group.__init__(self, atom)
        self.type = 'N32'
        self.residue_type = 'N32'
        _LOGGER.info('Found N32 group: %s', atom)

    def setup_atoms(self):
        """Set up atoms in this group."""
        # Identify the nitrogens
        PROTONATOR.protonate_atom(self.atom)
        the_hydrogens = self.atom.get_bonded_elements('H')
        # set the center using the nitrogen
        self.set_center([self.atom])
        self.set_interaction_atoms(the_hydrogens+[self.atom], [self.atom])


class N33Group(Group):
    """Unknown group.

    TODO - identify this group.
    """

    def __init__(self, atom):
        Group.__init__(self, atom)
        self.type = 'N33'
        self.residue_type = 'N33'
        _LOGGER.info('Found N33 group: %s', atom)

    def setup_atoms(self):
        """Set up atoms in this group."""
        # Identify the nitrogens
        PROTONATOR.protonate_atom(self.atom)
        the_hydrogens = self.atom.get_bonded_elements('H')
        # set the center using the nitrogen
        self.set_center([self.atom])
        self.set_interaction_atoms(the_hydrogens+[self.atom], [self.atom])


class NP1Group(Group):
    """Unknown group.

    TODO - identify this group.
    """

    def __init__(self, atom):
        Group.__init__(self, atom)
        self.type = 'NP1'
        self.residue_type = 'NP1'
        _LOGGER.info('Found NP1 group: %s', atom)

    def setup_atoms(self):
        """Set up atoms in group."""
        # Identify the nitrogens
        PROTONATOR.protonate_atom(self.atom)
        the_hydrogens = self.atom.get_bonded_elements('H')
        # set the center using the nitrogen
        self.set_center([self.atom])
        self.set_interaction_atoms(the_hydrogens+[self.atom], [self.atom])


class N1Group(Group):
    """Unknown group.

    TODO - identify this group.
    """

    def __init__(self, atom):
        Group.__init__(self, atom)
        self.type = 'N1'
        self.residue_type = 'N1'
        _LOGGER.info('Found N1 group: %s', atom)


class IonGroup(Group):
    """Ion group."""

    def __init__(self, atom):
        Group.__init__(self, atom)
        self.type = 'ION'
        self.residue_type = atom.res_name.strip()
        _LOGGER.info('Found ion group: %s', atom)


class NonTitratableLigandGroup(Group):
    """Non-titratable ligand group."""

    def __init__(self, atom):
        Group.__init__(self, atom)
        self.type = 'LG'
        self.residue_type = 'LG'


class TitratableLigandGroup(Group):
    """Titratable ligand group."""

    def __init__(self, atom):
        Group.__init__(self, atom)
        # set the charge and determine type (acid or base)
        self.charge = atom.charge
        if self.charge < 0:
            self.type = 'ALG'
            self.residue_type = 'ALG'
        elif self.charge > 0:
            self.type = 'BLG'
            self.residue_type = 'BLG'
        else:
            raise Exception('Unable to determine type of ligand group - '
                            'charge not set?')
        # check if marvin model pka has been calculated
        # this is not true if we are reading an input file
        if atom.marvin_pka:
            self.model_pka = atom.marvin_pka
            _LOGGER.info(
                'Titratable ligand group     %s %s %s',
                atom, self.model_pka, self.charge
            )
        self.model_pka_set = True


def is_group(parameters: Parameters, atom: Atom) -> Optional[Group]:
    """Identify whether the atom belongs to a group.

    Args:
        parameters:  parameters for check
        atom:  atom to check
    Returns:
        group for atom or None
    """
    atom.groups_extracted = True
    # check if this atom belongs to a protein group
    protein_group = is_protein_group(parameters, atom)
    if protein_group:
        return protein_group
    # check if this atom belongs to a ion group
    ion_group = is_ion_group(parameters, atom)
    if ion_group:
        return ion_group
    # check if this atom belongs to a ligand group
    if parameters.ligand_typing == 'marvin':
        ligand_group = is_ligand_group_by_marvin_pkas(parameters, atom)
    elif parameters.ligand_typing == 'sybyl':
        ligand_group = None
    elif parameters.ligand_typing == 'groups':
        ligand_group = is_ligand_group_by_groups(parameters, atom)
    else:
        raise Exception(
            'Unknown ligand typing method \'{0:s}\''.format(
                parameters.ligand_typing))
    if ligand_group:
        return ligand_group
    return None


def is_protein_group(parameters, atom: Atom) -> Optional[Group]:
    """Identify whether the atom belongs to a protein group.

    Args:
        parameters:  parameters for check
        atom:  atom to check
    Returns:
        group for atom or None
    """
    if atom.type != 'atom':
        return None
    # Check for termial groups
    if atom.terminal == 'N+':
        return NtermGroup(atom)
    elif atom.terminal == 'C-':
        return CtermGroup(atom)
    # Backbone
    if atom.type == 'atom' and atom.name == 'N':
        # ignore proline backbone nitrogens
        if atom.res_name != 'PRO':
            return BBNGroup(atom)
    if atom.type == 'atom' and atom.name == 'C':
        # ignore C- carboxyl
        if atom.count_bonded_elements('O') == 1:
            return BBCGroup(atom)
    # Filters for side chains based on PDB protein atom names
    key = '{0:s}-{1:s}'.format(atom.res_name, atom.name)
    if key in parameters.protein_group_mapping.keys():
        class_str = "{0:s}Group".format(parameters.protein_group_mapping[key])
        group_class = globals()[class_str]
        return group_class(atom)
    return None


def is_ligand_group_by_groups(_, atom: Atom) -> Optional[Group]:
    """Identify whether the atom belongs to a ligand group by checking groups.

    Args:
        _:  parameters for check
        atom:  atom to check
    Returns:
        group for atom or None
    """
    # Ligand group filters
    if atom.type != 'hetatm':
        return None
    PROTONATOR.protonate_atom(atom)
    if atom.sybyl_type == 'N.ar':
        if len(atom.get_bonded_heavy_atoms()) == 2:
            return NARGroup(atom)
    if atom.sybyl_type == 'N.am':
        return NAMGroup(atom)
    if atom.sybyl_type in ['N.3', 'N.4']:
        heavy_bonded = atom.get_bonded_heavy_atoms()
        if len(heavy_bonded) == 0:
            return N30Group(atom)
        elif len(heavy_bonded) == 1:
            return N31Group(atom)
        elif len(heavy_bonded) == 2:
            return N32Group(atom)
        elif len(heavy_bonded) == 3:
            return N33Group(atom)
    if atom.sybyl_type == 'N.1':
        return N1Group(atom)
    if atom.sybyl_type == 'N.pl3':
        # make sure that this atom is not part of a guadinium or amidinium
        # group
        bonded_carbons = atom.get_bonded_elements('C')
        if len(bonded_carbons) == 1:
            bonded_nitrogens = bonded_carbons[0].get_bonded_elements('N')
            if len(bonded_nitrogens) == 1:
                return NP1Group(atom)
    if atom.sybyl_type == 'C.2':
        # Guadinium and amidinium groups
        bonded_nitrogens = atom.get_bonded_elements('N')
        npls = [
            n for n in bonded_nitrogens
            if (n.sybyl_type == 'N.pl3'
                and len(n.get_bonded_heavy_atoms()) == 1)]
        if len(npls) == 2:
            n_with_max_two_heavy_atom_bonds = [
                n for n in bonded_nitrogens
                if len(n.get_bonded_heavy_atoms()) < 3]
            if len(n_with_max_two_heavy_atom_bonds) == 2:
                return C2NGroup(atom)
            if len(bonded_nitrogens) == 3:
                return CGGroup(atom)
        # carboxyl group
        bonded_oxygens = atom.get_bonded_elements('O')
        bonded_oxygens = [b for b in bonded_oxygens if 'O.co2' in b.sybyl_type]
        if len(bonded_oxygens) == 2:
            return OCOGroup(atom)
    if atom.sybyl_type == 'F':
        return FGroup(atom)
    if atom.sybyl_type == 'Cl':
        return ClGroup(atom)
    if atom.sybyl_type == 'O.3':
        if len(atom.get_bonded_heavy_atoms()) == 1:
            # phosphate group
            if atom.count_bonded_elements('P') == 1:
                return OPGroup(atom)
            # hydroxyl group
            else:
                return OHGroup(atom)
        # other SP3 oxygen
        else:
            return O3Group(atom)
    if atom.sybyl_type == 'O.2':
        return O2Group(atom)
    if atom.sybyl_type == 'S.3':
        # sulphide group
        if len(atom.get_bonded_heavy_atoms()) == 1:
            return SHGroup(atom)
    return None


def is_ligand_group_by_marvin_pkas(parameters: Parameters, atom: Atom) -> Optional[Group]:
    """Identify whether the atom belongs to a ligand group by calculating
    'Marvin pKas'.

    Args:
        parameters:  parameters for check
        atom:  atom to check
    Returns:
        group for atom or None
    """
    if atom.type != 'hetatm':
        return None
    # calculate Marvin ligand pkas for this conformation container
    # if not already done
    # TODO - double-check testing coverage of these functions.
    assert atom.molecular_container is not None
    assert atom.conformation_container is not None
    if not atom.conformation_container.marvin_pkas_calculated:
        lpka = LigandPkaValues(parameters)
        lpka.get_marvin_pkas_for_molecular_container(
            atom.molecular_container, min_ph=parameters.min_ligand_model_pka,
            max_ph=parameters.max_ligand_model_pka)
    if atom.marvin_pka:
        return TitratableLigandGroup(atom)
    # Special case of oxygen in carboxyl group not assigned pka value by marvin
    if atom.sybyl_type == 'O.co2':
        atom.charge = -1.0
        other_oxygen = [
            o for o
            in atom.get_bonded_elements('C')[0].get_bonded_elements('O')
            if not o == atom][0]
        atom.marvin_pka = other_oxygen.marvin_pka
        return TitratableLigandGroup(atom)
    raise NotImplementedError("hydrogen_bonds")
    if atom.element in parameters.hydrogen_bonds.elements:
        return NonTitratableLigandGroup(atom)
    return None


def is_ion_group(parameters, atom: Atom) -> Optional[Group]:
    """Identify whether the atom belongs to an ion group.

    Args:
        parameters:  parameters for check
        atom:  atom to check
    Returns:
        group for atom or None
    """
    if atom.res_name.strip() in parameters.ions.keys():
        return IonGroup(atom)
    return None
