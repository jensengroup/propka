"""
Molecular data structures
=========================

Container data structure for molecular conformations.
"""
import logging
import functools
from typing import Callable, Dict, Iterable, Iterator, List, NoReturn, Optional, TYPE_CHECKING, Set

from propka.lib import Options
from propka.version import Version

if TYPE_CHECKING:
    from propka.atom import Atom
    from propka.molecular_container import MolecularContainer

import propka.ligand
from propka.output import make_interaction_map
from propka.determinant import Determinant
from propka.coupled_groups import NCCG
from propka.determinants import set_backbone_determinants, set_ion_determinants
from propka.determinants import set_determinants
from propka.group import Group, is_group
from propka.parameters import Parameters


_LOGGER = logging.getLogger(__name__)

CallableGroupToGroups = Callable[[Group], List[Group]]


#: A large number that gets multipled with the integer obtained from applying
#: :func:`ord` to the atom chain ID.  Used in calculating atom keys for
#: sorting.
UNICODE_MULTIPLIER = 1e7

#: A number that gets mutiplied with an atom's residue number.  Used in
#: calculating keys for atom sorting.
RESIDUE_MULTIPLIER = 1000


class ConformationContainer:
    """Container for molecular conformations


    .. versionchanged:: 3.4.0
       Removed :meth:`additional_setup_when_reading_input_files` as reading
       PROPKA inputs is no longer supported.
    """

    def __init__(self,
                 name: str,
                 parameters: Parameters,
                 molecular_container: "MolecularContainer"):
        """Initialize conformation container.

        Args:
            name:  name for conformation
            parameters:  parmameters for conformation
            molecular_container:  container for molecule
        """
        self.molecular_container = molecular_container
        self.name = name
        self.parameters = parameters
        self.atoms: List["Atom"] = []
        self.groups: List[Group] = []
        self.chains: List[str] = []
        self.current_iter_item = 0
        self.marvin_pkas_calculated = False
        self.non_covalently_coupled_groups = False

    def extract_groups(self):
        """Generate molecular groups needed for calculating pKa values."""
        for atom in self.get_non_hydrogen_atoms():
            # has this atom been checked for groups?
            if atom.groups_extracted == 0:
                group = is_group(self.parameters, atom)
            else:
                group = atom.group
                # if the atom has been checked in a another conformation, check
                # if it has a group that should be used in this conformation
                # as well
            if group:
                self.setup_and_add_group(group)

    def set_common_charge_centres(self):
        """Assign charge centers to groups."""
        for system in self.get_coupled_systems(
                self.get_covalently_coupled_groups(),
                Group.get_covalently_coupled_groups):
            # make a list of the charge centre coordinates
            all_coordinates = list(map(lambda g: [g.x, g.y, g.z], system))
            # find the common charge center
            ccc = functools.reduce(
                lambda g1, g2: [g1[0]+g2[0], g1[1]+g2[1], g1[2]+g2[2]],
                all_coordinates)
            ccc = list(map(lambda c: c/len(system), ccc))
            # set the ccc for all coupled groups in this system
            for group in system:
                [group.x, group.y, group.z] = ccc
                group.common_charge_centre = True

    def find_covalently_coupled_groups(self):
        """Find covalently coupled groups and set common charge centres."""
        for group in self.get_titratable_groups():
            # Find covalently bonded groups
            bonded_groups = self.find_bonded_titratable_groups(
                group.atom, 1, group.atom)
            # coupled groups
            for bond_group in bonded_groups:
                if bond_group in group.covalently_coupled_groups:
                    continue
                if bond_group.atom.sybyl_type == group.atom.sybyl_type:
                    group.couple_covalently(bond_group)
        # check if we should set a common charge centre as well
        if self.parameters.common_charge_centre:
            self.set_common_charge_centres()
        # print coupling map
        map_ = make_interaction_map(
            'Covalent coupling map for {0:s}'.format(str(self)),
            self.get_covalently_coupled_groups(),
            lambda g1, g2: g1 in g2.covalently_coupled_groups)
        _LOGGER.info("Coupling map:\n%s", map_)

    def find_non_covalently_coupled_groups(self, verbose=False):
        """Find non-covalently coupled groups and set common charge centres.

        Args:
            verbose:  verbose output
        """
        # check if non-covalent coupling has already been set up in an input
        # file
        if len(list(filter(lambda g: len(g.non_covalently_coupled_groups) > 0,
                           self.get_titratable_groups()))) > 0:
            self.non_covalently_coupled_groups = True
        NCCG.identify_non_covalently_coupled_groups(self, verbose=verbose)
        # re-do the check
        if len(list(filter(lambda g: len(g.non_covalently_coupled_groups) > 0,
                           self.get_titratable_groups()))) > 0:
            self.non_covalently_coupled_groups = True

    def find_bonded_titratable_groups(self, atom: "Atom", num_bonds: int,
                                      original_atom: "Atom"):
        """Find bonded titrable groups.

        Args:
            atom:  atom to check for bonds
            num_bonds:  number of bonds for coupling
            original_atom:  another atom to check for bonds
        Returns:
            a set of bonded atom groups
        """
        assert self.parameters is not None
        res: Set[Group] = set()
        for bond_atom in atom.bonded_atoms:
            # skip the original atom
            if bond_atom == original_atom:
                continue
            # check if this atom has a titratable group
            if (bond_atom.group and bond_atom.group.titratable
                    and num_bonds
                    <= self.parameters.coupling_max_number_of_bonds):
                res.add(bond_atom.group)
            # check for titratable groups bonded to this atom
            if num_bonds < self.parameters.coupling_max_number_of_bonds:
                res |= self.find_bonded_titratable_groups(
                    bond_atom, num_bonds+1, original_atom)
        return res

    def setup_and_add_group(self, group: Optional[Group]):
        """Check if we want to include this group in the calculations.

        Args:
            group:  group to check
        """
        # Is it recognized as a group at all?
        if not group:
            return
        # Other checks (include ligands, which chains etc.)
        # if all ok, accept the group
        self.init_group(group)
        self.groups.append(group)

    def init_group(self, group: Group):
        """Initialize the given Group object.

        Args:
            group:  group object to initialize
        """
        # set up the group
        group.parameters = self.parameters
        group.setup()

        # If --titrate_only option is set, make non-specified residues
        # un-titratable:
        assert self.molecular_container.options is not None
        titrate_only = self.molecular_container.options.titrate_only
        if titrate_only is not None:
            atom = group.atom
            if (atom.chain_id, atom.res_num, atom.icode) not in titrate_only:
                group.titratable = False
                if group.residue_type == 'CYS':
                    group.exclude_cys_from_results = True

    def calculate_pka(self, version: Version, options: Options):
        """Calculate pKas for conformation container.

        Args:
            version:  version object
            options:  option object
        """
        _LOGGER.info('Calculating pKas for %s', self)
        # calculate desolvation
        for group in self.get_titratable_groups() + self.get_ions():
            version.calculate_desolvation(group)
        # calculate backbone interactions
        set_backbone_determinants(
            self.get_titratable_groups(), self.get_backbone_groups(), version)
        # setting ion determinants
        set_ion_determinants(self, version)
        # calculating the back-bone reorganization/desolvation term
        version.calculate_backbone_reorganization(self)
        # setting remaining non-iterative and iterative side-chain & Coulomb
        # interaction determinants
        set_determinants(
            self.get_sidechain_groups(), version=version, options=options)
        # calculating the total pKa values
        for group in self.groups:
            group.calculate_total_pka()
        # take coupling effects into account
        penalised_labels = self.coupling_effects()
        if (self.parameters.remove_penalised_group
                and len(penalised_labels) > 0):
            _LOGGER.info('Removing penalised groups!!!')
            for group in self.get_titratable_groups():
                group.remove_determinants(penalised_labels)
            # re-calculating the total pKa values
            for group in self.groups:
                group.calculate_total_pka()

    def coupling_effects(self):
        """Penalize groups based on coupling effects.

        Bases: The group with the highest pKa (the most stable one in the
        charged form) will be the first one to adopt a proton as pH is lowered
        and this group is allowed to titrate. The remaining groups are
        penalised.

        Acids: The group with the highest pKa (the least stable one in the
        charged form) will be the last group to loose the proton as pH is
        raised and will be penalised. The remaining groups are allowed to
        titrate.
        """
        penalised_labels = []
        for all_groups in self.get_coupled_systems(
                self.get_covalently_coupled_groups(),
                Group.get_covalently_coupled_groups):
            # check if we should share determinants
            if self.parameters.shared_determinants:
                self.share_determinants(all_groups)
            # find the group that has the highest pKa value
            first_group = max(all_groups, key=lambda g: g.pka_value)
            # In case of acids
            if first_group.charge < 0:
                first_group.coupled_titrating_group = min(
                    all_groups, key=lambda g: g.pka_value)
                # group with the highest pKa is penalised
                penalised_labels.append(first_group.label)
            # In case of bases
            else:
                for group in all_groups:
                    if group == first_group:
                        # group with the highest pKa is allowed to titrate...
                        continue
                    group.coupled_titrating_group = first_group
                    # ... and the rest are penalised
                    penalised_labels.append(group.label)
        return penalised_labels

    @staticmethod
    def share_determinants(groups: Iterable[Group]):
        """Share sidechain, backbone, and Coloumb determinants between groups.

        Args:
            groups:  groups to share between
        """
        # make a list of the determinants to share
        types = ['sidechain', 'backbone', 'coulomb']
        for type_ in types:
            # find maximum value for each determinant
            max_dets: Dict[Group, float] = {}
            for group in groups:
                for det in group.determinants[type_]:
                    # update max dets
                    if det.group not in max_dets.keys():
                        max_dets[det.group] = det.value
                    else:
                        max_dets[det.group] = max(det.value,
                                                  max_dets[det.group],
                                                  key=lambda v: abs(v))
            # overwrite/add maximum value for each determinant
            for det_group in max_dets:
                new_determinant = Determinant(det_group, max_dets[det_group])
                for group in groups:
                    group.set_determinant(new_determinant, type_)

    def get_coupled_systems(
        self,
        groups: Iterable[Group],
        get_coupled_groups: CallableGroupToGroups,
    ) -> Iterator[Set[Group]]:
        """A generator that yields covalently coupled systems.

        Args:
            groups:  groups for generating coupled systems
            get_coupled_groups:  TODO - I don't know what this is
        Yields:
            covalently coupled systems
        """
        groups = set(groups)
        while len(groups) > 0:
            # extract a system of coupled groups ...
            system: Set[Group] = set()
            self.get_a_coupled_system_of_groups(
                groups.pop(), system, get_coupled_groups)
            # ... and remove them from the list
            groups -= system
            yield system

    def get_a_coupled_system_of_groups(self, new_group: Group,
                                       coupled_groups: Set[Group],
                                       get_coupled_groups: CallableGroupToGroups):
        """Set up coupled systems of groups.

        Args:
            new_group:  added to coupled_groups
            coupled_groups:  existing coupled groups
            get_coupled_groups:  TODO - I don't know what this
        """
        coupled_groups.add(new_group)
        for coupled_group in get_coupled_groups(new_group):
            if coupled_group not in coupled_groups:
                self.get_a_coupled_system_of_groups(coupled_group,
                                                    coupled_groups,
                                                    get_coupled_groups)

    def calculate_folding_energy(self, ph=None, reference=None):
        """Calculate folding energy over all groups in conformation container.

        Args:
            ph:  pH for calculation
            reference:  reference state
        Returns:
            folding energy
            TODO - need units
        """
        ddg = 0.0
        for group in self.groups:
            ddg += group.calculate_folding_energy(self.parameters, ph=ph,
                                                  reference=reference)
        return ddg

    def calculate_charge(self, parameters: Parameters, ph: float):
        """Calculate charge for folded and unfolded states.

        Args:
            parameters:  parameters for calculation
            ph:  pH for calculation
        Returns:
            1. charge for unfolded state
            2. charge for folded state
        """
        unfolded = folded = 0.0
        for group in self.get_titratable_groups():
            unfolded += group.calculate_charge(parameters, ph=ph,
                                               state='unfolded')
            folded += group.calculate_charge(parameters, ph=ph,
                                             state='folded')
        return unfolded, folded

    def get_backbone_groups(self) -> List[Group]:
        """Get backbone groups needed for the pKa calculations.

        Returns:
            list of groups
        """
        return [group for group in self.groups if 'BB' in group.type]

    def get_sidechain_groups(self):
        """Get sidechain groups needed for the pKa calculations.

        Returns:
            list of groups
        """
        return [
            group for group in self.groups
            if ('BB' not in group.type and not group.atom.cysteine_bridge)]

    def get_covalently_coupled_groups(self):
        """Get covalently coupled groups needed for pKa calculations.

        Returns:
            list of groups
        """
        return [
            g for g in self.groups
            if len(g.covalently_coupled_groups) > 0]

    def get_non_covalently_coupled_groups(self):
        """Get non-covalently coupled groups needed for pKa calculations.

        Returns:
            list of groups
        """
        return [
            g for g in self.groups
            if len(g.non_covalently_coupled_groups) > 0]

    def get_backbone_nh_groups(self):
        """Get NH backbone groups needed for pKa calculations.

        Returns:
            list of groups
        """
        return [group for group in self.groups if group.type == 'BBN']

    def get_backbone_co_groups(self):
        """Get CO backbone groups needed for pKa calculations.

        Returns:
            list of groups
        """
        return [group for group in self.groups if group.type == 'BBC']

    def get_groups_in_residue(self, residue):
        """Get residue groups needed for pKa calculations.

        Args:
            residue:  specific residue with groups
        Returns:
            list of groups
        """
        return [
            group for group in self.groups if group.residue_type == residue]

    def get_titratable_groups(self):
        """Get all titratable groups needed for pKa calculations.

        Returns:
            list of groups
        """
        return [group for group in self.groups if group.titratable]

    def get_groups_for_calculations(self):
        """Get a list of groups that should be included in results report.

        If --titrate_only option is specified, only residues that are
        titratable and are in that list are included; otherwise all titratable
        residues and CYS residues are included.

        Returns:
            list of groups
        """
        return [group for group in self.groups if group.use_in_calculations()]

    def get_acids(self):
        """Get acid groups needed for pKa calculations.

        Returns:
            list of groups
        """
        return [
            group for group in self.groups
            if (group.residue_type in self.parameters.acid_list
                and not group.atom.cysteine_bridge)]

    def get_backbone_reorganisation_groups(self):
        """Get groups involved with backbone reorganization.

        Returns:
            list of groups
        """
        return [
            group for group in self.groups
            if (group.residue_type
                in self.parameters.backbone_reorganisation_list
                and not group.atom.cysteine_bridge)]

    def get_ions(self):
        """Get ion groups.

        Returns:
            list of groups
        """
        return [
            group for group in self.groups
            if group.residue_type in self.parameters.ions.keys()]

    def get_group_names(self, group_list: NoReturn) -> NoReturn:  # FIXME unused?
        """Get names of groups in list.

        Args:
            group_list:  list to check
        Returns:
            list of groups
        """
        if TYPE_CHECKING:
            assert False
        return [group for group in self.groups if group.type in group_list]

    def get_ligand_atoms(self) -> List["Atom"]:
        """Get atoms associated with ligands.

        Returns:
            list of atoms
        """
        return [atom for atom in self.atoms if atom.type == 'hetatm']

    def get_heavy_ligand_atoms(self) -> List["Atom"]:
        """Get heavy atoms associated with ligands.

        Returns:
            list of atoms
        """
        return [
            atom for atom in self.atoms
            if atom.type == 'hetatm' and atom.element != 'H']

    def get_chain(self, chain: str) -> List["Atom"]:
        """Get atoms associated with a specific chain.

        Args:
            chain:  chain to select
        Returns:
            list of atoms
        """
        return [atom for atom in self.atoms if atom.chain_id != chain]

    def add_atom(self, atom: "Atom"):
        """Add atom to container.

        Args:
            atom:  atom to add
        """
        self.atoms.append(atom)
        if not atom.conformation_container:
            atom.conformation_container = self
        if not atom.molecular_container:
            atom.molecular_container = self.molecular_container
        # store chain id for bookkeeping
        if atom.chain_id not in self.chains:
            self.chains.append(atom.chain_id)

    def copy_atom(self, atom):
        """Add a copy of the atom to container.

        Args:
            atom:  atom to copy and add
        """
        new_atom = atom.make_copy()
        self.atoms.append(new_atom)
        new_atom.conformation_container = self

    def get_non_hydrogen_atoms(self):
        """Get atoms that are not hydrogens.

        Returns:
            list of atoms
        """
        return [atom for atom in self.atoms if atom.element != 'H']

    def top_up(self, other):
        """Adds any atoms found in `other` but not in this container.

        Tops up self with all atoms found in other but not in self.

        Args:
            other:  conformation container with atoms to add
        """
        self.top_up_from_atoms(other.atoms)

    def top_up_from_atoms(self, other_atoms: Iterable["Atom"]):
        """Adds atoms which are missing from this container.

        Args:
            other_atoms: Reference atoms
        """
        my_residue_labels = {a.residue_label for a in self.atoms}
        res_names = {(a.chain_id, a.res_num): a.res_name for a in self.atoms}
        for atom in other_atoms:
            if atom.residue_label not in my_residue_labels:
                if res_names.setdefault((atom.chain_id, atom.res_num),
                                        atom.res_name) != atom.res_name:
                    # don't merge different residue types, e.g. alt-loc mutant
                    continue
                self.copy_atom(atom)

    def find_group(self, group):
        """Find a group in the container.

        Args:
            group:  group to find
        Returns:
            False (if group not found) or group
        """
        for group_ in self.groups:
            if group_.atom.residue_label == group.atom.residue_label:
                if group_.type == group.type:
                    return group_
        return False

    def set_ligand_atom_names(self):
        """Set names for atoms in ligands."""
        for atom in self.get_ligand_atoms():
            propka.ligand.assign_sybyl_type(atom)

    def __str__(self):
        """String that lists statistics of atoms and groups."""
        fmt = (
            "Conformation container {name} with {natoms:d} atoms and "
            "{ngroups:d} groups")
        str_ = fmt.format(
            name=self.name, natoms=len(self), ngroups=len(self.groups))
        return str_

    def __len__(self):
        """Number of atoms in container."""
        return len(self.atoms)

    def sort_atoms(self):
        """Sort atoms by `self.sort_atoms_key()` and renumber."""
        # sort the atoms ...
        self.atoms.sort(key=self.sort_atoms_key)
        # ... and re-number them
        for i in range(len(self.atoms)):
            self.atoms[i].numb = i+1

    @staticmethod
    def sort_atoms_key(atom: "Atom") -> float:
        """Generate key for atom sorting.

        Args:
            atom:  atom for key generation.
        Returns:
            key for atom
        """
        key = ord(atom.chain_id) * UNICODE_MULTIPLIER
        key += atom.res_num * RESIDUE_MULTIPLIER
        if len(atom.name) > len(atom.element):
            key += ord(atom.name[len(atom.element)])
        return key
