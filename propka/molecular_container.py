"""
PDB molecular container
=======================

Molecular container for storing all contents of PDB files.
"""
import logging
import os
from typing import Dict, List, Optional, Tuple
from propka.parameters import Parameters

import propka.version
from propka.output import write_pka, print_header, print_result
from propka.conformation_container import ConformationContainer
from propka.lib import make_grid, Options


_LOGGER = logging.getLogger(__name__)


class MolecularContainer:
    """Container for storing molecular contents of PDB files.

    TODO - this class name does not conform to PEP8 but has external use.
    We should deprecate and change eventually.


    .. versionchanged:: 3.4.0
       Removed :meth:`write_propka` and
       :meth:`additional_setup_when_reading_input_file` as reading and writing
       PROPKA input files is no longer supported.
    """

    conformation_names: List[str]
    conformations: Dict[str, ConformationContainer]
    name: Optional[str]
    version: propka.version.Version

    def __init__(self, parameters: Parameters, options: Options) -> None:
        """Initialize molecular container.

        Args:
            parameters:  Parameters() object
            options:  options object
        """
        # printing out header before parsing input
        print_header()
        self.conformation_names = []
        self.conformations = {}
        self.options = options
        self.name = None
        try:
            version_class = getattr(propka.version, parameters.version)
            self.version = version_class(parameters)
        except AttributeError as err:
            print(err)
            errstr = 'Error: Version {0:s} does not exist'.format(
                parameters.version)
            raise Exception(errstr)

    def top_up_conformations(self) -> None:
        """Makes sure that all atoms are present in all conformations."""
        ref_atoms = {
            atom.residue_label: atom
            for name in reversed(self.conformation_names)
            for atom in self.conformations[name].atoms
        }
        for conf in self.conformations.values():
            conf.top_up_from_atoms(ref_atoms.values())

    def find_covalently_coupled_groups(self) -> None:
        """Find covalently coupled groups."""
        for name in self.conformation_names:
            self.conformations[name].find_covalently_coupled_groups()

    def find_non_covalently_coupled_groups(self) -> None:
        """Find non-covalently coupled groups."""
        verbose = self.options.display_coupled_residues
        for name in self.conformation_names:
            self.conformations[name].find_non_covalently_coupled_groups(
                verbose=verbose)

    def extract_groups(self) -> None:
        """Identify the groups needed for pKa calculation."""
        for name in self.conformation_names:
            self.conformations[name].extract_groups()

    def calculate_pka(self) -> None:
        """Calculate pKa values."""
        # calculate for each conformation
        for name in self.conformation_names:
            self.conformations[name].calculate_pka(
                self.version, self.options)
        # find non-covalently coupled groups
        self.find_non_covalently_coupled_groups()
        # find the average of the conformations
        self.average_of_conformations()
        # print out the conformation-average results
        print_result(self, 'AVR', self.version.parameters)

    def average_of_conformations(self) -> None:
        """Generate an average of conformations."""
        parameters = self.conformations[self.conformation_names[0]].parameters
        # make a new configuration to hold the average values
        avr_conformation = ConformationContainer(
            name='average', parameters=parameters, molecular_container=self)
        container = self.conformations[self.conformation_names[0]]
        for group in container.get_groups_for_calculations():
            # new group to hold average values
            avr_group = group.clone()
            # sum up all groups ...
            for name in self.conformation_names:
                group_to_add = self.conformations[name].find_group(group)
                if group_to_add:
                    avr_group += group_to_add
                else:
                    str_ = (
                        'Group {0:s} could not be found in '
                        'conformation {1:s}.'.format(
                            group.atom.residue_label, name))
                    _LOGGER.warning(str_)
            # ... and store the average value
            avr_group = avr_group / len(self.conformation_names)
            avr_conformation.groups.append(avr_group)
        # store information on coupling in the average container
        if len(list(filter(lambda c: c.non_covalently_coupled_groups,
                           self.conformations.values()))):
            avr_conformation.non_covalently_coupled_groups = True
        # store chain info
        avr_conformation.chains = self.conformations[
            self.conformation_names[0]].chains
        self.conformations['AVR'] = avr_conformation

    def write_pka(self, filename=None, reference="neutral",
                  direction="folding", options=None) -> None:
        """Write pKa information to a file.

        Args:
            filename:  file to write to
            reference:  reference state
            direction:  folding vs. unfolding
            options:  options object
        """
        if filename is None:
            filename = os.path.join('{0:s}.pka'.format(self.name))
        # if the display_coupled_residues option is true, write the results out
        # to an alternative pka file
        if self.options.display_coupled_residues:
            filename = os.path.join('{0:s}_alt_state.pka'.format(self.name))
        if (hasattr(self.version.parameters, 'output_file_tag')
                and len(self.version.parameters.output_file_tag) > 0):
            filename = os.path.join(
                '{0:s}_{1:s}.pka'.format(
                    self.name, self.version.parameters.output_file_tag))
        write_pka(
            self, self.version.parameters, filename=filename,
            conformation='AVR', reference=reference)

    def get_folding_profile(self, conformation='AVR', reference="neutral",
                            grid: Tuple[float, float, float] = (0., 14., 0.1)):
        """Get a folding profile.

        Args:
            conformation:  conformation to select
            reference:  reference state
            direction:  folding direction (folding)
            grid:  the grid of pH values [min, max, step_size]
            options:  options object
        Returns:
            TODO - figure out what these are
            1. profile
            2. opt
            3. range_80pct
            4. stability_range
        """
        # calculate stability profile
        profile: List[Tuple[float, float]] = []
        for ph in make_grid(*grid):
            conf = self.conformations[conformation]
            ddg = conf.calculate_folding_energy(ph=ph, reference=reference)
            profile.append((ph, ddg))
        # find optimum
        opt: Tuple[Optional[float], float] = (None, 1e6)
        for point in profile:
            opt = min(opt, point, key=lambda v: v[1])
        # find values within 80 % of optimum
        range_80pct: Tuple[Optional[float], Optional[float]] = (None, None)
        values_within_80pct = [p[0] for p in profile if p[1] < 0.8*opt[1]]
        if len(values_within_80pct) > 0:
            range_80pct = (min(values_within_80pct), max(values_within_80pct))
        # find stability range
        stability_range: Tuple[Optional[float], Optional[float]] = (None, None)
        stable_values = [p[0] for p in profile if p[1] < 0.0]
        if len(stable_values) > 0:
            stability_range = (min(stable_values), max(stable_values))
        return profile, opt, range_80pct, stability_range

    def get_charge_profile(self, conformation: str = 'AVR', grid=[0., 14., .1]):
        """Get charge profile for conformation as function of pH.

        Args:
            conformation:  conformation to test
            grid:  grid of pH values [min, max, step]
        Returns:
            list of charge state values
        """
        charge_profile: List[List[float]] = []
        for ph in make_grid(*grid):
            conf = self.conformations[conformation]
            q_unfolded, q_folded = conf.calculate_charge(
                self.version.parameters, ph=ph)
            charge_profile.append([ph, q_unfolded, q_folded])
        return charge_profile

    def get_pi(self, conformation: str = 'AVR', grid=[0., 14., 1], *,
               precision: float = 1e-4) -> Tuple[float, float]:
        """Get the isoelectric points for folded and unfolded states.

        Args:
            conformation:  conformation to test
            grid:  grid of pH values [min, max, step]
            precision:  Compute pI up to this precision
        Returns:
            1. Folded state PI
            2. Unfolded state PI
        """
        conf = self.conformations[conformation]

        WHICH_UNFOLDED = 0
        WHICH_FOLDED = 1

        def pi(which, pH, min_, max_):
            charge = conf.calculate_charge(
                self.version.parameters, ph=pH)[which]
            if max_ - min_ > precision:
                if charge > 0.0:
                    min_ = pH
                else:
                    max_ = pH
                next_pH = (min_ + max_) / 2
                return pi(which, next_pH, min_, max_)
            return pH

        start = (grid[0] + grid[1]) / 2, grid[0], grid[1]

        return (
            pi(WHICH_FOLDED, *start),
            pi(WHICH_UNFOLDED, *start),
        )
