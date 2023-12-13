"""
Coupling between groups
=======================

Describe and analyze energetic coupling between groups.
"""
import logging
import itertools
from typing import Optional
import propka.lib
from propka.group import Group
from propka.output import make_interaction_map
from propka.parameters import Parameters


_LOGGER = logging.getLogger(__name__)


class NonCovalentlyCoupledGroups:
    """Groups that are coupled without covalent bonding."""
    parameters: Optional[Parameters] = None
    do_prot_stat = True

    def is_coupled_protonation_state_probability(self, group1, group2,
                                                 energy_method,
                                                 return_on_fail=True):
        """Check whether two groups are energetically coupled.

        Args:
            group1:  first group for interaction
            group2:  second group for interaction
            energy_method:  function for calculating energy
            return_on_fail:  return if part of the calculation fails
        Returns:
            dictionary describing coupling
        """
        assert self.parameters is not None
        # check if the interaction energy is high enough
        interaction_energy = max(self.get_interaction(group1, group2),
                                 self.get_interaction(group2, group1))
        if (interaction_energy <= self.parameters.min_interaction_energy
                and return_on_fail):
            return {'coupling_factor': -1.0}
        # calculate intrinsic pKa's, if not already done
        for group in [group1, group2]:
            if group.intrinsic_pka is None:
                group.calculate_intrinsic_pka()
        use_ph = self.parameters.pH
        if self.parameters.pH == 'variable':
            use_ph = min(group1.pka_value, group2.pka_value)
        default_energy = energy_method(ph=use_ph,
                                       reference=self.parameters.reference)
        default_pka1 = group1.pka_value
        default_pka2 = group2.pka_value
        # check that pka values are within relevant limits
        if (max(default_pka1, default_pka2) < self.parameters.min_pka
                or min(default_pka1, default_pka2) > self.parameters.max_pka):
            if return_on_fail:
                return {'coupling_factor': -1.0}
        # Swap interactions and re-calculate pKa values
        self.swap_interactions([group1], [group2])
        group1.calculate_total_pka()
        group2.calculate_total_pka()
        # store swapped energy and pka's
        swapped_energy = energy_method(
            ph=use_ph, reference=self.parameters.reference)
        swapped_pka1 = group1.pka_value
        swapped_pka2 = group2.pka_value
        pka_shift1 = swapped_pka1 - default_pka1
        pka_shift2 = swapped_pka2 - default_pka2
        # Swap back to original protonation state
        self.swap_interactions([group1], [group2])
        group1.calculate_total_pka()
        group2.calculate_total_pka()
        # check difference in free energy
        if (abs(default_energy - swapped_energy)
                > self.parameters.max_free_energy_diff and return_on_fail):
            return {'coupling_factor': -1.0}
        # check pka shift
        if (max(abs(pka_shift1), abs(pka_shift2))
                < self.parameters.min_swap_pka_shift and return_on_fail):
            return {'coupling_factor': -1.0}
        # check intrinsic pka diff
        if (abs(group1.intrinsic_pka - group2.intrinsic_pka)
                > self.parameters.max_intrinsic_pka_diff and return_on_fail):
            return {'coupling_factor': -1.0}
        # if everything is OK, calculate the coupling factor and return all
        # info
        factor = (
            self.get_free_energy_diff_factor(default_energy, swapped_energy)
            * self.get_pka_diff_factor(group1.intrinsic_pka,
                                       group2.intrinsic_pka)
            * self.get_interaction_factor(interaction_energy))
        return {'coupling_factor': factor, 'default_energy': default_energy,
                'swapped_energy': swapped_energy,
                'interaction_energy': interaction_energy,
                'swapped_pka1': swapped_pka1, 'swapped_pka2': swapped_pka2,
                'pka_shift1': pka_shift1, 'pka_shift2': pka_shift2,
                'pH': use_ph}

    def get_pka_diff_factor(self, pka1, pka2):
        """Get scaling factor for difference between intrinsic pKa values.

        Args:
            pka1:  first pKa to compare
            pka2:  second pKa to compare
        Returns:
            float value of scaling factor
        """
        assert self.parameters is not None
        intrinsic_pka_diff = abs(pka1-pka2)
        res = 0.0
        if intrinsic_pka_diff <= self.parameters.max_intrinsic_pka_diff:
            res = (
                1-(intrinsic_pka_diff
                   / self.parameters.max_intrinsic_pka_diff)**2)
        return res

    def get_free_energy_diff_factor(self, energy1, energy2):
        """Get scaling factor for difference between free energies.

        Args:
            energy1:  first energy to compare
            energy2:  second energy to compare
        Returns:
            float value of scaling factor
        """
        assert self.parameters is not None
        free_energy_diff = abs(energy1-energy2)
        res = 0.0
        if free_energy_diff <= self.parameters.max_free_energy_diff:
            res = 1-(free_energy_diff/self.parameters.max_free_energy_diff)**2
        return res

    def get_interaction_factor(self, interaction_energy):
        """Get scaling factor related to interaction energy.

        Args:
            interaction_energy:  interaction energy
        Returns:
            float value of scaling factor
        """
        assert self.parameters is not None
        res = 0.0
        interaction_energy = abs(interaction_energy)
        if interaction_energy >= self.parameters.min_interaction_energy:
            res = (
                (interaction_energy-self.parameters.min_interaction_energy)
                / (1.0+interaction_energy
                   - self.parameters.min_interaction_energy))
        return res

    def identify_non_covalently_coupled_groups(self, conformation,
                                               verbose=True):
        """Find coupled residues in protein.

        Args:
            conformation:  protein conformation to test
            verbose:  verbose output (boolean)
        """
        self.parameters = conformation.parameters
        if verbose:
            info_fmt = (
                '\n'
                ' Warning: When using the -d option, pKa values based on \n'
                '\'swapped\' interactions\n'
                '          will be writting to the output .pka file\n'
                '\n'
                '{sep}\n'
                '\n'
                ' Detecting non-covalently coupled residues\n'
                '{sep}\n'
                '   Maximum pKa difference:     '
                '{c.max_intrinsic_pka_diff:>4.2f} pKa units\n'
                '   Minimum interaction energy: '
                '{c.min_interaction_energy:>4.2f} pKa units\n'
                '   Maximum free energy diff.:  '
                '{c.max_free_energy_diff:>4.2f} pKa units\n'
                '   Minimum swap pKa shift:     '
                '{c.min_swap_pka_shift:>4.2f} pKa units\n'
                '   pH:                         '
                '{c.pH:>6} \n'
                '   Reference:                  '
                '{c.reference}\n'
                '   Min pKa:                    '
                '{c.min_pka:>4.2f}\n'
                '   Max pKa:                    '
                '{c.max_pka:>4.2f}\n'
                '\n')
            sep = "-" * 103
            _LOGGER.info(info_fmt.format(sep=sep, c=self.parameters))
        # find coupled residues
        titratable_groups = conformation.get_titratable_groups()
        if not conformation.non_covalently_coupled_groups:
            for group1 in titratable_groups:
                for group2 in titratable_groups:
                    if group1 == group2:
                        break
                    if (group1 not in group2.non_covalently_coupled_groups
                            and self.do_prot_stat):
                        data = (
                            self
                            .is_coupled_protonation_state_probability(
                                group1, group2,
                                conformation.calculate_folding_energy))
                        if data['coupling_factor'] > 0.0:
                            group1.couple_non_covalently(group2)
        if verbose:
            self.print_out_swaps(conformation)

    def print_out_swaps(self, conformation):
        """Print out something having to do with coupling interactions.

        Args:
            conformation:  conformation to print
        """
        map_ = make_interaction_map(
            'Non-covalent coupling map for {0:s}'.format(str(conformation)),
            conformation.get_non_covalently_coupled_groups(),
            lambda g1, g2: g1 in g2.non_covalently_coupled_groups)
        _LOGGER.info(map_)
        for system in conformation.get_coupled_systems(
                conformation.get_non_covalently_coupled_groups(),
                Group.get_non_covalently_coupled_groups):
            self.print_system(conformation, list(system))

    def print_system(self, conformation, system):
        """Print out something about the system.

        Args:
            conformation:  conformation to print
            system:  system to print
        """
        _LOGGER.info(
            'System containing {0:d} groups:'.format(len(system)))
        # make list of interactions within this system
        interactions = list(itertools.combinations(system, 2))
        # print out coupling info for each interaction
        coup_info = ''
        for interaction in interactions:
            data = (
                self.is_coupled_protonation_state_probability(
                    interaction[0], interaction[1],
                    conformation.calculate_folding_energy,
                    return_on_fail=False))
            coup_info += (
                self.make_data_to_string(data, interaction[0], interaction[1])
                + '\n\n')
        _LOGGER.info(coup_info)
        # make list of possible combinations of swap to try out
        combinations = propka.lib.generate_combinations(interactions)
        # Make possible swap combinations
        swap_info = ''
        swap_info += self.print_determinants_section(system, 'Original')
        for combination in combinations:
            # Tell the user what is swap in this combination
            swap_info += 'Swapping the following interactions:\n'
            for interaction in combination:
                swap_info += ' {0:s} {1:s}\n'.format(
                    interaction[0].label, interaction[1].label)
            # swap...
            for interaction in combination:
                self.swap_interactions([interaction[0]], [interaction[1]])
            swap_info += self.print_determinants_section(system, 'Swapped')
        _LOGGER.info(swap_info)

    @staticmethod
    def get_interaction(group1: Group, group2: Group, include_side_chain_hbs=True):
        """Get interaction energy between two groups.

        Args:
            group1:  first group for interaction
            group2:  second group for interaction
            include_side_chain_hbs:  include sidechain hydrogen bonds in energy
        Returns:
            interaction energy (float)
        """
        determinants = group1.determinants['coulomb']
        if include_side_chain_hbs:
            determinants = (
                group1.determinants['sidechain']
                + group1.determinants['coulomb'])
        interaction_energy = 0.0
        for det in determinants:
            if group2 == det.group:
                interaction_energy += det.value
        return interaction_energy

    def print_determinants_section(self, system, tag):
        """Print determinants of system.

        Args:
            system:  set of groups
            tag:  something to add to output
        Returns:
            string with summary
        """
        all_labels = [g.label for g in system]
        str_ = ' ' + '-' * 113 + '\n'
        for group in system:
            str_ += self.tagged_format(
                ' {0:<8s}|'.format(tag), group.get_determinant_string(),
                all_labels)
        return str_ + '\n'

    def swap_interactions(self, groups1, groups2, include_side_chain_hbs=True):
        """Swap interactions between two groups.

        Args:
            group1:  first group to swap
            group2:  second group to swap
            """
        for i, group1 in enumerate(groups1):
            group2 = groups2[i]
            # swap the interactions!
            self.transfer_determinant(group1.determinants['coulomb'],
                                      group2.determinants['coulomb'],
                                      group1.label, group2.label)
            if include_side_chain_hbs:
                self.transfer_determinant(group1.determinants['sidechain'],
                                          group2.determinants['sidechain'],
                                          group1.label, group2.label)
                # re-calculate pKa values
                group1.calculate_total_pka()
                group2.calculate_total_pka()

    @staticmethod
    def transfer_determinant(determinants1, determinants2,
                             label1, label2):
        """Transfer information between two sets of determinants.

        Args:
            determinants1:  determinant list
            determinants2:  determinant list
            label1:  label for list 1
            label2:  label for list 2
        """
        # find out what to transfer...
        from1to2 = []
        from2to1 = []
        for det in determinants1:
            if det.label == label2:
                from1to2.append(det)
        for det in determinants2:
            if det.label == label1:
                from2to1.append(det)
        # ...and transfer it!
        for det in from1to2:
            det.label = label1
            determinants2.append(det)
            determinants1.remove(det)
        for det in from2to1:
            det.label = label2
            determinants1.append(det)
            determinants2.remove(det)

    @staticmethod
    def tagged_format(tag, str_, labels):
        """Tag a string.

        Args:
            tag:  tag to add
            str_:  string to tag
            labels:  labels to replace
        Returns:
            tagged string
        """
        str_ = "{0:s} {1:s}".format(tag, str_)
        str_ = str_.replace('\n', '\n{0:s} '.format(tag))
        for label in labels:
            str_ = str_.replace(label, '\033[31m{0:s}\033[30m'.format(label))
        return str_ + '\n'

    @staticmethod
    def make_data_to_string(data, group1, group2):
        """Describe interaction between groups.

        Args:
            data:  data about interactions
            group1:  first group
            group2:  second group
        Returns:
            formatted string with information.
        """
        str_ = (
            " {label1} and {label2} coupled (prot.state): {coupl_fact:>5.2f}\n"
            " Energy levels:       {def_energy:>6.2f}, {swap_energy:>6.2f}  "
            "(difference: {diff_energy:>6.2f}) at pH {ph:>6.2f}\n"
            " Interaction energy:  {int_energy:>6.2f}\n"
            " Intrinsic pka's:     {pka1:>6.2f}, {pka2:>6.2f}  "
            "(difference: {diff_pka:>6.2f})\n"
            " Swapped pKa's:       {swap1:>6.2f}, {swap2:>6.2f}  "
            "(difference: {shift1:>6.2f}, {shift2:>6.2f})"
        ).format(
            label1=group1.label, label2=group2.label,
            coupl_fact=data['coupling_factor'],
            def_energy=data['default_energy'],
            swap_energy=data['swapped_energy'],
            diff_energy=data['default_energy']-data['swapped_energy'],
            ph=data['pH'], int_energy=data['interaction_energy'],
            pka1=group1.intrinsic_pka, pka2=group2.intrinsic_pka,
            diff_pka=group1.intrinsic_pka-group2.intrinsic_pka,
            swap1=data['swapped_pka1'], swap2=data['swapped_pka2'],
            shift1=data['pka_shift1'], shift2=data['pka_shift2'])
        return str_


NCCG = NonCovalentlyCoupledGroups()
