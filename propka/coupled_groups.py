
from __future__ import division
from __future__ import print_function

import math, propka.output, propka.group, propka.lib, itertools
from propka.lib import info, warning


class non_covalently_couple_groups:
    def __init__(self):

        self.parameters = None

        # self.do_intrinsic = False
        # self.do_pair_wise = False
        self.do_prot_stat = True

        return


    #
    # Methods for finding coupled groups
    #
    def is_coupled_protonation_state_probability(self, group1, group2, energy_method, return_on_fail=True):

        # check if the interaction energy is high enough
        interaction_energy = max(self.get_interaction(group1,group2), self.get_interaction(group2,group1))

        if interaction_energy<=self.parameters.min_interaction_energy and return_on_fail:
            return {'coupling_factor':-1.0}

         # calculate intrinsic pKa's, if not already done
        for group in [group1, group2]:
            if not hasattr(group, 'intrinsic_pKa'):
                group.calculate_intrinsic_pka()

        use_pH = self.parameters.pH
        if self.parameters.pH == 'variable':
            use_pH = min(group1.pka_value, group2.pka_value)

        default_energy = energy_method(pH=use_pH, reference=self.parameters.reference)
        default_pka1 = group1.pka_value
        default_pka2 = group2.pka_value

        # check that pka values are within relevant limits
        if max(default_pka1, default_pka2) < self.parameters.min_pka or \
                min(default_pka1, default_pka2) > self.parameters.max_pka:
            if return_on_fail:
                return {'coupling_factor':-1.0}

        # Swap interactions and re-calculate pKa values
        self.swap_interactions([group1], [group2])
        group1.calculate_total_pka()
        group2.calculate_total_pka()

        # store swapped energy and pka's
        swapped_energy = energy_method(pH=use_pH, reference=self.parameters.reference)
        swapped_pka1 = group1.pka_value
        swapped_pka2 = group2.pka_value

        pka_shift1 = swapped_pka1 - default_pka1
        pka_shift2 = swapped_pka2 - default_pka2

        # Swap back to original protonation state
        self.swap_interactions([group1], [group2])
        group1.calculate_total_pka()
        group2.calculate_total_pka()

        # check difference in free energy
        if abs(default_energy - swapped_energy) > self.parameters.max_free_energy_diff and return_on_fail:
            return {'coupling_factor':-1.0}

        # check pka shift
        if max(abs(pka_shift1), abs(pka_shift2)) < self.parameters.min_swap_pka_shift and return_on_fail:
            return {'coupling_factor':-1.0}

        # check intrinsic pka diff
        if abs(group1.intrinsic_pKa - group2.intrinsic_pKa) > self.parameters.max_intrinsic_pKa_diff and return_on_fail:
            return {'coupling_factor':-1.0}

        # if everything is OK, calculate the coupling factor and return all info
        factor = self.get_free_energy_diff_factor(default_energy, swapped_energy)*\
            self.get_pka_diff_factor(group1.intrinsic_pKa, group2.intrinsic_pKa)*\
            self.get_interaction_factor(interaction_energy)

        return {'coupling_factor':factor,
                'default_energy':default_energy,
                'swapped_energy':swapped_energy,
                'interaction_energy':interaction_energy,
                'swapped_pka1':swapped_pka1,
                'swapped_pka2':swapped_pka2,
                'pka_shift1':pka_shift1,
                'pka_shift2':pka_shift2,
                'pH':use_pH}



    #
    # Methods for calculating the coupling factor
    #

    def get_pka_diff_factor(self, pka1, pka2):
        intrinsic_pka_diff = abs(pka1-pka2)
        res = 0.0
        if intrinsic_pka_diff <= self.parameters.max_intrinsic_pKa_diff:
            res = 1-(intrinsic_pka_diff/self.parameters.max_intrinsic_pKa_diff)**2

        return res

    def get_free_energy_diff_factor(self, energy1, energy2):
        free_energy_diff = abs(energy1-energy2)
        res = 0.0
        if free_energy_diff <= self.parameters.max_free_energy_diff:
            res = 1-(free_energy_diff/self.parameters.max_free_energy_diff)**2
        return res

    def get_interaction_factor(self, interaction_energy):
        res = 0.0
        interaction_energy = abs(interaction_energy)
        if interaction_energy >= self.parameters.min_interaction_energy:
            res = (interaction_energy-self.parameters.min_interaction_energy)/(1.0+interaction_energy-self.parameters.min_interaction_energy)

        return res




    def identify_non_covalently_coupled_groups(self, conformation, verbose=True):
        """ Finds coupled residues in protein """

        self.parameters = conformation.parameters
        if verbose:
            info('')
            info(' Warning: When using the -d option, pKa values based on \'swapped\' interactions')
            info('          will be writting to the output .pka file')
            info('')
            info('-' * 103)
            info(' Detecting non-covalently coupled residues')
            info('-' * 103)
            info('   Maximum pKa difference:     %4.2f pKa units' % self.parameters.max_intrinsic_pKa_diff)
            info('   Minimum interaction energy: %4.2f pKa units' % self.parameters.min_interaction_energy)
            info('   Maximum free energy diff.:  %4.2f pKa units' % self.parameters.max_free_energy_diff)
            info('   Minimum swap pKa shift:     %4.2f pKa units' % self.parameters.min_swap_pka_shift)
            info('   pH:                         %6s ' % str(self.parameters.pH))
            info('   Reference:                  %s' % self.parameters.reference)
            info('   Min pKa:                    %4.2f' % self.parameters.min_pka)
            info('   Max pKa:                    %4.2f' % self.parameters.max_pka)
            info('')

        # find coupled residues
        titratable_groups = conformation.get_titratable_groups()

        if not conformation.non_covalently_coupled_groups:
            for g1 in titratable_groups:
                for g2 in titratable_groups:
                    if g1==g2:
                        break

                    if not g1 in g2.non_covalently_coupled_groups and self.do_prot_stat:
                        data = self.is_coupled_protonation_state_probability(g1, g2,conformation.calculate_folding_energy)
                        if data['coupling_factor'] >0.0:
                            g1.couple_non_covalently(g2)

        if verbose:
            self.print_out_swaps(conformation)

        return

    def print_out_swaps(self, conformation, verbose=True):
        map = propka.output.make_interaction_map('Non-covalent coupling map for %s'%conformation,
                                                 conformation.get_non_covalently_coupled_groups(),
                                                 lambda g1,g2: g1 in g2.non_covalently_coupled_groups)
        info(map)

        for system in conformation.get_coupled_systems(conformation.get_non_covalently_coupled_groups(),propka.group.Group.get_non_covalently_coupled_groups):
            self.print_system(conformation, list(system))
        return

    def print_system(self, conformation, system):

        info('System containing %d groups:' % len(system))

        # make list of interactions withi this system
        interactions = list(itertools.combinations(system,2))

        # print out coupling info for each interaction
        coup_info = ''
        for interaction in interactions:
            data = self.is_coupled_protonation_state_probability(interaction[0], interaction[1],conformation.calculate_folding_energy, return_on_fail=False)
            coup_info += self.make_data_to_string(data,interaction[0],interaction[1])+'\n\n'
        info(coup_info)

        # make list of possible combinations of swap to try out
        combinations = propka.lib.generate_combinations(interactions)

        # Make possible swap combinations
        swap_info = ''
        swap_info += self.print_determinants_section(system,'Original')

        for combination in combinations:
            # Tell the user what is swap in this combination
            swap_info += 'Swapping the following interactions:\n'
            for interaction in combination:
                swap_info += ' %s %s\n'%(interaction[0].label, interaction[1].label)

            # swap...
            for interaction in combination:
                self.swap_interactions([interaction[0]],[interaction[1]])

            swap_info += self.print_determinants_section(system,'Swapped')
            # ...and swap back
            #for interaction in combination:
            #    self.swap_interactions([interaction[0]], [interaction[1]])

        info(swap_info)

        return
    #
    # Interaction and swapping methods
    #

    def get_interaction(self, group1, group2, include_side_chain_hbs = True):
        determinants = group1.determinants['coulomb']
        if include_side_chain_hbs:
            determinants = group1.determinants['sidechain'] + group1.determinants['coulomb']

        interaction_energy = 0.0
        for det in determinants:
            if group2 == det.group:
                interaction_energy += det.value

        return interaction_energy



    def print_determinants_section(self, system, tag):
        all_labels = [g.label for g in system]
        s = ' '+'-'*113+'\n'
        for group in system:
            s += self.tagged_format(' %-8s|' % tag, group.getDeterminantString(), all_labels)

        return s+'\n'

    def swap_interactions(self, groups1, groups2, include_side_chain_hbs = True):

        for i in range(len(groups1)):
            group1 = groups1[i]
            group2 = groups2[i]

            # swap the interactions!
            self.transfer_determinant(group1.determinants['coulomb'], group2.determinants['coulomb'], group1.label, group2.label)
            if include_side_chain_hbs:
                self.transfer_determinant(group1.determinants['sidechain'], group2.determinants['sidechain'], group1.label, group2.label)

            # re-calculate pKa values
                group1.calculate_total_pka()
                group2.calculate_total_pka()

        return


    def transfer_determinant(self, determinants1, determinants2, label1, label2):
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

        return

    #
    # Output methods
    #

    def tagged_format(self, tag, s, labels):
        s = "%s %s"%(tag,s)
        s = s.replace('\n','\n%s '%tag)
        for label in labels:
            s = s.replace(label, '\033[31m%s\033[30m'%label)
        return s+'\n'


    def make_data_to_string(self, data, group1, group2):
        s = """ %s and %s coupled (prot.state): %5.2f
 Energy levels:       %6.2f, %6.2f  (difference: %6.2f) at pH %6.2f
 Interaction energy:  %6.2f
 Intrinsic pka's:     %6.2f, %6.2f  (difference: %6.2f)
 Swapped pKa's:       %6.2f, %6.2f  (difference: %6.2f, %6.2f)"""%(group1.label,
                                                                   group2.label,
                                                                   data['coupling_factor'],
                                                                   data['default_energy'], data['swapped_energy'],
                                                                   data['default_energy'] - data['swapped_energy'],
                                                                   data['pH'],
                                                                   data['interaction_energy'],
                                                                   group1.intrinsic_pKa,
                                                                   group2.intrinsic_pKa,
                                                                   group1.intrinsic_pKa-group2.intrinsic_pKa,
                                                                   data['swapped_pka1'],
                                                                   data['swapped_pka2'],
                                                                   data['pka_shift1'],
                                                                   data['pka_shift2'])

        return s


nccg = non_covalently_couple_groups()
