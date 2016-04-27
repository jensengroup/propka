#!/usr/bin/python
#
# Molecular container for storing all contents of pdb files
#
#
from __future__ import division
from __future__ import print_function

import os, sys

import propka.pdb, propka.version, propka.output, propka.conformation_container, propka.group, propka.lib
from propka.lib import info, warning

class Molecular_container:
    def __init__(self, input_file, options=None):
        # printing out header before parsing input
        propka.output.printHeader()

        # set up some values
        self.options = options
        self.input_file = input_file
        self.dir = os.path.split(input_file)[0]
        self.file = os.path.split(input_file)[1]
        self.name = self.file[0:self.file.rfind('.')]
        input_file_extension = input_file[input_file.rfind('.'):]

        # set the version
        if options:
            parameters = propka.parameters.Parameters(self.options.parameters)
        else:
            parameters = propka.parameters.Parameters('propka.cfg')
        try:
            exec('self.version = propka.version.%s(parameters)'%parameters.version)
        except:
            raise Exception('Error: Version %s does not exist'%parameters.version)

        # read the input file
        if input_file_extension[0:4] == '.pdb':
            # input is a pdb file
            # read in atoms and top up containers to make sure that all atoms are present in all conformations
            [self.conformations, self.conformation_names] = propka.pdb.read_pdb(input_file, self.version.parameters,self)
            if len(self.conformations)==0:
                info('Error: The pdb file does not seems to contain any molecular conformations')
                sys.exit(-1)

            self.top_up_conformations()

            # make a structure precheck
            propka.pdb.protein_precheck(self.conformations, self.conformation_names)

            # set up atom bonding and protonation
            self.version.setup_bonding_and_protonation(self)

            # Extract groups
            self.extract_groups()

            # sort atoms
            for name in self.conformation_names:
                self.conformations[name].sort_atoms()

            # find coupled groups
            self.find_covalently_coupled_groups()

            # write out the input file
            filename = self.file.replace(input_file_extension,'.propka_input')
            propka.pdb.write_input(self, filename)


        elif input_file_extension == '.propka_input':
            #input is a propka_input file
            [self.conformations, self.conformation_names] = propka.pdb.read_input(input_file, self.version.parameters, self)

            # Extract groups - this merely sets up the groups found in the input file
            self.extract_groups()

            # do some additional set up
            self.additional_setup_when_reading_input_file()

        else:
            info('Unrecognized input file:%s' % input_file)
            sys.exit(-1)


        return
    def top_up_conformations(self):
        """ Makes sure that all atoms are present in all conformations """
        for name in self.conformation_names:
            if name!='1A' and (len(self.conformations[name]) < len(self.conformations['1A'])):
                self.conformations[name].top_up(self.conformations['1A'])

        return

    def find_covalently_coupled_groups(self):
        info('-' * 103)
        for name in self.conformation_names:
            self.conformations[name].find_covalently_coupled_groups()

        return

    def find_non_covalently_coupled_groups(self):
        info('-' * 103)
        for name in self.conformation_names:
            self.conformations[name].find_non_covalently_coupled_groups(verbose=self.options.display_coupled_residues)

        return

    def extract_groups(self):
        """ Identify the groups needed for pKa calculation """
        for name in self.conformation_names:
            self.conformations[name].extract_groups()

        return

    def additional_setup_when_reading_input_file(self):
        for name in self.conformation_names:
            self.conformations[name].additional_setup_when_reading_input_file()

        return


    def calculate_pka(self):
        # calculate for each conformation
        for name in self.conformation_names:
            self.conformations[name].calculate_pka(self.version, self.options)

        # find non-covalently coupled groups
        self.find_non_covalently_coupled_groups()

        # find the average of the conformations
        self.average_of_conformations()

        # print out the conformation-average results
        propka.output.printResult(self, 'AVR', self.version.parameters)

        return

    def average_of_conformations(self):
        # make a new configuration to hold the average values
        avr_conformation = propka.conformation_container.Conformation_container(name='average',
                                                                                parameters=self.conformations[self.conformation_names[0]].parameters,
                                                                                molecular_container=self)

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
                    warning('Group %s could not be found in conformation %s.' % (group.atom.residue_label, name))
            # ... and store the average value
            avr_group = avr_group / len(self.conformation_names)
            avr_conformation.groups.append(avr_group)

        # store information on coupling in the average container
        if len(list(filter(lambda c: c.non_covalently_coupled_groups, self.conformations.values()))):
            avr_conformation.non_covalently_coupled_groups = True

        # store chain info
        avr_conformation.chains = self.conformations[self.conformation_names[0]].chains

        self.conformations['AVR'] = avr_conformation
        return

    def write_pka(self, filename=None, reference="neutral", direction="folding", options=None):
        #for name in self.conformation_names:
        #    propka.output.writePKA(self, self.version.parameters, filename='%s_3.1_%s.pka'%(self.name, name),
        #                           conformation=name,reference=reference,
        #                           direction=direction, options=options)

        # write out the average conformation
        filename=os.path.join('%s.pka'%(self.name))

        # if the display_coupled_residues option is true,
        # write the results out to an alternative pka file
        if self.options.display_coupled_residues:
            filename=os.path.join('%s_alt_state.pka'%(self.name))

        if hasattr(self.version.parameters, 'output_file_tag') and len(self.version.parameters.output_file_tag)>0:
            filename=os.path.join('%s_%s.pka'%(self.name,self.version.parameters.output_file_tag))

        propka.output.writePKA(self, self.version.parameters, filename=filename,
                               conformation='AVR',reference=reference,
                               direction=direction, options=options)

        return

    def getFoldingProfile(self, conformation='AVR',reference="neutral", direction="folding", grid=[0., 14., 0.1], options=None):
        # calculate stability profile
        profile = []
        for ph in propka.lib.make_grid(*grid):
            ddg = self.conformations[conformation].calculate_folding_energy( pH=ph, reference=reference)
            #info(ph,ddg)
            profile.append([ph, ddg])

        # find optimum
        opt =[None, 1e6]
        for point in profile:
            opt = min(opt, point, key=lambda v:v[1])

        # find values within 80 % of optimum
        range_80pct = [None, None]
        values_within_80pct = [p[0] for p in profile if p[1]< 0.8*opt[1]]
        if len(values_within_80pct)>0:
            range_80pct = [min(values_within_80pct), max(values_within_80pct)]

        # find stability range
        stability_range = [None, None]
        stable_values = [p[0] for p in profile if p[1]< 0.0]

        if len(stable_values)>0:
            stability_range = [min(stable_values), max(stable_values)]

        return profile, opt, range_80pct, stability_range


    def getChargeProfile(self, conformation='AVR', grid=[0., 14., .1]):
        charge_profile = []
        for ph in propka.lib.make_grid(*grid):
            q_unfolded, q_folded = self.conformations[conformation].calculate_charge(self.version.parameters, pH=ph)
            charge_profile.append([ph, q_unfolded, q_folded])

        return charge_profile

    def getPI(self, conformation='AVR', grid=[0., 14., 1], iteration=0):
        #info('staring',grid, iteration)
        # search
        charge_profile = self.getChargeProfile(conformation=conformation, grid=grid)
        pi = []
        pi_folded = pi_unfolded = [None, 1e6,1e6]
        for point in charge_profile:
            pi_folded = min(pi_folded, point, key=lambda v:abs(v[2]))
            pi_unfolded = min(pi_unfolded, point, key=lambda v:abs(v[1]))

        # If results are not good enough, do it again with a higher sampling resolution
        pi_folded_value   = pi_folded[0]
        pi_unfolded_value = pi_unfolded[0]
        step = grid[2]
        if (pi_folded[2] > 0.01 or pi_unfolded[1] > 0.01) and iteration<4:
            pi_folded_value, x   = self.getPI(conformation=conformation, grid=[pi_folded[0]-step,   pi_folded[0]+step,   step/10.0], iteration=iteration+1)
            x, pi_unfolded_value = self.getPI(conformation=conformation, grid=[pi_unfolded[0]-step, pi_unfolded[0]+step, step/10.0], iteration=iteration+1)


        return pi_folded_value, pi_unfolded_value



if __name__ == '__main__':
    input_file = sys.argv[1]
    mc = Molecular_container(input_file)
