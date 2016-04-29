#!/usr/bin/env python

from __future__ import division
from __future__ import print_function
from propka.lib import info, warning

import propka.molecular_container, propka.calculations, propka.calculations, propka.parameters, propka.pdb, propka.lib, os, subprocess, sys

class ligand_pka_values:
    def __init__(self, parameters):
        self.parameters = parameters

        # attempt to find Marvin executables in the path
        self.molconvert = self.find_in_path('molconvert')
        self.cxcalc = self.find_in_path('cxcalc')
        info('Found Marvin executables:')
        info(self.cxcalc)
        info(self.molconvert)

        return


    def find_in_path(self, program):
        path = os.environ.get('PATH').split(os.pathsep)

        l = [i for i in filter(lambda loc: os.access(loc, os.F_OK),
                               map(lambda dir: os.path.join(dir, program),path))]

        if len(l) == 0:
            info('Error: Could not find %s. Please make sure that it is found in the path.' % program)
            sys.exit(-1)

        return l[0]

    def get_marvin_pkas_for_pdb_file(self, file, no_pkas=10, min_pH =-10, max_pH=20):
        molecule = propka.molecular_container.Molecular_container(file)
        self.get_marvin_pkas_for_molecular_container(molecule, no_pkas=no_pkas, min_pH =min_pH, max_pH=max_pH)
        return

    def get_marvin_pkas_for_molecular_container(self, molecule, no_pkas=10, min_pH =-10, max_pH=20):
        for name in molecule.conformation_names:
            filename = '%s_%s'%(molecule.name,name)
            self.get_marvin_pkas_for_conformation_container(molecule.conformations[name], name=filename, reuse=molecule.options.reuse_ligand_mol2_file,
                                                            no_pkas=no_pkas, min_pH =min_pH, max_pH=max_pH)

        return

    def get_marvin_pkas_for_conformation_container(self, conformation, name = 'temp', reuse=False, no_pkas=10, min_pH =-10, max_pH=20):
        conformation.marvin_pkas_calculated = True
        self.get_marvin_pkas_for_atoms(conformation.get_heavy_ligand_atoms(), name=name, reuse=reuse,
                                       no_pkas=no_pkas, min_pH =min_pH, max_pH=max_pH)

        return

    def get_marvin_pkas_for_atoms(self, atoms, name='temp', reuse=False, no_pkas=10, min_pH =-10, max_pH=20):
        # do one molecule at the time so we don't confuse marvin
        molecules = propka.lib.split_atoms_into_molecules(atoms)
        for i in range(len(molecules)):
            filename = '%s_%d.mol2'%(name, i+1)
            self.get_marvin_pkas_for_molecule(molecules[i], filename=filename, reuse=reuse, no_pkas=no_pkas, min_pH =min_pH, max_pH=max_pH)

        return


    def get_marvin_pkas_for_molecule(self, atoms, filename='__tmp_ligand.mol2', reuse=False, no_pkas=10, min_pH =-10, max_pH=20):
        # print out structure unless we are using user-modified structure
        if not reuse:
            propka.pdb.write_mol2_for_atoms(atoms, filename)
        # check that we actually have a file to work with
        if not os.path.isfile(filename):
            warning('Didn\'t find a user-modified file \'%s\' - generating one' % filename)
            propka.pdb.write_mol2_for_atoms(atoms, filename)



        # Marvin
        # calculate pKa values
        options = 'pka -a %d -b %d --min %f --max %f -d large'%(no_pkas, no_pkas, min_pH, max_pH)
        (output,errors) = subprocess.Popen([self.cxcalc, filename]+options.split(),
                                  stdout=subprocess.PIPE, stderr=subprocess.PIPE).communicate()

        if len(errors)>0:
            info('********************************************************************************************************')
            info('* Warning: Marvin execution failed:                                                                    *')
            info('* %-100s *' % errors)
            info('*                                                                                                      *')
            info('* Please edit the ligand mol2 file and re-run PropKa with the -l option: %29s *' % filename)
            info('********************************************************************************************************')
            sys.exit(-1)

        # extract calculated pkas
        indices,pkas,types = self.extract_pkas(output)

        # store calculated pka values
        for i in range(len(indices)):
            atoms[indices[i]].marvin_pka = pkas[i]
            atoms[indices[i]].charge = {'a':-1,'b':+1}[types[i]]
            info('%s model pKa: %.2f' % (atoms[indices[i]], pkas[i]))

        return

    def extract_pkas(self, output):
        # split output
        [tags, values,empty_line] = output.decode().split('\n')
        #info(tags)
        #info(values)
        tags = tags.split('\t')
        values = values.split('\t')

        # format values
        types = [tags[i][0] for i in range(1,len(tags)-1) if len(values)>i and values[i] != '']
        indices = [int(a)-1 for a in values[-1].split(',') if a !='']
        values = [float(v.replace(',','.')) for v in values[1:-1] if v != '']

        if len(indices) != len(values) != len(types):
            raise Exception('Lengths of atoms and pka values mismatch')

        return indices, values, types




