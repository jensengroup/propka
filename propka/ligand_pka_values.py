"""Ligand pKa values"""
import os
import subprocess
import sys
import propka.molecular_container
import propka.calculations
import propka.parameters
import propka.pdb
import propka.lib
from propka.lib import info, warning


class LigandPkaValues:
    """Ligand pKa value class."""

    def __init__(self, parameters):
        """Initialize object with parameters.

        Args:
            parameters:  parameters
        """
        self.parameters = parameters
        # attempt to find Marvin executables in the path
        self.molconvert = self.find_in_path('molconvert')
        self.cxcalc = self.find_in_path('cxcalc')
        info('Found Marvin executables:')
        info(self.cxcalc)
        info(self.molconvert)

    @staticmethod
    def find_in_path(program):
        """Find a program in the system path.

        Args:
            program:  program to find
        Returns:
            location of program
        """
        path = os.environ.get('PATH').split(os.pathsep)
        locs = [
            i for i in filter(lambda loc: os.access(loc, os.F_OK),
                              map(lambda dir: os.path.join(dir, program),
                                  path))]
        if len(locs) == 0:
            str_ = "'Error: Could not find %s." % program
            str_ += ' Please make sure that it is found in the path.'
            info(str_)
            sys.exit(-1)
        return locs[0]

    def get_marvin_pkas_for_pdb_file(self, pdbfile, num_pkas=10, min_ph=-10,
                                     max_ph=20):
        """Use Marvin executables to get pKas for a PDB file.

        Args:
            pdbfile:  PDB file
            num_pkas:  number of pKas to get
            min_ph:  minimum pH value
            max_ph:  maximum pH value
        """
        molecule = propka.molecular_container.Molecular_container(pdbfile)
        self.get_marvin_pkas_for_molecular_container(
            molecule, num_pkas=num_pkas, min_ph=min_ph, max_ph=max_ph)

    def get_marvin_pkas_for_molecular_container(self, molecule, num_pkas=10,
                                                min_ph=-10, max_ph=20):
        """Use Marvin executables to calculate pKas for a molecular container.

        Args:
            molecule:  molecular container
            num_pkas:  number of pKas to calculate
            min_ph:  minimum pH value
            max_ph:  maximum pH value
        """
        for name in molecule.conformation_names:
            filename = '%s_%s' % (molecule.name, name)
            self.get_marvin_pkas_for_conformation_container(
                molecule.conformations[name], name=filename,
                reuse=molecule.options.reuse_ligand_mol2_file,
                num_pkas=num_pkas, min_ph=min_ph, max_ph=max_ph)

    def get_marvin_pkas_for_conformation_container(self, conformation,
                                                   name='temp', reuse=False,
                                                   num_pkas=10, min_ph=-10,
                                                   max_ph=20):
        """Use Marvin executables to calculate pKas for a conformation container.

        Args:
            conformation:  conformation container
            name:  filename
            reuse:  flag to reuse the structure files
            num_pkas:  number of pKas to calculate
            min_ph:  minimum pH value
            max_ph:  maximum pH value
        """
        conformation.marvin_pkas_calculated = True
        self.get_marvin_pkas_for_atoms(
            conformation.get_heavy_ligand_atoms(), name=name, reuse=reuse,
            num_pkas=num_pkas, min_ph=min_ph, max_ph=max_ph)

    def get_marvin_pkas_for_atoms(self, atoms, name='temp', reuse=False,
                                  num_pkas=10, min_ph=-10, max_ph=20):
        """Use Marvin executables to calculate pKas for a list of atoms.

        Args:
            atoms:  list of atoms
            name:  filename
            reuse:  flag to reuse the structure files
            num_pkas:  number of pKas to calculate
            min_ph:  minimum pH value
            max_ph:  maximum pH value
        """
        # do one molecule at the time so we don't confuse marvin
        molecules = propka.lib.split_atoms_into_molecules(atoms)
        for i, molecule in enumerate(molecules):
            filename = '%s_%d.mol2'%(name, i+1)
            self.get_marvin_pkas_for_molecule(
                molecule, filename=filename, reuse=reuse, num_pkas=num_pkas,
                min_ph=min_ph, max_ph=max_ph)

    def get_marvin_pkas_for_molecule(self, atoms, filename='__tmp_ligand.mol2',
                                     reuse=False, num_pkas=10, min_ph=-10,
                                     max_ph=20):
        """Use Marvin executables to calculate pKas for a molecule.

        Args:
            molecule:  the molecule
            name:  filename
            reuse:  flag to reuse the structure files
            num_pkas:  number of pKas to calculate
            min_ph:  minimum pH value
            max_ph:  maximum pH value
        """
        # print out structure unless we are using user-modified structure
        if not reuse:
            propka.pdb.write_mol2_for_atoms(atoms, filename)
        # check that we actually have a file to work with
        if not os.path.isfile(filename):
            errstr = ("Didn't find a user-modified file '%s' - generating one"
                      % filename)
            warning(errstr)
            propka.pdb.write_mol2_for_atoms(atoms, filename)
        # Marvin calculate pKa values
        options = ('pka -a %d -b %d --min %f --max %f -d large'
                   % (num_pkas, num_pkas, min_ph, max_ph))
        (output, errors) = subprocess.Popen(
            [self.cxcalc, filename]+options.split(), stdout=subprocess.PIPE,
            stderr=subprocess.PIPE).communicate()
        if len(errors) > 0:
            info('***********************************************************'
                 '*********************************************')
            info('* Warning: Marvin execution failed:                        '
                 '                                            *')
            info('* %-100s *' % errors)
            info('*                                                          '
                 '                                            *')
            info('* Please edit the ligand mol2 file and re-run PropKa with '
                 'the -l option: %29s *' % filename)
            info('***********************************************************'
                 '*********************************************')
            sys.exit(-1)
        # extract calculated pkas
        indices, pkas, types = self.extract_pkas(output)
        # store calculated pka values
        for i, index in enumerate(indices):
            atoms[index].marvin_pka = pkas[i]
            atoms[index].charge = {'a': -1, 'b': 1}[types[i]]
            info('%s model pKa: %.2f' % (atoms[index], pkas[i]))

    @staticmethod
    def extract_pkas(output):
        """Extract pKa value from output.

        Args:
            output:  output string to parse
        Returns:
            1. Indices
            2. Values
            3. Types
        """
        # split output
        [tags, values, _] = output.decode().split('\n')
        tags = tags.split('\t')
        values = values.split('\t')
        # format values
        types = [
            tags[i][0] for i in range(1, len(tags)-1) 
            if len(values) > i and values[i] != '']
        indices = [int(a)-1 for a in values[-1].split(',') if a != '']
        values = [float(v.replace(',', '.')) for v in values[1:-1] if v != '']
        if len(indices) != len(values) != len(types):
            raise Exception('Lengths of atoms and pka values mismatch')
        return indices, values, types
