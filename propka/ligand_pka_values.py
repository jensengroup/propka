"""
Ligand pKa values from Marvin
=============================

Ligand pKa values can be obtained from the commercial `Marvin`_
software (namely, the :program:`cxcalc` and :program:`molconvert`
programs are required).

.. _Marvin: https://chemaxon.com/products/marvin

"""
import logging
import os
import shutil
import subprocess
import sys
import warnings
from typing import TYPE_CHECKING, NoReturn

from propka.output import write_mol2_for_atoms
from propka.lib import split_atoms_into_molecules
from propka.parameters import Parameters

if TYPE_CHECKING:
    from propka.molecular_container import MolecularContainer

_LOGGER = logging.getLogger(__name__)


class LigandPkaValues:
    """Ligand pKa value class."""

    def __init__(self, parameters: Parameters):
        """Initialize object with parameters.

        Args:
            parameters:  parameters
        """
        self.parameters = parameters
        # attempt to find Marvin executables in the path
        self.molconvert = self.find_in_path('molconvert')
        self.cxcalc = self.find_in_path('cxcalc')
        _LOGGER.info('Found Marvin executables:')
        _LOGGER.info(self.cxcalc)
        _LOGGER.info(self.molconvert)

    @staticmethod
    def find_in_path(program):
        """Find a program in the system path.

        Args:
            program:  program to find
        Returns:
            location of program
        """
        loc = shutil.which(program)
        if loc is None:
            str_ = "'Error: Could not find {0:s}.".format(program)
            str_ += ' Please make sure that it is found in the path.'
            _LOGGER.info(str_)
            sys.exit(-1)
        return loc

    def get_marvin_pkas_for_pdb_file(
            self, molecule: "MolecularContainer", parameters: NoReturn,
            num_pkas=10, min_ph=-10.0, max_ph=20.0):
        """Use Marvin executables to get pKas for a PDB file.

        Args:
            pdbfile:  PDB file
            molecule:  MolecularContainer object
            num_pkas:  number of pKas to get
            min_ph:  minimum pH value
            max_ph:  maximum pH value
        """
        warnings.warn("unused and untested by propka")
        self.get_marvin_pkas_for_molecular_container(
            molecule, num_pkas=num_pkas, min_ph=min_ph, max_ph=max_ph)

    def get_marvin_pkas_for_molecular_container(self, molecule: "MolecularContainer", num_pkas=10,
                                                min_ph=-10.0, max_ph=20.0):
        """Use Marvin executables to calculate pKas for a molecular container.

        Args:
            molecule:  molecular container
            num_pkas:  number of pKas to calculate
            min_ph:  minimum pH value
            max_ph:  maximum pH value
        """
        for name in molecule.conformation_names:
            filename = '{0:s}_{1:s}'.format(molecule.name, name)
            self.get_marvin_pkas_for_conformation_container(
                molecule.conformations[name], name=filename,
                reuse=molecule.options.reuse_ligand_mol2_file,
                num_pkas=num_pkas, min_ph=min_ph, max_ph=max_ph)

    def get_marvin_pkas_for_conformation_container(self, conformation,
                                                   name='temp', reuse=False,
                                                   num_pkas=10, min_ph=-10.0,
                                                   max_ph=20.0):
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
                                  num_pkas=10, min_ph=-10.0, max_ph=20.0):
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
        molecules = split_atoms_into_molecules(atoms)
        for i, molecule in enumerate(molecules):
            filename = '{0:s}_{1:d}.mol2'.format(name, i+1)
            self.get_marvin_pkas_for_molecule(
                molecule, filename=filename, reuse=reuse, num_pkas=num_pkas,
                min_ph=min_ph, max_ph=max_ph)

    def get_marvin_pkas_for_molecule(self, atoms, filename='__tmp_ligand.mol2',
                                     reuse=False, num_pkas=10, min_ph=-10.0,
                                     max_ph=20.0):
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
            write_mol2_for_atoms(atoms, filename)
        # check that we actually have a file to work with
        if not os.path.isfile(filename):
            errstr = (
                "Didn't find a user-modified file '{0:s}' "
                "- generating one".format(
                    filename))
            _LOGGER.warning(errstr)
            write_mol2_for_atoms(atoms, filename)
        # Marvin calculate pKa values
        fmt = (
            'pka -a {num1} -b {num2} --min {min_ph} '
            '--max {max_ph} -d large')
        options = (
            fmt.format(
                num1=num_pkas, num2=num_pkas, min_ph=min_ph, max_ph=max_ph))
        (output, errors) = subprocess.Popen(
            [self.cxcalc, filename]+options.split(), stdout=subprocess.PIPE,
            stderr=subprocess.PIPE).communicate()
        if len(errors) > 0:
            err = (
                f'Error: Marvin execution failed:  {errors}\n'
                f'Please edit the ligand mol2 file and re-run PropKa with '
                f'the -l option: {filename}'
            )
            _LOGGER.error(err)
            raise OSError(err)
        # extract calculated pkas
        indices, pkas, types = self.extract_pkas(output)
        # store calculated pka values
        for i, index in enumerate(indices):
            atoms[index].marvin_pka = pkas[i]
            atoms[index].charge = {'a': -1, 'b': 1}[types[i]]
            _LOGGER.info(
                '{0:s} model pKa: {1:<.2f}'.format(atoms[index], pkas[i])
            )

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
