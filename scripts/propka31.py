#!/usr/bin/env python
"""PROPKA script.

This is the original propka script. However, this distribute-based
installation moved the main() function into propka.run.main and just
generates a script called propka31 from the setup.py installation
script. You should not need to use this script.

(Also note that there can be import problems because the script name
is the same as the module name; that's why the new script is called
propka31.)
"""
from propka.lib import loadOptions
from propka.input import read_parameter_file, read_molecule_file
from propka.parameters import Parameters
from propka.molecular_container import MolecularContainer


def main():
    """Read in structure files, calculates pKa values, and prints pKa files."""
    # loading options, flaggs and arguments
    options = loadOptions([])
    pdbfiles = options.filenames
    parameters = read_parameter_file(options.parameters, Parameters())

    for pdbfile in pdbfiles:
        my_molecule = MolecularContainer(parameters, options)
        my_molecule = read_molecule_file(pdbfile, my_molecule)
        my_molecule.calculate_pka()
        my_molecule.write_pka()


if __name__ == '__main__':
    main()
