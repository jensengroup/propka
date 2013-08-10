#!/usr/bin/env python

# This is the original propka script. However, this distribute-based
# installation moved the main() function into propka.run.main and just
# generates a script called propka31 from the setup.py installation
# script. You should not need to use this script.
#
# (Also note that there can be import problems because the script name
# is the same as the module name; that's why the new script is called
# propka31.)

import propka.lib, propka.molecular_container

def main():
    """
    Reads in structure files, calculates pKa values, and prints pKa files
    """
    # loading options, flaggs and arguments
    options, pdbfiles = propka.lib.loadOptions()

    for pdbfile in pdbfiles:
        my_molecule = propka.molecular_container.Molecular_container(pdbfile, options)
        my_molecule.calculate_pka()
        my_molecule.write_pka()


if __name__ == '__main__':
    #import cProfile
    #cProfile.run('main()',sort=1)
    main()

