# entry point for propka script

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
