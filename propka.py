#!/usr/bin/python3

import string,re,sys,os,math
import Source.lib, Source.molecular_container 

def main():
    """
    Reads in structure files, calculates pKa values, and prints pKa files
    """
    # loading options, flaggs and arguments
    options, pdbfiles = Source.lib.loadOptions()

    for pdbfile in pdbfiles:
        my_molecule = Source.molecular_container.Molecular_container(pdbfile, options)
        my_molecule.calculate_pka()
        my_molecule.write_pka()


if __name__ == '__main__': 
    #import cProfile
    #cProfile.run('main()',sort=1)
    main()

