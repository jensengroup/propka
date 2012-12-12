#!/usr/bin/python

from __future__ import division
from __future__ import print_function

import string, sys, copy, math, os


#
# file I/O
#
def open_file_for_reading(filename):
    if not os.path.isfile(filename):
        raise Exception('Cannot find file %s' %filename)
    
    return open(filename,'r')

def open_file_for_writing(filename):
    res = open(filename,'w')
    if not res:
        raise Exception('Could not open %s'%filename)
    return res
    
#
# bookkeeping etc.
#
def conformation_sorter(conf):
    model = int(conf[:-1])
    altloc = conf[-1:]
    return model*100+ord(altloc)

def split_atoms_into_molecules(atoms):
    molecules = []

    while len(atoms)>0:
        initial_atom = atoms.pop()
        molecules.append( make_molecule(initial_atom,atoms))

    return molecules

def make_molecule(atom, atoms):
    bonded_atoms = [a for a in atoms if atom in a.bonded_atoms]
    res_atoms = [atom,]

    for ba in bonded_atoms:
        if ba in atoms:
            atoms.remove(ba)
            res_atoms.extend(make_molecule(ba, atoms))
         
    return res_atoms


def make_grid(min,max,step):
    x = min
    while x <= max:
        yield x
        x += step
    return

def generate_combinations(interactions):
    res = [[]]
    for interaction in interactions:
        res = make_combination(res, interaction)
    res.remove([])

    return res
    

def make_combination(combis, interaction):
    res = []
    for combi in combis:
        res.append(combi+[interaction])
        res.append(combi)
    return res



def loadOptions():
    """
    load the arguments parser with options
    """
    from optparse import OptionParser

    # defining a 'usage' message
    usage = "usage: %prog [options] filename"

    # creating a parser
    parser = OptionParser(usage)

    # loading the parser
    parser.add_option("-f", "--file", action="append", dest="filenames", 
           help="read data from <filename>, i.e. <filename> is added to arguments")
    parser.add_option("-r", "--reference", dest="reference", default="neutral", 
           help="setting which reference to use for stability calculations [neutral/low-pH]")
    parser.add_option("-c", "--chain", action="append", dest="chains", 
           help="creating the protein with only a specified chain, note, chains without ID are labeled 'A' [all]")
    parser.add_option("-t", "--thermophile", action="append", dest="thermophiles", 
           help="defining a thermophile filename; usually used in 'alignment-mutations'")
    parser.add_option("-a", "--alignment", action="append", dest="alignment", 
           help="alignment file connecting <filename> and <thermophile> [<thermophile>.pir]")
    parser.add_option("-m", "--mutation", action="append", dest="mutations", 
           help="specifying mutation labels which is used to modify <filename> according to, e.g. N25R/N181D")
    parser.add_option("-v", "--version", dest="version_label", default="Jan15", 
           help="specifying the sub-version of propka [Jan15/Dec19]")
    parser.add_option("-p", "--parameters",dest="parameters", default="propka.cfg",
                      help="set the parameter file")
    parser.add_option("-z", "--verbose", dest="verbose", action="store_true", default=True, 
           help="sleep during calculations")
    parser.add_option("-q", "--quiet", dest="verbose", action="store_false",
           help="sleep during calculations")
    parser.add_option("-s", "--silent",  dest="verbose", action="store_false", 
           help="not activated yet")
    parser.add_option("--verbosity",  dest="verbosity", action="store_const", 
           help="level of printout - not activated yet")
    parser.add_option("-o", "--pH", dest="pH", type="float", default=7.0, 
           help="setting pH-value used in e.g. stability calculations [7.0]")
    parser.add_option("-w", "--window", dest="window", nargs=3, type="float", default=(0.0, 14.0, 1.0),
           help="setting the pH-window to show e.g. stability profiles [0.0, 14.0, 1.0]")
    parser.add_option("-g", "--grid",   dest="grid",   nargs=3, type="float", default=(0.0, 14.0, 0.1),
           help="setting the pH-grid to calculate e.g. stability related properties [0.0, 14.0, 0.1]")
    parser.add_option("--mutator", dest="mutator", 
           help="setting approach for mutating <filename> [alignment/scwrl/jackal]")
    parser.add_option("--mutator-option", dest="mutator_options", action="append",
           help="setting property for mutator [e.g. type=\"side-chain\"]")

    parser.add_option("-d","--display-coupled-residues", dest="display_coupled_residues", action="store_true",
           help="Displays alternative pKa values due to coupling of titratable groups")
    parser.add_option("-l","--reuse-ligand-mol2-files", dest="reuse_ligand_mol2_file", action="store_true",
           help="Reuses the ligand mol2 files allowing the user to alter ligand bond orders", default=False)
    parser.add_option("-k","--keep-protons", dest="keep_protons", action="store_true",
           help="Keep protons in input file", default=False)
    parser.add_option("--protonate-all", dest="protonate_all", action="store_true",
           help="Protonate all atoms (will not influence pKa calculation)", default=False)


    # parsing and returning options and arguments
    options, args = parser.parse_args()

    # adding specified filenames to arguments
    if options.filenames:
      for filename in options.filenames:
        args.append(filename)

    # checking at early stage that there is at least one pdbfile to work with
    if len(args) == 0:
      print("Warning: no pdbfile provided")
      #sys.exit(9)



    # done!
    return options, args





def makeTidyAtomLabel(name,element):
    """
    Returns a 'tidier' atom label for printing the new pdbfile
    """
    
    if len(name)>4:# if longer than 4, just truncate the name
        label=name[0:4]
    elif len(name)==4:# if lenght is 4, otherwise use the name as it is
        label = name
    else: # if less than 4 characters long, insert white space as needed
        if len(element)==1:
            label = ' %-3s'%name
        else: # The element shoul occupy the two first chars
            label = '%-4s'%name

    return label



def get_sorted_configurations(configuration_keys):
    """
    extract and sort configurations
    """
    configurations = list(configuration_keys)
    configurations.sort(key=configuration_compare)
    return configurations

def configuration_compare(conf):
    return 100*int(conf[1:-2]) + ord(conf[-1])




def writeFile(filename, lines):
    """
    Writes a new file
    """
    file = open(filename, 'w')

    for line in lines:
        file.write( "%s\n" % (line) )
    file.close()


