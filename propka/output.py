
from __future__ import division
from __future__ import print_function

import sys

import propka.lib
from propka.lib import info, warning


def printHeader():
    """
    prints the header section
    """
    str  = "%s\n" % ( getPropkaHeader() )
    str += "%s\n" % ( getReferencesHeader() )
    str += "%s\n" % ( getWarningHeader() )

    info(str)


def writePDB(protein, file=None, filename=None, include_hydrogens=False, options=None):
    """
    Write the residue to the new pdbfile
    """

    if file == None:
      # opening file if not given
      if filename == None:
        filename = "%s.pdb" % (protein.name)
      file = open(filename, 'w')
      info("writing pdbfile %s" % (filename))
      close_file = True
    else:
      # don't close the file, it was opened in a different place
      close_file = False

    numb = 0
    for chain in protein.chains:
      for residue in chain.residues:
        if residue.resName not in ["N+ ", "C- "]:
          for atom in residue.atoms:
            if include_hydrogens == False and atom.name[0] == "H":
              """ don't print """
            else:
              numb += 1
              line = atom.makePDBLine(numb=numb)
              line += "\n"
              file.write(line)

    if close_file == True:
      file.close()


def writePKA(protein, parameters, filename=None, conformation ='1A',reference="neutral", direction="folding", verbose=False, options=None):
    """
    Write the pka-file based on the given protein
    """
    verbose = True
    if filename == None:
      filename = "%s.pka" % (protein.name)
    file = open(filename, 'w')
    if verbose == True:
      info("Writing %s" % (filename))

    # writing propka header
    str  = "%s\n" % ( getPropkaHeader() )
    str += "%s\n" % ( getReferencesHeader() )
    str += "%s\n" % ( getWarningHeader() )

    # writing pKa determinant section
    str += getDeterminantSection(protein,conformation, parameters)

    # writing pKa summary section
    str += getSummarySection(protein,conformation,parameters)
    str += "%s\n" % ( getTheLine() )

    # printing Folding Profile
    str += getFoldingProfileSection(protein, conformation=conformation, reference=reference, direction=direction, window=[0., 14., 1.0], options=options)

    # printing Protein Charge Profile
    str += getChargeProfileSection(protein, conformation=conformation)

    # now, writing the pka text to file
    file.write(str)

    file.close()


def printTmProfile(protein, reference="neutral", window=[0., 14., 1.], Tm=[0.,0.], Tms=None, ref=None, verbose=False, options=None):
    """
    prints Tm profile
    """
    profile = protein.getTmProfile(reference=reference, grid=[0., 14., 0.1], Tms=Tms, ref=ref, options=options)
    if profile == None:
      str  = "Could not determine Tm-profile\n"
    else:
      str  = " suggested Tm-profile for %s\n" % (protein.name)
      for (pH, Tm) in profile:
        if pH >= window[0] and pH <= window[1] and (pH%window[2] < 0.01 or pH%window[2] > 0.99*window[2]):
          str += "%6.2lf%10.2lf\n" % (pH, Tm)
      info(str)


def printResult(protein, conformation, parameters):
    """
    prints all resulting output from determinants and down
    """
    printPKASection(protein, conformation, parameters)


def printPKASection(protein, conformation, parameters):
    """
    prints out the pka-section of the result
    """
    # geting the determinants section
    str = getDeterminantSection(protein, conformation, parameters)
    info(str)

    str = getSummarySection(protein,conformation,parameters)
    info(str)


def getDeterminantSection(protein, conformation, parameters):
    """
    prints out the pka-section of the result
    """
    # getting the same order as in propka2.0
    str  = "%s\n" % ( getDeterminantsHeader() )
    # printing determinants
    for chain in protein.conformations[conformation].chains:
        for residue_type in parameters.write_out_order:
            groups = [g for g in protein.conformations[conformation].groups if g.atom.chainID == chain]
            for group in groups:
                if group.residue_type == residue_type:
                    str += "%s" % ( group.getDeterminantString(parameters.remove_penalised_group) )

    # Add a warning in case of coupled residues
    if protein.conformations[conformation].non_covalently_coupled_groups and not protein.options.display_coupled_residues:
        str += 'Coupled residues (marked *) were detected. Please rerun PropKa with the --display-coupled-residues \nor -d option for detailed information.\n'

    return str


def getSummarySection(protein, conformation, parameters):
    """
    prints out the pka-section of the result
    """
    str  = "%s\n" % ( getSummaryHeader() )
    # printing pKa summary
    for residue_type in parameters.write_out_order:
        for group in protein.conformations[conformation].groups:
          if group.residue_type == residue_type:
            str += "%s" % ( group.getSummaryString(parameters.remove_penalised_group) )

    return str


def getFoldingProfileSection(protein, conformation='AVR', direction="folding", reference="neutral", window=[0., 14., 1.0], verbose=False, options=None):
    """
    returns the protein-folding-profile section
    """
    str  = getTheLine()
    str += "\n"
    str += "Free energy of %9s (kcal/mol) as a function of pH (using %s reference)\n" % (direction, reference)

    profile, [pH_opt, dG_opt], [dG_min, dG_max], [pH_min, pH_max] = protein.getFoldingProfile(conformation=conformation,
                                                                                 reference=reference,
                                                                                 direction=direction, grid=[0., 14., 0.1], options=options)
    if profile == None:
      str += "Could not determine folding profile\n"
    else:
      for (pH, dG) in profile:
        if pH >= window[0] and pH <= window[1] and (pH%window[2] < 0.05 or pH%window[2] > 0.95):
          str += "%6.2lf%10.2lf\n" % (pH, dG)
      str += "\n"

    if pH_opt == None or dG_opt == None:
      str += "Could not determine pH optimum\n"
    else:
      str += "The pH of optimum stability is %4.1lf for which the free energy is%6.1lf kcal/mol at 298K\n" % (pH_opt, dG_opt)

    if dG_min == None or dG_max == None:
      str += "Could not determine pH values where the free energy is within 80 %s of minimum\n" % ("%")
    else:
      str += "The free energy is within 80 %s of maximum at pH %4.1lf to %4.1lf\n" % ("%", dG_min, dG_max)

    if pH_min == None or pH_max == None:
      str += "Could not determine the pH-range where the free energy is negative\n\n"
    else:
      str += "The free energy is negative in the range %4.1lf - %4.1lf\n\n" % (pH_min, pH_max)


    return str



def getChargeProfileSection(protein, conformation='AVR', options=None):
    """
    returns the protein-folding-profile section
    """
    str  = "Protein charge of folded and unfolded state as a function of pH\n"

    profile = protein.getChargeProfile(conformation=conformation,grid=[0., 14., 1.])
    if profile == None:
      str += "Could not determine charge profile\n"
    else:
      str += "%6s%10s%8s\n" % ("pH", "unfolded", "folded")
      for (pH, Q_mod, Q_pro) in profile:
        str += "%6.2lf%10.2lf%8.2lf\n" % (pH, Q_mod, Q_pro)


    pI_pro, pI_mod = protein.getPI(conformation=conformation)
    if pI_pro == None or pI_mod == None:
      str += "Could not determine the pI\n\n"
    else:
      str += "The pI is %5.2lf (folded) and %5.2lf (unfolded)\n" % (pI_pro, pI_mod)


    return str


def writeJackalScapFile(mutationData=None, filename="1xxx_scap.list", options=None):
    """
    writing a scap file for, i.e., generating a mutated protein
    """
    file = open(filename, 'w')

    for chainID, code1, resNumb, code2 in mutationData:
      str = "%s, %d, %s\n" % (chainID, resNumb, code2)
      file.write(str)
    file.close()


def writeScwrlSequenceFile(sequence, filename="x-ray.seq", options=None):
    """
    writing a scwrl sequence file for, e.g.,  generating a mutated protein
    """
    file = open(filename, 'w')

    start = 0
    while len(sequence[start:]) > 60:
      file.write( "%s\n" % (sequence[start:start+60]) )
      start += 60
    file.write( "%s\n" % (sequence[start:]) )

    file.close()



# --- various header text --- #


def getPropkaHeader():
    """
    Creates the header
    """
    from datetime import date
    today = date.today()



    str  = "propka3.1 %93s\n" % (today)
    str += "-------------------------------------------------------------------------------------------------------\n"
    str += "--                                                                                                   --\n"
    str += "--                                   PROPKA: A PROTEIN PKA PREDICTOR                                 --\n"
    str += "--                                                                                                   --\n"
    str += "--                                 VERSION 1.0,  04/25/2004, IOWA CITY                               --\n"
    str += "--                                             BY HUI LI                                             --\n"
    str += "--                                                                                                   --\n"
    str += "--                            VERSION 2.0,  11/05/2007, IOWA CITY/COPENHAGEN                         --\n"
    str += "--                                BY DELPHINE C. BAS AND DAVID M. ROGERS                             --\n"
    str += "--                                                                                                   --\n"
    str += "--                                VERSION 3.0,  01/06/2011, COPENHAGEN                               --\n"
    str += "--                            BY MATS H.M. OLSSON AND CHRESTEN R. SONDERGARD                         --\n"
    str += "--                                                                                                   --\n"
    str += "--                                VERSION 3.1,  07/01/2011, COPENHAGEN                               --\n"
    str += "--                            BY CHRESTEN R. SONDERGARD AND MATS H.M. OLSSON                         --\n"
    str += "-------------------------------------------------------------------------------------------------------\n"
    str += "\n"



    return str


def getReferencesHeader():
    """
    Returns the 'references' part in output file
    """

    str  = ""
    str += "-------------------------------------------------------------------------------------------------------\n"
    str += " References:\n"
    str += "\n"
    str += "   Very Fast Empirical Prediction and Rationalization of Protein pKa Values\n"
    str += "   Hui Li, Andrew D. Robertson and Jan H. Jensen\n"
    str += "   PROTEINS: Structure, Function, and Bioinformatics 61:704-721 (2005)\n"
    str += "   \n"
    str += "   Very Fast Prediction and Rationalization of pKa Values for Protein-Ligand Complexes\n"
    str += "   Delphine C. Bas, David M. Rogers and Jan H. Jensen\n"
    str += "   PROTEINS: Structure, Function, and Bioinformatics 73:765-783 (2008)\n"
    str += "   \n"
    str += "   PROPKA3: Consistent Treatment of Internal and Surface Residues in Empirical pKa predictions\n"
    str += "   Mats H.M. Olsson, Chresten R. Sondergard, Michal Rostkowski, and Jan H. Jensen\n"
    str += "   Journal of Chemical Theory and Computation, 7(2):525-537 (2011)\n"
    str += "   \n"
    str += "   Improved Treatment of Ligands and Coupling Effects in Empirical Calculation\n"
    str += "    and Rationalization of pKa Values\n"
    str += "   Chresten R. Sondergaard, Mats H.M. Olsson, Michal Rostkowski, and Jan H. Jensen\n"
    str += "   Journal of Chemical Theory and Computation, (2011)\n"
    str += "   \n"
    str += "-------------------------------------------------------------------------------------------------------\n"

    return str


def getWarningHeader():
    """
    Returns the 'warning' part in output file
    """

    str  = ""

    return str


def getDeterminantsHeader():
    """
    Creates the Determinant header
    """
    str  = ""
    str += "---------  -----   ------   ---------------------    --------------    --------------    --------------\n"
    str += "                            DESOLVATION  EFFECTS       SIDECHAIN          BACKBONE        COULOMBIC    \n"
    str += " RESIDUE    pKa    BURIED     REGULAR      RE        HYDROGEN BOND     HYDROGEN BOND      INTERACTION  \n"
    str += "---------  -----   ------   ---------   ---------    --------------    --------------    --------------\n"

    return str


def getSummaryHeader():
    """
    returns the summary header
    """
    str  = getTheLine()
    str += "\n"
    str += "SUMMARY OF THIS PREDICTION\n"
    str += "       Group      pKa  model-pKa   ligand atom-type"

    return str


def getTheLine():
    """
    draw the line - Johnny Cash would have been proud - or actually Aerosmith!
    """
    str  = ""
    for i in range(0, 104):
      str += "-"

    return str


# Interaction maps
def make_interaction_map(name, list, interaction):
    """ Print out an interaction map named 'name' of the groups in 'list'
    based on the function 'interaction' """

    # return an empty string, if the list is empty
    if len(list)==0:
        return ''

    # for long list, use condensed formatting
    if len(list)>10:
        res = 'Condensed form:\n'
        for i in range(len(list)):
            for j in range(i,len(list)):
                if interaction(list[i],list[j]):
                    res += 'Coupling: %9s - %9s\n'%(list[i].label,list[j].label)
        return res

    # Name and map header
    res = '%s\n%12s'%(name,'')
    for g in list:
        res += '%9s | '%g.label

    # do the map
    for g1 in list:
        res += '\n%-12s'%(g1.label)
        for g2 in list:
            tag = ''
            if interaction(g1, g2):
                tag = '    X     '
            res += '%10s| '%tag

    return res

