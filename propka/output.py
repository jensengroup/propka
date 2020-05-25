"""Output routines."""
from datetime import date
from propka.lib import info


def print_header():
    """Print header section of output."""
    str_ = "%s\n" % get_propka_header()
    str_ += "%s\n" % get_references_header()
    str_ += "%s\n" % get_warning_header()
    info(str_)


def write_pdb(protein, pdbfile=None, filename=None, include_hydrogens=False,
              _=None):
    """Write a residue to the new PDB file.

    Args:
        protein:  protein object
        pdbfile:  PDB file
        filename:  file to write to
        include_hydrogens:  Boolean indicating whether to include hydrogens
        options:  options object
    """
    if pdbfile is None:
        # opening file if not given
        if filename is None:
            filename = "%s.pdb" % (protein.name)
        # TODO - this would be better as a context manager
        pdbfile = open(filename, 'w')
        info("writing pdbfile %s" % (filename))
        close_file = True
    else:
        # don't close the file, it was opened in a different place
        close_file = False
    numb = 0
    for chain in protein.chains:
        for residue in chain.residues:
            if residue.res_name not in ["N+ ", "C- "]:
                for atom in residue.atoms:
                    if (not include_hydrogens) and atom.name[0] == "H":
                        # don't print
                        pass
                    else:
                        numb += 1
                        line = atom.make_pdb_line2(numb=numb)
                        line += "\n"
                        pdbfile.write(line)
    if close_file:
        pdbfile.close()


def write_pka(protein, parameters, filename=None, conformation='1A',
              reference="neutral", _="folding", verbose=False,
              __=None):
    """Write the pKa-file based on the given protein.

    Args:
        protein:  protein object
        filename:  output file name
        conformation:  TODO - figure this out
        reference:  reference state
        _:  "folding" or other
        verbose:  Boolean flag for verbosity
        __:  options object
    """
    # TODO - the code immediately overrides the verbose argument; why?
    verbose = True
    if filename is None:
        filename = "%s.pka" % (protein.name)
    # TODO - this would be much better with a context manager
    file_ = open(filename, 'w')
    if verbose:
        info("Writing %s" % (filename))
    # writing propka header
    str_ = "%s\n" % get_propka_header()
    str_ += "%s\n" % get_references_header()
    str_ += "%s\n" % get_warning_header()
    # writing pKa determinant section
    str_ += get_determinant_section(protein, conformation, parameters)
    # writing pKa summary section
    str_ += get_summary_section(protein, conformation, parameters)
    str_ += "%s\n" % get_the_line()
    # printing Folding Profile
    str_ += get_folding_profile_section(protein, conformation=conformation,
                                        reference=reference,
                                        window=[0., 14., 1.0])
    # printing Protein Charge Profile
    str_ += get_charge_profile_section(protein, conformation=conformation)
    # now, writing the pka text to file
    file_.write(str_)
    file_.close()


def print_tm_profile(protein, reference="neutral", window=[0., 14., 1.],
                     __=[0., 0.], tms=None, ref=None, _=False,
                     options=None):
    """Print Tm profile.

    I think Tm refers to the denaturation temperature.

    Args:
        protein:  protein object
        reference:  reference state
        window:  pH window [min, max, step]
        __:  temperature range [min, max]
        tms:  TODO - figure this out
        ref:  TODO - figure this out (probably reference state?)
        _:  Boolean for verbosity
        options:  options object
    """
    profile = protein.getTmProfile(reference=reference, grid=[0., 14., 0.1],
                                   tms=tms, ref=ref, options=options)
    if profile is None:
        str_ = "Could not determine Tm-profile\n"
    else:
        str_ = " suggested Tm-profile for %s\n" % (protein.name)
        for (ph, tm_) in profile:
            if ph >= window[0] and ph <= window[1] and (ph%window[2] < 0.01 \
                or ph%window[2] > 0.99*window[2]):
                str_ += "%6.2lf%10.2lf\n" % (ph, tm_)
        info(str_)


def print_result(protein, conformation, parameters):
    """Prints all resulting output from determinants and down.

    Args:
        protein:  protein object
        conformation:  specific conformation
        parameters:  parameters
    """
    print_pka_section(protein, conformation, parameters)


def print_pka_section(protein, conformation, parameters):
    """Prints out pKa section of results.

    Args:
        protein:  protein object
        conformation:  specific conformation
        parameters:  parameters
    """
    # geting the determinants section
    str_ = get_determinant_section(protein, conformation, parameters)
    info(str_)
    str_ = get_summary_section(protein, conformation, parameters)
    info(str_)


def get_determinant_section(protein, conformation, parameters):
    """Returns string with determinant section of results.

    Args:
        protein:  protein object
        conformation:  specific conformation
        parameters:  parameters
    Returns:
        string
    """
    # getting the same order as in propka2.0
    str_ = "%s\n" % get_determinants_header()
    # printing determinants
    for chain in protein.conformations[conformation].chains:
        for residue_type in parameters.write_out_order:
            groups = [g for g in protein.conformations[conformation].groups \
                if g.atom.chain_id == chain]
            for group in groups:
                if group.residue_type == residue_type:
                    str_ += "%s" \
                        % group.get_determinant_string(parameters.remove_penalised_group)
    # Add a warning in case of coupled residues
    if protein.conformations[conformation].non_covalently_coupled_groups \
        and not protein.options.display_coupled_residues:
        str_ += 'Coupled residues (marked *) were detected.'
        str_ += 'Please rerun PropKa with the --display-coupled-residues \n'
        str_ += 'or -d option for detailed information.\n'
    return str_


def get_summary_section(protein, conformation, parameters):
    """Returns string with summary section of the results.

    Args:
        protein:  protein object
        conformation:  specific conformation
        parameters:  parameters
    Returns:
        string
    """
    str_ = "%s\n" % get_summary_header()
    # printing pKa summary
    for residue_type in parameters.write_out_order:
        for group in protein.conformations[conformation].groups:
            if group.residue_type == residue_type:
                str_ += "%s" \
                    % group.get_summary_string(parameters.remove_penalised_group)
    return str_


def get_folding_profile_section(protein, conformation='AVR',
                                direction="folding", reference="neutral",
                                window=[0., 14., 1.0], _=False,
                                __=None):
    """Returns string with the folding profile section of the results.

    Args:
        protein:  protein object
        conformation:  specific conformation
        direction:  'folding' or other
        reference:  reference state
        window:  pH window [min, max, step]
        _:  Boolean for verbose output
        __:  options object
    Returns:
        string
    """
    str_ = get_the_line()
    str_ += "\n"
    str_ += "Free energy of %9s (kcal/mol) as a function" % direction
    str_ += " of pH (using %s reference)\n" %  reference
    profile, [ph_opt, dg_opt], [dg_min, dg_max], [ph_min, ph_max] \
        = protein.get_folding_profile(conformation=conformation,
                                      reference=reference, grid=[0., 14., 0.1])
    if profile is None:
        str_ += "Could not determine folding profile\n"
    else:
        for (ph, dg) in profile:
            if ph >= window[0] and ph <= window[1]:
                if ph%window[2] < 0.05 or ph%window[2] > 0.95:
                    str_ += "%6.2lf%10.2lf\n" % (ph, dg)
        str_ += "\n"
    if ph_opt is None or dg_opt is None:
        str_ += "Could not determine pH optimum\n"
    else:
        str_ += "The pH of optimum stability is %4.1lf" % ph_opt
        str_ += " for which the free energy is %6.1lf kcal/mol at 298K\n" % dg_opt
    if dg_min is None or dg_max is None:
        str_ += "Could not determine pH values where the free energy"
        str_ += " is within 80 %s of minimum\n" % ("%")
    else:
        str_ += "The free energy is within 80 \% of maximum"
        str_ += " at pH %4.1lf to %4.1lf\n" % (dg_min, dg_max)
    if ph_min is None or ph_max is None:
        str_ += "Could not determine the pH-range where the free"
        str_ += " energy is negative\n\n"
    else:
        str_ += "The free energy is negative in the range"
        str_ += " %4.1lf - %4.1lf\n\n" % (ph_min, ph_max)
    return str_


def get_charge_profile_section(protein, conformation='AVR', _=None):
    """Returns string with the charge profile section of the results.

    Args:
        protein:  protein object
        conformation:  specific conformation
        _:  options object
    Returns:
        string
    """
    str_ = "Protein charge of folded and unfolded state as a function of pH\n"
    profile = protein.get_charge_profile(conformation=conformation,
                                         grid=[0., 14., 1.])
    if profile is None:
        str_ += "Could not determine charge profile\n"
    else:
        str_ += "%6s%10s%8s\n" % ("pH", "unfolded", "folded")
        for (ph, q_mod, q_pro) in profile:
            str_ += "%6.2lf%10.2lf%8.2lf\n" % (ph, q_mod, q_pro)
    pi_pro, pi_mod = protein.get_pi(conformation=conformation)
    if pi_pro is None or pi_mod is None:
        str_ += "Could not determine the pI\n\n"
    else:
        str_ += "The pI is %5.2lf (folded) and %5.2lf (unfolded)\n" % (pi_pro,
                                                                       pi_mod)
    return str_


def write_jackal_scap_file(mutation_data=None, filename="1xxx_scap.list",
                           _=None):
    """Write a scap file for, i.e., generating a mutated protein

    TODO - figure out what this is
    """
    with open(filename, 'w') as file_:
        for chain_id, _, res_num, code2 in mutation_data:
            str_ = "%s, %d, %s\n" % (chain_id, res_num, code2)
        file_.write(str_)


def write_scwrl_sequence_file(sequence, filename="x-ray.seq", _=None):
    """Write a scwrl sequence file for, e.g.,  generating a mutated protein

    TODO - figure out what this is
    """
    with open(filename, 'w') as file_:
        start = 0
        while len(sequence[start:]) > 60:
            file_.write("%s\n" % (sequence[start:start+60]))
            start += 60
        file_.write("%s\n" % (sequence[start:]))


def get_propka_header():
    """Create the header.

    Returns:
        string
    """
    today = date.today()
    str_ = "propka3.1 %93s\n" % (today)
    str_ += "-------------------------------------------------------------------------------------------------------\n"
    str_ += "--                                                                                                   --\n"
    str_ += "--                                   PROPKA: A PROTEIN PKA PREDICTOR                                 --\n"
    str_ += "--                                                                                                   --\n"
    str_ += "--                                 VERSION 1.0,  04/25/2004, IOWA CITY                               --\n"
    str_ += "--                                             BY HUI LI                                             --\n"
    str_ += "--                                                                                                   --\n"
    str_ += "--                            VERSION 2.0,  11/05/2007, IOWA CITY/COPENHAGEN                         --\n"
    str_ += "--                                BY DELPHINE C. BAS AND DAVID M. ROGERS                             --\n"
    str_ += "--                                                                                                   --\n"
    str_ += "--                                VERSION 3.0,  01/06/2011, COPENHAGEN                               --\n"
    str_ += "--                            BY MATS H.M. OLSSON AND CHRESTEN R. SONDERGARD                         --\n"
    str_ += "--                                                                                                   --\n"
    str_ += "--                                VERSION 3.1,  07/01/2011, COPENHAGEN                               --\n"
    str_ += "--                            BY CHRESTEN R. SONDERGARD AND MATS H.M. OLSSON                         --\n"
    str_ += "-------------------------------------------------------------------------------------------------------\n"
    str_ += "\n"
    return str_


def get_references_header():
    """Create the 'references' part of output file.

    Returns:
        string
    """
    str_ = ""
    str_ += "-------------------------------------------------------------------------------------------------------\n"
    str_ += " References:\n"
    str_ += "\n"
    str_ += "   Very Fast Empirical Prediction and Rationalization of Protein pKa Values\n"
    str_ += "   Hui Li, Andrew D. Robertson and Jan H. Jensen\n"
    str_ += "   PROTEINS: Structure, Function, and Bioinformatics 61:704-721 (2005)\n"
    str_ += "   \n"
    str_ += "   Very Fast Prediction and Rationalization of pKa Values for Protein-Ligand Complexes\n"
    str_ += "   Delphine C. Bas, David M. Rogers and Jan H. Jensen\n"
    str_ += "   PROTEINS: Structure, Function, and Bioinformatics 73:765-783 (2008)\n"
    str_ += "   \n"
    str_ += "   PROPKA3: Consistent Treatment of Internal and Surface Residues in Empirical pKa predictions\n"
    str_ += "   Mats H.M. Olsson, Chresten R. Sondergard, Michal Rostkowski, and Jan H. Jensen\n"
    str_ += "   Journal of Chemical Theory and Computation, 7(2):525-537 (2011)\n"
    str_ += "   \n"
    str_ += "   Improved Treatment of Ligands and Coupling Effects in Empirical Calculation\n"
    str_ += "    and Rationalization of pKa Values\n"
    str_ += "   Chresten R. Sondergaard, Mats H.M. Olsson, Michal Rostkowski, and Jan H. Jensen\n"
    str_ += "   Journal of Chemical Theory and Computation, (2011)\n"
    str_ += "   \n"
    str_ += "-------------------------------------------------------------------------------------------------------\n"
    return str_


def get_warning_header():
    """Create the 'warning' part of the output file.

    TODO - this function is essentially a no-op.

    Returns:
        string
    """
    str_ = ""
    return str_


def get_determinants_header():
    """Create the Determinant header.

    Returns:
        string
    """
    str_ = ""
    str_ += "---------  -----   ------   ---------------------    --------------    --------------    --------------\n"
    str_ += "                            DESOLVATION  EFFECTS       SIDECHAIN          BACKBONE        COULOMBIC    \n"
    str_ += " RESIDUE    pKa    BURIED     REGULAR      RE        HYDROGEN BOND     HYDROGEN BOND      INTERACTION  \n"
    str_ += "---------  -----   ------   ---------   ---------    --------------    --------------    --------------\n"
    return str_


def get_summary_header():
    """Create the summary header.

    Returns:
        string
    """
    str_ = get_the_line()
    str_ += "\n"
    str_ += "SUMMARY OF THIS PREDICTION\n"
    str_ += "       Group      pKa  model-pKa   ligand atom-type"
    return str_


def get_the_line():
    """Draw the line - Johnny Cash would have been proud - or actually Aerosmith!

    NOTE - Johnny Cash walked the line.

    Returns:
        string
    """
    str_ = ""
    str_ += ("-" * 104)
    return str_


def make_interaction_map(name, list_, interaction):
    """Print out an interaction map named 'name' of the groups in 'list'
    based on the function 'interaction'

    Args:
        list_:  list of groups
        interaction:  some sort of function
    Returns:
        string
    """
    # return an empty string, if the list is empty
    if len(list_) == 0:
        return ''
    # for long list, use condensed formatting
    if len(list_) > 10:
        res = 'Condensed form:\n'
        for i, group1 in enumerate(list_):
            for group2 in list_[i:]:
                if interaction(group1, group2):
                    res += 'Coupling: %9s - %9s\n' % (group1.label,
                                                      group2.label)
        return res
    # Name and map header
    res = '%s\n%12s' % (name, '')
    for group in list_:
        res += '%9s | ' % group.label
    # do the map
    for group1 in list_:
        res += '\n%-12s' % (group1.label)
        for group2 in list_:
            tag = ''
            if interaction(group1, group2):
                tag = '    X     '
            res += '%10s| '%tag
    return res
