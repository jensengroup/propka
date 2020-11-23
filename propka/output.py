"""
Output
======

Output routines.
"""
from datetime import date
from propka.lib import info
from . import __version__


def open_file_for_writing(input_file):
    """Open file or file-like stream for writing.

    TODO - convert this to a context manager.

    Args:
        input_file: path to file or file-like object. If file-like object,
        then will attempt to get file mode.
    """
    try:
        if not input_file.writable():
            raise IOError("File/stream not open for writing")
        return input_file
    except AttributeError:
        pass
    try:
        file_ = open(input_file, "wt")
    except FileNotFoundError:
        raise Exception(f"Could not open {input_file}")
    return file_


def write_file(filename, lines):
    """Writes a new file.

    Args:
        filename:  name of file
        lines:  lines to write to file
    """
    file_ = open_file_for_writing(filename)
    for line in lines:
        file_.write(f"{line}\n")
    file_.close()


def print_header():
    """Print header section of output."""
    info(
        f"{get_propka_header()}\n"
        f"{get_references_header()}\n"
        f"{get_warning_header()}\n"
    )


def write_pdb_for_protein(
    protein, pdbfile=None, filename=None, include_hydrogens=False, _=None
):
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
            filename = f"{protein.name}.pdb"
        # TODO - this would be better as a context manager
        pdbfile = open(filename, "w")
        info(f"writing pdbfile {filename}")
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


def write_pdb_for_conformation(conformation, filename):
    """Write PDB conformation to a file.

    Args:
        conformation:  conformation container
        filename:  filename for output
    """
    write_pdb_for_atoms(conformation.atoms, filename)


def write_pka(
    protein,
    parameters,
    filename=None,
    conformation="1A",
    reference="neutral",
    _="folding",
    verbose=False,
    __=None,
):
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
        filename = f"{protein.name}.pka"
    # TODO - this would be much better with a context manager
    file_ = open(filename, "w")
    if verbose:
        info(f"Writing {filename}")
    # writing propka header
    str_ = f"{get_propka_header()}\n"
    str_ += f"{get_references_header()}\n"
    str_ += f"{get_warning_header()}\n"
    # writing pKa determinant section
    str_ += get_determinant_section(protein, conformation, parameters)
    # writing pKa summary section
    str_ += get_summary_section(protein, conformation, parameters)
    str_ += f"{get_the_line()}\n"
    # printing Folding Profile
    str_ += get_folding_profile_section(
        protein,
        conformation=conformation,
        reference=reference,
        window=[0.0, 14.0, 1.0],
    )
    # printing Protein Charge Profile
    str_ += get_charge_profile_section(protein, conformation=conformation)
    # now, writing the pka text to file
    file_.write(str_)
    file_.close()


def print_tm_profile(
    protein,
    reference="neutral",
    window=[0.0, 14.0, 1.0],
    __=[0.0, 0.0],
    tms=None,
    ref=None,
    _=False,
    options=None,
):
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
    profile = protein.getTmProfile(
        reference=reference,
        grid=[0.0, 14.0, 0.1],
        tms=tms,
        ref=ref,
        options=options,
    )
    if profile is None:
        str_ = "Could not determine Tm-profile\n"
    else:
        str_ = f" suggested Tm-profile for {protein.name}\n"
        for (ph, tm_) in profile:
            if (
                ph >= window[0]
                and ph <= window[1]
                and (
                    ph % window[2] < 0.01 or ph % window[2] > 0.99 * window[2]
                )
            ):
                str_ += f"{ph:>6.2f}{tm_:>10.2f}\n"
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
    str_ = f"{get_determinants_header()}\n"
    # printing determinants
    for chain in protein.conformations[conformation].chains:
        for residue_type in parameters.write_out_order:
            groups = [
                g
                for g in protein.conformations[conformation].groups
                if g.atom.chain_id == chain
            ]
            for group in groups:
                if group.residue_type == residue_type:
                    v = group.get_determinant_string(
                        parameters.remove_penalised_group
                    )
                    str_ += f"{v}"
    # Add a warning in case of coupled residues
    if (
        protein.conformations[conformation].non_covalently_coupled_groups
        and not protein.options.display_coupled_residues
    ):
        str_ += "Coupled residues (marked *) were detected."
        str_ += "Please rerun PropKa with the --display-coupled-residues \n"
        str_ += "or -d option for detailed information.\n"
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
    str_ = f"{get_summary_header()}\n"
    # printing pKa summary
    for residue_type in parameters.write_out_order:
        for group in protein.conformations[conformation].groups:
            if group.residue_type == residue_type:
                v = group.get_summary_string(parameters.remove_penalised_group)
                str_ += f"{v}"
    return str_


def get_folding_profile_section(
    protein,
    conformation="AVR",
    direction="folding",
    reference="neutral",
    window=[0.0, 14.0, 1.0],
    _=False,
    __=None,
):
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
    str_ += f"Free energy of {direction:>9s} (kcal/mol) as a function"
    str_ += f" of pH (using {reference} reference)\n"
    (
        profile,
        [ph_opt, dg_opt],
        [dg_min, dg_max],
        [ph_min, ph_max],
    ) = protein.get_folding_profile(
        conformation=conformation, reference=reference, grid=[0.0, 14.0, 0.1]
    )
    if profile is None:
        str_ += "Could not determine folding profile\n"
    else:
        for (ph, dg) in profile:
            if (
                ph >= window[0]
                and ph <= window[1]
                and (ph % window[2] < 0.05 or ph % window[2] > 0.95)
            ):
                str_ += f"{ph:>6.2f}{dg:>10.2f}\n"
        str_ += "\n"
    if ph_opt is None or dg_opt is None:
        str_ += "Could not determine pH optimum\n"
    else:
        str_ += f"The pH of optimum stability is {ph_opt:>4.1f}"
        str_ += f" for which the free energy is {dg_opt:>6.1f} "
        str_ += "kcal/mol at 298K\n"
    if dg_min is None or dg_max is None:
        str_ += "Could not determine pH values where the free energy"
        str_ += " is within 80 % of minimum\n"
    else:
        str_ += "The free energy is within 80 % of maximum"
        str_ += f" at pH {dg_min:>4.1f} to {dg_max:>4.1f}\n"
    if ph_min is None or ph_max is None:
        str_ += "Could not determine the pH-range where the free"
        str_ += " energy is negative\n\n"
    else:
        str_ += "The free energy is negative in the range"
        str_ += f" {ph_min:>4.1f} - {ph_max:>4.1f}\n\n"
    return str_


def get_charge_profile_section(protein, conformation="AVR", _=None):
    """Returns string with the charge profile section of the results.

    Args:
        protein:  protein object
        conformation:  specific conformation
        _:  options object
    Returns:
        string
    """
    str_ = "Protein charge of folded and unfolded state as a function of pH\n"
    profile = protein.get_charge_profile(
        conformation=conformation, grid=[0.0, 14.0, 1.0]
    )
    if profile is None:
        str_ += "Could not determine charge profile\n"
    else:
        str_ += "    pH  unfolded  folded\n"
        for (ph, q_mod, q_pro) in profile:
            str_ += f"{ph:6.2f}{q_mod:10.2f}{q_pro:8.2f}\n"
    pi_pro, pi_mod = protein.get_pi(conformation=conformation)
    if pi_pro is None or pi_mod is None:
        str_ += "Could not determine the pI\n\n"
    else:
        str_ += f"The pI is {pi_pro:>5.2f} (folded) and {pi_mod:>5.2f} "
        str_ += "(unfolded)\n"
    return str_


def write_jackal_scap_file(
    mutation_data=None, filename="1xxx_scap.list", _=None
):
    """Write a scap file for, i.e., generating a mutated protein

    TODO - figure out what this is
    """
    with open(filename, "w") as file_:
        for chain_id, _, res_num, code2 in mutation_data:
            file_.write(f"{chain_id}, {res_num:d}, {code2}\n")


def write_scwrl_sequence_file(sequence, filename="x-ray.seq", _=None):
    """Write a scwrl sequence file for, e.g.,  generating a mutated protein

    TODO - figure out what this is
    """
    with open(filename, "w") as file_:
        start = 0
        while len(sequence[start:]) > 60:
            file_.write(f"{sequence[start : start + 60]}s\n")
            start += 60
        file_.write(f"{sequence[start:]}\n")


def get_propka_header():
    """Create the header.

    Returns:
        string
    """
    today = date.today()
    str_ = f"propka{__version__:<53} {today!s:>43}\n"
    str_ += """
-------------------------------------------------------------------------------
--                                                                           --
--  PROPKA: A PROTEIN PKA PREDICTOR                                          --
--                                                                           --
--  VERSION 1.0,  04/25/2004,  IOWA CITY                                     --
--  BY HUI LI                                                                --
--                                                                           --
--  VERSION 2.0,  11/05/2007, IOWA CITY/COPENHAGEN                           --
--  BY DELPHINE C. BAS AND DAVID M. ROGERS                                   --
--                                                                           --
--  VERSION 3.0,  01/06/2011, COPENHAGEN                                     --
--  BY MATS H.M. OLSSON AND CHRESTEN R. SONDERGARD                           --
--                                                                           --
--  VERSION 3.1,  07/01/2011, COPENHAGEN                                     --
--  BY CHRESTEN R. SONDERGARD AND MATS H.M. OLSSON                           --
--                                                                           --
--  VERSION 3.2,  06/17/2020, PLANET EARTH                                   --
--  SEE ABOVE FOR AUTHORS                                                    --
--                                                                           --
-------------------------------------------------------------------------------
"""
    return str_


def get_references_header():
    """Create the 'references' part of output file.

    Returns:
        string
    """
    str_ = """
-------------------------------------------------------------------------------
References:

Very Fast Empirical Prediction and Rationalization of Protein pKa Values.
Hui Li, Andrew D. Robertson and Jan H. Jensen. PROTEINS: Structure, Function,
and Bioinformatics. 61:704-721 (2005)

Very Fast Prediction and Rationalization of pKa Values for Protein-Ligand
Complexes.  Delphine C. Bas, David M. Rogers and Jan H. Jensen.  PROTEINS:
Structure, Function, and Bioinformatics 73:765-783 (2008)

PROPKA3: Consistent Treatment of Internal and Surface Residues in Empirical
pKa predictions.  Mats H.M. Olsson, Chresten R. Sondergard, Michal Rostkowski,
and Jan H. Jensen.  Journal of Chemical Theory and Computation, 7(2):525-537
(2011)

Improved Treatment of Ligands and Coupling Effects in Empirical Calculation
and Rationalization of pKa Values.  Chresten R. Sondergaard, Mats H.M. Olsson,
Michal Rostkowski, and Jan H. Jensen.  Journal of Chemical Theory and
Computation, (2011)
-------------------------------------------------------------------------------
"""
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
    str_ = """
---------  -----   ------   ---------------------    --------------    --------------    --------------
                            DESOLVATION  EFFECTS       SIDECHAIN          BACKBONE        COULOMBIC
 RESIDUE    pKa    BURIED     REGULAR      RE        HYDROGEN BOND     HYDROGEN BOND      INTERACTION
---------  -----   ------   ---------   ---------    --------------    --------------    --------------
"""
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
    """Draw the line-Johnny Cash would have been proud-or actually Aerosmith!

    NOTE - Johnny Cash walked the line.

    Returns:
        string
    """
    str_ = ""
    str_ += "-" * 104
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
        return ""
    # for long list, use condensed formatting
    if len(list_) > 10:
        res = "Condensed form:\n"
        for i, group1 in enumerate(list_):
            for group2 in list_[i:]:
                if interaction(group1, group2):
                    res += f"Coupling: {group1.label:>9s} - "
                    res += f"{group2.label:>9s}\n"
        return res
    # Name and map header
    res = f"{name:s}\n{'':>12s}"
    for group in list_:
        res += f"{group.label:>9s} | "
    # do the map
    for group1 in list_:
        res += f"\n{group1.label:<12s}"
        for group2 in list_:
            tag = ""
            if interaction(group1, group2):
                tag = "    X     "
            res += f"{tag:>10s}| "
    return res


def write_pdb_for_atoms(atoms, filename, make_conect_section=False):
    """Write out PDB file for atoms.

    Args:
        atoms:  list of atoms
        filename:  name of file
        make_conect_section:  generate a CONECT PDB section
    """
    out = open_file_for_writing(filename)
    for atom in atoms:
        out.write(atom.make_pdb_line())
    if make_conect_section:
        for atom in atoms:
            out.write(atom.make_conect_line())
    out.close()


def get_bond_order(atom1, atom2):
    """Get the order of a bond between two atoms.

    Args:
        atom1:  first atom in bond
        atom2:  second atom in bond
    Returns:
        string with bond type
    """
    type_ = "1"
    pi_electrons1 = atom1.num_pi_elec_2_3_bonds
    pi_electrons2 = atom2.num_pi_elec_2_3_bonds
    if ".ar" in atom1.sybyl_type:
        pi_electrons1 -= 1
    if ".ar" in atom2.sybyl_type:
        pi_electrons2 -= 1
    if pi_electrons1 > 0 and pi_electrons2 > 0:
        type_ = f"{min(pi_electrons1, pi_electrons2) + 1:d}"
    if ".ar" in atom1.sybyl_type and ".ar" in atom2.sybyl_type:
        type_ = "ar"
    return type_


def write_mol2_for_atoms(atoms, filename):
    """Write out MOL2 file for atoms.

    Args:
        atoms:  list of atoms
        filename:  name of file
    """
    # TODO - header needs to be converted to format string
    header = "@<TRIPOS>MOLECULE\n\n{natom:d} {id:d}\nSMALL\nUSER_CHARGES\n"
    atoms_section = "@<TRIPOS>ATOM\n"
    for i, atom in enumerate(atoms):
        atoms_section += atom.make_mol2_line(i + 1)
    bonds_section = "@<TRIPOS>BOND\n"
    id_ = 1
    for i, atom1 in enumerate(atoms):
        for j, atom2 in enumerate(atoms, i + 1):
            if atom1 in atom2.bonded_atoms:
                type_ = get_bond_order(atom1, atom2)
                bonds_section += f"{id_:>7d} {i + 1:>7d} "
                bonds_section += f"{j + 1:>7d} {type_:>7s}\n"
                id_ += 1
    substructure_section = "@<TRIPOS>SUBSTRUCTURE\n\n"
    if len(atoms) > 0:
        substructure_section = (
            f"@<TRIPOS>SUBSTRUCTURE\n{atoms[0].res_num:<7d} "
            f"{atoms[0].res_name:>10s} "
            f"{atoms[0].numb:>7d}\n"
        )
    out = open_file_for_writing(filename)
    out.write(header.format(natom=len(atoms), id=id_ - 1))
    out.write(atoms_section)
    out.write(bonds_section)
    out.write(substructure_section)
    out.close()


def write_propka(molecular_container, filename):
    """Write PROPKA input file for molecular container.

    Args:
        molecular_container:  molecular container
        filename:  output file name
    """
    out = open_file_for_writing(filename)
    for conformation_name in molecular_container.conformation_names:
        out.write(f"MODEL {conformation_name}\n")
        # write atoms
        for atom in molecular_container.conformations[conformation_name].atoms:
            out.write(atom.make_input_line())
        # write bonds
        for atom in molecular_container.conformations[conformation_name].atoms:
            out.write(atom.make_conect_line())
        # write covalently coupled groups
        for group in molecular_container.conformations[
            conformation_name
        ].groups:
            out.write(group.make_covalently_coupled_line())
        # write non-covalently coupled groups
        for group in molecular_container.conformations[
            conformation_name
        ].groups:
            out.write(group.make_non_covalently_coupled_line())
        out.write("ENDMDL\n")
    out.close()
