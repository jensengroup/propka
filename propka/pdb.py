
from __future__ import division
from __future__ import print_function

import string, sys, copy

import propka.lib
from propka.lib import info, warning

from propka.atom import Atom
from propka.conformation_container import Conformation_container

expected_atom_numbers = {'ALA':5,
                         'ARG':11,
                         'ASN':8,
                         'ASP':8,
                         'CYS':6,
                         'GLY':4,
                         'GLN':9,
                         'GLU':9,
                         'HIS':10,
                         'ILE':8,
                         'LEU':8,
                         'LYS':9,
                         'MET':8,
                         'PHE':11,
                         'PRO':7,
                         'SER':6,
                         'THR':7,
                         'TRP':14,
                         'TYR':12,
                         'VAL':7}


def read_pdb(pdb_file, parameters, molecule):
    conformations = {}

    # read in all atoms in the file
    lines = get_atom_lines_from_pdb(pdb_file, ignore_residues = parameters.ignore_residues, keep_protons = molecule.options.keep_protons, chains=molecule.options.chains)
    for (name, atom) in lines:
        if not name in conformations.keys():
            conformations[name] = Conformation_container(name=name, parameters=parameters, molecular_container=molecule)
        conformations[name].add_atom(atom)

    # make a sorted list of conformation names
    names = sorted(conformations.keys(), key=propka.lib.conformation_sorter)

    return [conformations, names]

def protein_precheck(conformations, names):

    for name in names:
        atoms = conformations[name].atoms

        # Group the atoms by their residue:
        atoms_by_residue = {}
        for a in atoms:
            if a.element != 'H':
                res_id = resid_from_atom(a)
                try:
                    atoms_by_residue[res_id].append(a)
                except KeyError:
                    atoms_by_residue[res_id] = [a]

        for res_id, res_atoms in atoms_by_residue.items():
            resname = res_atoms[0].resName
            residue_label = '%3s%5s'%(resname, res_id)

            # ignore ligand residues
            if resname not in expected_atom_numbers:
                continue

            # check for c-terminal
            if 'C-' in [a.terminal for a in res_atoms]:
                if len(res_atoms) != expected_atom_numbers[resname]+1:
                    warning('Unexpected number (%d) of atoms in residue %s in conformation %s' % (len(res_atoms), residue_label, name))
                continue

            # check number of atoms in residue
            if len(res_atoms) != expected_atom_numbers[resname]:
                warning('Unexpected number (%d) of atoms in residue %s in conformation %s' % (len(res_atoms), residue_label, name))

    return

def resid_from_atom(a):
    return '%4d %s %s'%(a.resNumb,a.chainID,a.icode)


def get_atom_lines_from_pdb(pdb_file, ignore_residues = [], keep_protons=False, tags = ['ATOM  ', 'HETATM'], chains=None):

    lines = propka.lib.open_file_for_reading(pdb_file).readlines()
    nterm_residue = 'next_residue'
    old_residue = None
    terminal = None
    model = 1


    for line in lines:
        tag = line[0:6]

        # set the model number
        if tag == 'MODEL ':
            model = int(line[6:])
            nterm_residue = 'next_residue'

        if tag == 'TER   ':
            nterm_residue = 'next_residue'

        if tag in tags:
            alt_conf_tag = line[16]
            residue_name = line[12:16]
            residue_number = line[22:26]

            # check if we want this residue
            if line[17:20] in ignore_residues:
                continue
            if chains and line[21] not in chains:
                continue

            # set the Nterm residue number - nessecary because we may need to
            # identify more than one N+ group for structures with alt_conf tags
            if nterm_residue == 'next_residue' and tag == 'ATOM  ':
                # make sure that we reached a new residue - nessecary if OXT is not the last atom in
                # the previous residue
                if old_residue != residue_number:
                    nterm_residue = residue_number
                    old_residue = None


            # Identify the configuration
            # convert digits to letters
            if alt_conf_tag in '123456789':
                alt_conf_tag = chr(ord(alt_conf_tag)+16)
            if alt_conf_tag == ' ':
                alt_conf_tag = 'A'
            conformation = '%d%s'%(model, alt_conf_tag)

            # set the terminal
            if  tag == 'ATOM  ':
                if residue_name.strip() == 'N' and nterm_residue == residue_number:
                    terminal = 'N+'
                if  residue_name.strip() in ['OXT','O\'\'']:
                    terminal = 'C-'
                    nterm_residue = 'next_residue'
                    old_residue = residue_number
            # and yield the atom
            atom = Atom(line=line)
            atom.terminal = terminal
            #if keep_protons:
            #    atom.is_protonated = True
            if not (atom.element == 'H' and not keep_protons): #ignore hydrogen
                yield (conformation, atom)

            terminal = None

    return


def write_pdb(conformation, filename):
    write_pdb_for_atoms(conformation.atoms, filename)
    return

def write_pdb_for_atoms(atoms, filename, make_conect_section=False):
    out = propka.lib.open_file_for_writing(filename)

    for atom in atoms:
        out.write(atom.make_pdb_line())

    if make_conect_section:
        for atom in atoms:
            out.write(atom.make_conect_line())


    out.close()

    return



def write_mol2_for_atoms(atoms, filename):

    header = '@<TRIPOS>MOLECULE\n\n%d %d\nSMALL\nUSER_CHARGES\n'

    atoms_section = '@<TRIPOS>ATOM\n'
    for i in range(len(atoms)):
        atoms_section += atoms[i].make_mol2_line(i+1)


    bonds_section = '@<TRIPOS>BOND\n'
    id = 1
    for i in range(len(atoms)):
        for j in range(i+1,len(atoms)):
            if atoms[i] in atoms[j].bonded_atoms:
                type = get_bond_order(atoms[i],atoms[j])
                bonds_section += '%7d %7d %7d %7s\n'%(id, i+1, j+1, type)
                id+=1

    substructure_section = '@<TRIPOS>SUBSTRUCTURE\n\n'
    if len(atoms)>0:
        substructure_section = '@<TRIPOS>SUBSTRUCTURE\n%-7d %10s %7d\n'%(atoms[0].resNumb,atoms[0].resName,atoms[0].numb)

    out = propka.lib.open_file_for_writing(filename)
    out.write(header%(len(atoms),id-1))
    out.write(atoms_section)
    out.write(bonds_section)
    out.write(substructure_section)
    out.close()

    return

def get_bond_order(atom1, atom2):
    type = '1'
    pi_electrons1 = atom1.number_of_pi_electrons_in_double_and_triple_bonds
    pi_electrons2 = atom2.number_of_pi_electrons_in_double_and_triple_bonds

    if '.ar' in atom1.sybyl_type:
        pi_electrons1 -=1
    if '.ar' in atom2.sybyl_type:
        pi_electrons2 -=1

    if pi_electrons1 > 0 and pi_electrons2 > 0:
        type = '%d'%(min(pi_electrons1, pi_electrons2)+1)

    if '.ar' in atom1.sybyl_type and '.ar' in atom2.sybyl_type:
        type = 'ar'


    return type



def write_input(molecular_container, filename):
    out = propka.lib.open_file_for_writing(filename)

    for conformation_name in molecular_container.conformation_names:
        out.write('MODEL %s\n'%conformation_name)
        # write atoms
        for atom in molecular_container.conformations[conformation_name].atoms:
            out.write(atom.make_input_line())
        # write bonds
        for atom in molecular_container.conformations[conformation_name].atoms:
            out.write(atom.make_conect_line())
        # write covalently coupled groups
        for group in molecular_container.conformations[conformation_name].groups:
            out.write(group.make_covalently_coupled_line())
        # write non-covalently coupled groups
        for group in molecular_container.conformations[conformation_name].groups:
            out.write(group.make_non_covalently_coupled_line())
        out.write('ENDMDL\n')

    out.close()

    return


def read_input(input_file, parameters,molecule):
    conformations = {}

    # read in all atoms in the input file
    lines = get_atom_lines_from_input(input_file)
    for (name, atom) in lines:
        if not name in conformations.keys():
            conformations[name] = Conformation_container(name=name, parameters=parameters, molecular_container=molecule)
        conformations[name].add_atom(atom)

    # make a sorted list of conformation names
    names = sorted(conformations.keys(), key=propka.lib.conformation_sorter)

    return [conformations, names]



def get_atom_lines_from_input(input_file, tags = ['ATOM  ','HETATM']):
    lines = propka.lib.open_file_for_reading(input_file).readlines()
    conformation = ''

    atoms = {}
    numbers = []

    for line in lines:
        tag = line[0:6]

        # set the conformation
        if tag == 'MODEL ':
            conformation = line[6:].strip()

        # found an atom - save it
        if tag in tags:
            atom = Atom(line=line)
            atom.get_input_parameters()
            atom.groups_extracted = 1
            atom.is_protonated = True
            atoms[atom.numb] = atom
            numbers.append(atom.numb)

        # found bonding information - apply it
        if tag == 'CONECT' and len(line)>14:
            conect_numbers = [line[i:i+5] for i in range(6, len(line)-1, 5)]
            center_atom = atoms[int(conect_numbers[0])]
            for n in conect_numbers[1:]:
                b = atoms[int(n)]
                # remember to check for cysteine bridges
                if center_atom.element == 'S' and b.element == 'S':
                        center_atom.cysteine_bridge = True
                        b.cysteine_bridge = True
                # set up bonding
                if not b in center_atom.bonded_atoms:
                    center_atom.bonded_atoms.append(b)
                if not center_atom in b.bonded_atoms:
                    b.bonded_atoms.append(center_atom)

        # found info on covalent coupling
        if tag == 'CCOUPL' and len(line)>14:
            conect_numbers = [line[i:i+5] for i in range(6, len(line)-1, 5)]
            center_atom = atoms[int(conect_numbers[0])]
            for n in conect_numbers[1:]:
                cg = atoms[int(n)]
                center_atom.group.couple_covalently(cg.group)

        # found info on non-covalent coupling
        if tag == 'NCOUPL' and len(line)>14:
            conect_numbers = [line[i:i+5] for i in range(6, len(line)-1, 5)]
            center_atom = atoms[int(conect_numbers[0])]
            for n in conect_numbers[1:]:
                cg = atoms[int(n)]
                center_atom.group.couple_non_covalently(cg.group)

        # this conformation is done - yield the atoms
        if tag == 'ENDMDL':
            for n in numbers:
                yield (conformation, atoms[n])
            # prepare for next conformation
            atoms = {}
            numbers = []


    return
