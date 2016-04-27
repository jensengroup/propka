#!/usr/bin/python

from __future__ import division
from __future__ import print_function

from propka.vector_algebra import *
import propka.bonds, propka.pdb, propka.atom
from propka.lib import info, warning, debug

class Protonate:
    """ Protonates atoms using VSEPR theory """

    def __init__(self, verbose=False):
        self.verbose=verbose

        self.valence_electrons = {'H': 1,
                                  'He':2,
                                  'Li':1,
                                  'Be':2,
                                  'B': 3,
                                  'C': 4,
                                  'N': 5,
                                  'O': 6,
                                  'F': 7,
                                  'Ne':8,
                                  'Na':1,
                                  'Mg':2,
                                  'Al':3,
                                  'Si':4,
                                  'P': 5,
                                  'S': 6,
                                  'Cl':7,
                                  'Ar':8,
                                  'K': 1,
                                  'Ca':2,
                                  'Sc':2,
                                  'Ti':2,
                                  'Va':2,
                                  'Cr':1,
                                  'Mn':2,
                                  'Fe':2,
                                  'Co':2,
                                  'Ni':2,
                                  'Cu':1,
                                  'Zn':2,
                                  'Ga':3,
                                  'Ge':4,
                                  'As':5,
                                  'Se':6,
                                  'Br':7,
                                  'Kr':8,
                                  'I':7,
                                }





        self.standard_charges= {'ARG-NH1':1.0,
                                'ASP-OD2':-1.0,
                                'GLU-OE2':-1.0,
                                'HIS-ND1':1.0,
                                'LYS-NZ':1.0,
                                'N+':1.0,
                                'C-':-1.0}


        self.sybyl_charges = {'N.pl3':+1,
                              'N.3':+1,
                              'N.4':+1,
                              'N.ar':+1,
                              'O.co2-':-1}


        self.bond_lengths = {'C':1.09,
                             'N':1.01,
                             'O':0.96,
                             'F':0.92,
                             'Cl':1.27,
                             'Br':1.41,
                             'I':1.61,
                             'S':1.35}


        # protonation_methods[steric_number] = method
        self.protonation_methods = {4:self.tetrahedral,
                                    3:self.trigonal}


        return




    def protonate(self, molecules):
        """ Will protonate all atoms in the molecular container """

        debug('----- Protonation started -----')
        # Remove all currently present hydrogen atoms
        self.remove_all_hydrogen_atoms(molecules)

        # protonate all atoms
        for name in molecules.conformation_names:
            non_H_atoms = molecules.conformations[name].get_non_hydrogen_atoms()

            for atom in non_H_atoms:
                self.protonate_atom(atom)

            # fix hydrogen names
            #self.set_proton_names(non_H_atoms)

        return


    def remove_all_hydrogen_atoms(self, molecular_container):
        for name in molecular_container.conformation_names:
            molecular_container.conformations[name].atoms = molecular_container.conformations[name].get_non_hydrogen_atoms()
        return


    def set_charge(self, atom):
        # atom is a protein atom
        if atom.type=='atom':
            key = '%3s-%s'%(atom.resName, atom.name)
            if atom.terminal:
                debug(atom.terminal)
                key=atom.terminal
            if key in list(self.standard_charges.keys()):
                atom.charge = self.standard_charges[key]
                debug('Charge', atom, atom.charge)
                atom.charge_set = True
        # atom is a ligand atom
        elif atom.type=='hetatm':
            if atom.sybyl_type in list(self.sybyl_charges.keys()):
                atom.charge = self.sybyl_charges[atom.sybyl_type]
                atom.sybyl_type = atom.sybyl_type.replace('-','')
                atom.charge_set = True

        return

    def protonate_atom(self, atom):
        if atom.is_protonated: return
        if atom.element == 'H': return

        self.set_charge(atom)
        self.set_number_of_protons_to_add(atom)
        self.set_steric_number_and_lone_pairs(atom)
        self.add_protons(atom)
        atom.is_protonated = True
        return

    def set_proton_names(self, heavy_atoms):

        for heavy_atom in heavy_atoms:
            i = 1
            for bonded in heavy_atom.bonded_atoms:

                if bonded.element == 'H':
                    bonded.name+='%d'%i
                    i+=1


        return


    def set_number_of_protons_to_add(self, atom):
        debug('*'*10)
        debug('Setting number of protons to add for',atom)
        atom.number_of_protons_to_add  = 8
        debug('                  %4d'%8)
        atom.number_of_protons_to_add -= self.valence_electrons[atom.element]
        debug('Valence eletrons: %4d'%-self.valence_electrons[atom.element])
        atom.number_of_protons_to_add -= len(atom.bonded_atoms)
        debug('Number of bonds:  %4d'%- len(atom.bonded_atoms))
        atom.number_of_protons_to_add -= atom.number_of_pi_electrons_in_double_and_triple_bonds
        debug('Pi electrons:     %4d'%-atom.number_of_pi_electrons_in_double_and_triple_bonds)
        atom.number_of_protons_to_add += int(atom.charge)
        debug('Charge:           %4.1f'%atom.charge)

        debug('-'*10)
        debug(atom.number_of_protons_to_add)

        return

    def set_steric_number_and_lone_pairs(self, atom):

        # If we already did this, there is no reason to do it again
        if atom.steric_number_and_lone_pairs_set:
            return

        debug('='*10)
        debug('Setting steric number and lone pairs for',atom)

        # costumly set the N backbone atoms up for peptide bond trigonal planer shape
        #if atom.name == 'N' and len(atom.bonded_atoms) == 2:
        #    atom.steric_number = 3
        #    atom.number_of_lone_pairs = 0
        #    self.display 'Peptide bond: steric number is %d and number of lone pairs is %s'%(atom.steric_number,
         #                                                                             atom.number_of_lone_pairs)
        #    return


        atom.steric_number = 0

        debug('%65s: %4d'%('Valence electrons',self.valence_electrons[atom.element]))
        atom.steric_number += self.valence_electrons[atom.element]

        debug('%65s: %4d'%('Number of bonds',len(atom.bonded_atoms)))
        atom.steric_number += len(atom.bonded_atoms)

        debug('%65s: %4d'%('Number of hydrogen atoms to add',atom.number_of_protons_to_add))
        atom.steric_number += atom.number_of_protons_to_add

        debug('%65s: %4d'%('Number of pi-electrons in double and triple bonds(-)',atom.number_of_pi_electrons_in_double_and_triple_bonds))
        atom.steric_number -= atom.number_of_pi_electrons_in_double_and_triple_bonds

        debug('%65s: %4d'%('Number of pi-electrons in conjugated double and triple bonds(-)',atom.number_of_pi_electrons_in_conjugate_double_and_triple_bonds))
        atom.steric_number -= atom.number_of_pi_electrons_in_conjugate_double_and_triple_bonds

        debug('%65s: %4d'%('Number of donated co-ordinated bonds',0))
        atom.steric_number += 0

        debug('%65s: %4.1f'%('Charge(-)',atom.charge))
        atom.steric_number -= atom.charge

        atom.steric_number = math.floor(atom.steric_number/2.0)

        atom.number_of_lone_pairs = atom.steric_number - len(atom.bonded_atoms) - atom.number_of_protons_to_add

        debug('-'*70)
        debug('%65s: %4d'%('Steric number',atom.steric_number))
        debug('%65s: %4d'%('Number of lone pairs',atom.number_of_lone_pairs))

        atom.steric_number_and_lone_pairs_set = True

        return


    def add_protons(self, atom):
        # decide which method to use
        debug('PROTONATING',atom)
        if atom.steric_number in list(self.protonation_methods.keys()):
            self.protonation_methods[atom.steric_number](atom)
        else:
            warning('Do not have a method for protonating', atom, '(steric number: %d)' % atom.steric_number)

        return


    def trigonal(self, atom):
        debug('TRIGONAL - %d bonded atoms'%(len(atom.bonded_atoms)))
        rot_angle = math.radians(120.0)

        c = vector(atom1 = atom)

        # 0 bonds
        if len(atom.bonded_atoms) == 0:
            pass

        # 1 bond
        if len(atom.bonded_atoms) == 1 and atom.number_of_protons_to_add > 0:
            # Add another atom with the right angle to the first one
            a = vector(atom1 = atom, atom2 = atom.bonded_atoms[0])
            # use plane of bonded trigonal atom - e.g. arg

            self.set_steric_number_and_lone_pairs(atom.bonded_atoms[0])
            if atom.bonded_atoms[0].steric_number == 3 and len(atom.bonded_atoms[0].bonded_atoms)>1:
                # use other atoms bonded to the neighbour to establish the plane, if possible
                other_atom_indices = []
                for i in range(len(atom.bonded_atoms[0].bonded_atoms)):
                    if atom.bonded_atoms[0].bonded_atoms[i] != atom:
                        other_atom_indices.append(i)


                v1 = vector(atom1 = atom, atom2 = atom.bonded_atoms[0])
                v2 = vector(atom1 = atom.bonded_atoms[0],
                            atom2 = atom.bonded_atoms[0].bonded_atoms[other_atom_indices[0]])

                axis = v1**v2

                 # this is a trick to make sure that the order of atoms doesn't influence
                 # the final postions of added protons
                if len(other_atom_indices)>1:
                    v3 = vector(atom1 = atom.bonded_atoms[0],
                                atom2 = atom.bonded_atoms[0].bonded_atoms[other_atom_indices[1]])

                    axis2 = v1**v3

                    if axis * axis2>0:
                        axis = axis+axis2
                    else:
                        axis = axis-axis2

            else:
                axis = a.orthogonal()

            a = rotate_vector_around_an_axis(rot_angle, axis, a)
            a = self.set_bond_distance(a, atom.element)
            self.add_proton(atom, c+a)

        # 2 bonds
        if len(atom.bonded_atoms) == 2 and atom.number_of_protons_to_add > 0:
            # Add another atom with the right angle to the first two
            a1 = vector(atom1 = atom, atom2 = atom.bonded_atoms[0]).rescale(1.0)
            a2 = vector(atom1 = atom, atom2 = atom.bonded_atoms[1]).rescale(1.0)

            new_a = -a1 - a2
            new_a = self.set_bond_distance(new_a, atom.element)
            self.add_proton(atom, c+new_a)


        return


    def tetrahedral(self, atom):
        debug('TETRAHEDRAL - %d bonded atoms'%(len(atom.bonded_atoms)))
        rot_angle = math.radians(109.5)

        # sanity check
        # if atom.number_of_protons_to_add + len(atom.bonded_atoms) != 4:
        # self.display 'Error: Attempting tetrahedral structure with %d bonds'%(atom.number_of_protons_to_add +
        #                                                                len(atom.bonded_atoms))

        c = vector(atom1 = atom)

        # 0 bonds
        if len(atom.bonded_atoms) == 0:
            pass

        # 1 bond
        if len(atom.bonded_atoms) == 1 and atom.number_of_protons_to_add > 0:
            # Add another atom with the right angle to the first one
            a = vector(atom1 = atom, atom2 = atom.bonded_atoms[0])
            axis = a.orthogonal()
            a = rotate_vector_around_an_axis(rot_angle, axis, a)
            a = self.set_bond_distance(a, atom.element)
            self.add_proton(atom, c+a)

        # 2 bonds
        if len(atom.bonded_atoms) == 2 and atom.number_of_protons_to_add > 0:
            # Add another atom with the right angle to the first two
            a1 = vector(atom1 = atom, atom2 = atom.bonded_atoms[0]).rescale(1.0)
            a2 = vector(atom1 = atom, atom2 = atom.bonded_atoms[1]).rescale(1.0)

            axis = a1 + a2

            new_a =  rotate_vector_around_an_axis(math.radians(90), axis, -a1)
            new_a = self.set_bond_distance(new_a, atom.element)
            self.add_proton(atom, c+new_a)

        # 3 bonds
        if len(atom.bonded_atoms) == 3 and atom.number_of_protons_to_add > 0:
            a1 = vector(atom1 = atom, atom2 = atom.bonded_atoms[0]).rescale(1.0)
            a2 = vector(atom1 = atom, atom2 = atom.bonded_atoms[1]).rescale(1.0)
            a3 = vector(atom1 = atom, atom2 = atom.bonded_atoms[2]).rescale(1.0)

            new_a =  -a1-a2-a3
            new_a = self.set_bond_distance(new_a, atom.element)
            self.add_proton(atom, c+new_a)

        return


    def add_proton(self, atom, position):
        # Create the new proton
        new_H = propka.atom.Atom()
        new_H.setProperty(numb    = None,
                          name    = 'H%s'%atom.name[1:],
                          resName = atom.resName,
                          chainID = atom.chainID,
                          resNumb = atom.resNumb,
                          x       = round(position.x,3), # round of to three digimal points
                          y       = round(position.y,3), # to avoid round-off differences
                          z       = round(position.z,3), # when input file
                          occ     = None,
                          beta    = None)
        new_H.element = 'H'
        new_H.type = atom.type

        new_H.bonded_atoms = [atom]
        new_H.charge = 0
        new_H.steric_number = 0
        new_H.number_of_lone_pairs = 0
        new_H.number_of_protons_to_add = 0
        new_H.number_of_pi_electrons_in_double_and_triple_bonds = 0
        new_H.is_protonates = True

        atom.bonded_atoms.append(new_H)
        atom.number_of_protons_to_add -=1
        atom.conformation_container.add_atom(new_H)

        # update names of all protons on this atom
        new_H.residue_label = "%-3s%4d%2s" % (new_H.name,new_H.resNumb, new_H.chainID)
        no_protons = atom.count_bonded_elements('H')
        if no_protons > 1:
            i = 1
            for proton in atom.get_bonded_elements('H'):
                proton.name = 'H%s%d'%(atom.name[1:],i)
                proton.residue_label = "%-3s%4d%2s" % (proton.name,proton.resNumb, proton.chainID)
                i+=1


        debug('added',new_H, 'to',atom)
        return

    def set_bond_distance(self, a, element):
        d = 1.0
        if element in list(self.bond_lengths.keys()):
            d = self.bond_lengths[element]
        else:
            warning('Bond length for %s not found, using the standard value of %f' % (element, d))

        a = a.rescale(d)

        return a


if __name__ == '__main__':
    import protein, pdb, sys,os
    arguments = sys.argv
    if len(arguments) != 2:
        info('Usage: protonate.py <pdb_file>')
        sys.exit(0)

    filename = arguments[1]
    if not os.path.isfile(filename):
        info('Error: Could not find \"%s\"' % filename)
        sys.exit(1)


    p = Protonate()
    pdblist = pdb.readPDB(filename)
    my_protein = protein.Protein(pdblist,'test.pdb')

    p.remove_all_hydrogen_atoms_from_protein(my_protein)
    my_protein.writePDB('before_protonation.pdb')

    p.protonate_protein(my_protein)

    ## write out protonated file
    my_protein.writePDB('protonated.pdb')
