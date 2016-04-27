#!/usr/bin/python

from __future__ import division
from __future__ import print_function

import sys

import propka.calculations
from propka.vector_algebra import *
from propka.lib import info, warning


all_sybyl_types = [
    'C.3',   #  carbon sp3
    'H',     #  hydrogen
    'C.2',   #  carbon sp2
    'H.spc', #  hydrogen in Single Point Charge (SPC) water model
    'C.1',   #  carbon sp
    'H.t3p', #  hydrogen in Transferable intermolecular Potential (TIP3P) water model
    'C.ar',  #  carbon aromatic
    'LP',    #  lone pair
    'C.cat', #  carbocation (C+) used only in a guadinium group
    'Du',    #  dummy atom
    'N.3',   #  nitrogen sp3
    'Du.C',  #  dummy carbon
    'N.2',   #  nitrogen sp2
    'Any',   #  any atom
    'N.1',   #  nitrogen sp
    'Hal',   #  halogen
    'N.ar',  #  nitrogen aromatic
    'Het',   #  heteroatom = N, O, S, P
    'N.am',  #  nitrogen amide
    'Hev',   #  heavy atom (non hydrogen)
    'N.pl3', #  nitrogen trigonal planar
    'Li',    #  lithium
    'N.4',   #  nitrogen sp3 positively charged
    'Na',    #  sodium
    'O.3',   #  oxygen sp3
    'Mg',    #  magnesium
    'O.2',   #  oxygen sp2
    'Al',    #  aluminum
    'O.co2', #  oxygen in carboxylate and phosphate groups
    'Si',    #  silicon
    'O.spc', #  oxygen in Single Point Charge (SPC) water model
    'K',     #  potassium
    'O.t3p', #  oxygen in Transferable Intermolecular Potential (TIP3P) water model
    'Ca',    #  calcium
    'S.3',   #  sulfur sp3
    'Cr.th', #  chromium (tetrahedral)
    'S.2',   #  sulfur sp2
    'Cr.oh', #  chromium (octahedral)
    'S.O',   #  sulfoxide sulfur
    'Mn',    #  manganese
    'S.O2',  #  sulfone sulfur
    'Fe',    #  iron
    'P.3',   #  phosphorous sp3
    'Co.oh', #  cobalt (octahedral)
    'F',     #  fluorine
    'Cu',    #  copper
    'Cl',    #  chlorine
    'Zn',    #  zinc
    'Br',    #  bromine
    'Se',    #  selenium
    'I',     #  iodine
    'Mo',    #  molybdenum
    'Sn']    #  tin


#propka_input_types = ['1P','1N','2P','2N']
#for type in all_sybyl_types:
#    temp = type.replace('.','')
#    if len(temp)>3:
#        temp = temp[0:3]
#    propka_input_types.append(temp)
#
#for t in propka_input_types:
#    print (t)


propka_input_types = [
    '1P',
    '1N',
    '2P',
    '2N',
    'C3',
    'H',
    'C2',
    'Hsp',
    'C1',
    'Ht3',
    'Car',
    'LP',
    'Cca',
    'Du',
    'N3',
    'DuC',
    'N2',
    'Any',
    'N1',
    'Hal',
    'Nar',
    'Het',
    'Nam',
    'Hev',
    'Npl',
    'Li',
    'N4',
    'Na',
    'O3',
    'Mg',
    'O2',
    'Al',
    'Oco',
    'Si',
    'Osp',
    'K',
    'Ot3',
    'Ca',
    'S3',
    'Crt',
    'S2',
    'Cro',
    'SO',
    'Mn',
    'SO2',
    'Fe',
    'P3',
    'Coo',
    'F',
    'Cu',
    'Cl',
    'Zn',
    'Br',
    'Se',
    'I',
    'Mo',
    'Sn']


max_C_double_bond = 1.3
max_C_triple_bond = 1.2

max_C_double_bond_squared = max_C_double_bond*max_C_double_bond
max_C_triple_bond_squared = max_C_triple_bond*max_C_triple_bond






def assign_sybyl_type(atom):
    # check if we already have assigned a name to this atom
    if atom.sybyl_assigned:
        #info(atom.name,'already assigned')
        return

    # find some properties of the atom
    ring_atoms = is_ring_member(atom)
    aromatic   = is_aromatic_ring(ring_atoms)
    planar     = is_planar(atom)
    bonded_elements = {}
    for i in range(len(atom.bonded_atoms)):
        bonded_elements[i]=atom.bonded_atoms[i].element



    # Aromatic carbon/nitrogen
    if aromatic:
        for ra in ring_atoms:
            if ra.element in ['C','N']:
                set_type(ra, ra.element+'.ar')
        return

    # check for amide
    if atom.element in ['O','N','C']:
        O=None
        N=None
        C=None
        # oxygen, nitrogen
        if atom.element in ['O','N']:
            for b in atom.get_bonded_elements('C'):
                for bb in b.bonded_atoms:
                    if (bb.element =='N' and atom.element == 'O'):
                        O=atom
                        C=b
                        N=bb
                    elif (bb.element =='O' and atom.element == 'N'):
                        N=atom
                        C=b
                        O=bb
        # carbon
        if atom.element == 'C':
            nitrogens = atom.get_bonded_elements('N')
            oxygens = atom.get_bonded_elements('O')
            if len(nitrogens)==1 and len(oxygens)==1:
                C = atom
                N = nitrogens[0]
                O = oxygens[0]


        if C and N and O:
            # make sure that the Nitrogen is not aromatic and that it has two heavy atom bonds
            if not is_aromatic_ring(is_ring_member(N)) and len(N.get_bonded_heavy_atoms())==2:
                set_type(N,'N.am')
                set_type(C,'C.2')
                set_type(O,'O.2')
                return


    if atom.element=='C':
        # check for carboxyl
        if len(atom.bonded_atoms)==3 and list(bonded_elements.values()).count('O')==2:
            i1 = list(bonded_elements.values()).index('O')
            i2 = list(bonded_elements.values()).index('O',i1+1)
            if len(atom.bonded_atoms[i1].bonded_atoms)==1 and len(atom.bonded_atoms[i2].bonded_atoms)==1:
                set_type(atom.bonded_atoms[i1],'O.co2-')
                set_type(atom.bonded_atoms[i2],'O.co2')
                set_type(atom,'C.2')
                return



        # sp carbon
        if len(atom.bonded_atoms)<=2:
            for b in atom.bonded_atoms:
                if propka.calculations.squared_distance(atom, b) < max_C_triple_bond_squared:
                    set_type(atom,'C.1')
                    set_type(b,b.element+'.1')
            if atom.sybyl_assigned:
                return

        # sp2 carbon
        if planar:
            set_type(atom,'C.2')
            # check for N.pl3
            for b in atom.bonded_atoms:
                if b.element=='N':
                    if len(b.bonded_atoms)<3 or is_planar(b):
                        set_type(b,'N.pl3')
            return

        # sp3 carbon
        set_type(atom, 'C.3')
        return

    # Nitrogen
    if atom.element == 'N':
        # check for planar N
        if len(atom.bonded_atoms)==1:
            if is_planar(atom.bonded_atoms[0]):
                set_type(atom,'N.pl3')
                return

        if planar:
            set_type(atom,'N.pl3')
            return

        set_type(atom,'N.3')
        return

    # Oxygen
    if atom.element == 'O':
        set_type(atom,'O.3')

        if len(atom.bonded_atoms) == 1:
            # check for carboxyl
            if atom.bonded_atoms[0].element == 'C':
                the_carbon = atom.bonded_atoms[0]
                if len(the_carbon.bonded_atoms)==3 and the_carbon.count_bonded_elements('O')==2:
                    [O1,O2] = the_carbon.get_bonded_elements('O')

                    if len(O1.bonded_atoms)==1 and len(O2.bonded_atoms)==1:
                        set_type(O1,'O.co2-')
                        set_type(O2,'O.co2')
                        set_type(the_carbon,'C.2')
                        return

            # check for X=O
            if propka.calculations.squared_distance(atom, atom.bonded_atoms[0]) < max_C_double_bond_squared:
                set_type(atom,'O.2')
                if atom.bonded_atoms[0].element=='C':
                    set_type(atom.bonded_atoms[0],'C.2')
        return


    # Sulphur
    if atom.element == 'S':
        # check for SO2
        if list(bonded_elements.values()).count('O')==2:
            i1 = list(bonded_elements.values()).index('O')
            i2 = list(bonded_elements.values()).index('O',i1+1)
            set_type(atom.bonded_atoms[i1],'O.2')
            set_type(atom.bonded_atoms[i2],'O.2')
            set_type(atom,'S.o2')
            return

        # check for SO4
        if list(bonded_elements.values()).count('O')==4:
            no_O2 = 0
            for i in range(len(atom.bonded_atoms)):
                if len(atom.bonded_atoms[i].bonded_atoms)==1 and no_O2<2:
                    set_type(atom.bonded_atoms[i],'O.2')
                    no_O2+=1
                else:
                    set_type(atom.bonded_atoms[i],'O.3')

        set_type(atom,'S.3')


        return


    # Phosphorus
    if atom.element == 'P':
        set_type(atom,'P.3')

        # check for phosphate group
        bonded_oxygens = atom.get_bonded_elements('O')
        for o in bonded_oxygens: set_type(o,'O.3')
#         if len(bonded_oxygens)>=3:
#             # find oxygens only bonded to current phosphorus
#             bonded_oxygens_1 = [o for o in bonded_oxygens if len(o.get_bonded_heavy_atoms())==1]
#             # find the closest oxygen ...
#             closest_oxygen = min(bonded_oxygens_1,
#                                  key= lambda o:propka.calculations.squared_distance(atom,o))
#             # ... and set it to O.2
#             set_type(closest_oxygen,'O.2')

        return



    element = atom.element.capitalize()
    set_type(atom,element)
    # info('Using element as type for %s'%atom.element)

    return


def is_ring_member(atom):
    return identify_ring(atom,atom,0,[])

def identify_ring(this_atom, original_atom, number, past_atoms):
    number+=1
    past_atoms=past_atoms+[this_atom]
    return_atoms = []
    if number > 10:
        return return_atoms

    for atom in this_atom.get_bonded_heavy_atoms():
        if atom == original_atom and number>2:
            return past_atoms

        if atom not in past_atoms:
            these_return_atoms = identify_ring(atom, original_atom, number, past_atoms)
            if len(these_return_atoms) > 0:
                if len(return_atoms)>len(these_return_atoms) or len(return_atoms)==0:
                    return_atoms = these_return_atoms

    return return_atoms




def set_type(atom,type):
    #info(atom, '->',type)
    atom.sybyl_type = type
    atom.sybyl_assigned=True
    return




def is_planar(atom):
    """ Finds out if atom forms a plane together with its bonded atoms"""
    atoms = [atom]+atom.bonded_atoms
    return are_atoms_planar(atoms)

def are_atoms_planar(atoms):
    if len(atoms)==0:
        return False
    if len(atoms)<4:
        return False
    v1 = vector(atom1=atoms[0], atom2=atoms[1])
    v2 = vector(atom1=atoms[0], atom2=atoms[2])
    n = (v1**v2).rescale(1.0)

    margin = 0.20
    for b in atoms[3:]:
        v = vector(atom1=atoms[0], atom2=b).rescale(1.0)
        #info(atoms[0],abs(v*n) )
        if abs(v*n)>margin:
            return False

    return True

def is_aromatic_ring(atoms):
    if len(atoms)<5:
        return False

    for i in range(len(atoms)):
        if not are_atoms_planar(atoms[i:]+atoms[:i]):
            return False

    return True




