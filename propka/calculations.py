
from __future__ import division
from __future__ import print_function

import math, propka.protonate, propka.bonds,copy, sys
from propka.lib import info, warning


#
# methods for setting up atom naming, bonding, and protonation
#

def setup_bonding_and_protonation(parameters, molecular_container):
    # make bonds
    my_bond_maker = setup_bonding(parameters, molecular_container)

    # set up ligand atom names
    set_ligand_atom_names(molecular_container)

    # apply information on pi electrons
    my_bond_maker.add_pi_electron_information(molecular_container)

    # Protonate atoms
    if molecular_container.options.protonate_all:
        my_protonator = propka.protonate.Protonate(verbose=False)
        my_protonator.protonate(molecular_container)


    return



def setup_bonding(parameters, molecular_container):
    # make bonds
    my_bond_maker = propka.bonds.bondmaker()
    my_bond_maker.find_bonds_for_molecules_using_boxes(molecular_container)

    return my_bond_maker


def set_ligand_atom_names(molecular_container):
    for name in molecular_container.conformation_names:
        molecular_container.conformations[name].set_ligand_atom_names()
    return



#
# The following methods imitates propka3.0 protonation!
#
def setup_bonding_and_protonation_30_style(parameters, molecular_container):
    # Protonate atoms
    protonate_30_style(molecular_container)

    # make bonds
    my_bond_maker = propka.bonds.bondmaker()
    my_bond_maker.find_bonds_for_molecules_using_boxes(molecular_container)

    return


def protonate_30_style(molecular_container):
    for name in molecular_container.conformation_names:
        info('Now protonating', name)
        # split atom into residues
        curres = -1000000
        residue = []
        O=None
        C=None
        for atom in molecular_container.conformations[name].atoms:
            if atom.resNumb != curres:
                curres = atom.resNumb
                if len(residue)>0:
                    #backbone
                    [O, C]= addBackBoneHydrogen(residue,O,C)
                    #arginine
                    if residue[0].resName == 'ARG':
                        addArgHydrogen(residue)
                    #histidine
                    if residue[0].resName == 'HIS':
                        addHisHydrogen(residue)
                    #tryptophan
                    if residue[0].resName == 'TRP':
                        addTrpHydrogen(residue)
                    #amides
                    if residue[0].resName in ['GLN','ASN']:
                        addAmdHydrogen(residue)


                    residue = []
            if atom.type=='atom':
                residue.append(atom)

    return

def addArgHydrogen(residue):
    """
    Adds Arg hydrogen atoms to residues according to the 'old way'.
    """
    #info('Adding arg H',residue)
    for atom in residue:
        if   atom.name == "CD":
            CD  = atom
        elif atom.name == "CZ":
            CZ  = atom
        elif atom.name == "NE":
            NE  = atom
        elif atom.name == "NH1":
            NH1 = atom
        elif atom.name == "NH2":
            NH2 = atom

    H1 = protonateSP2([CD, NE, CZ])
    H1.name = "HE"
    H2 = protonateDirection([NH1, NE, CZ])
    H2.name = "HN1"
    H3 = protonateDirection([NH1, NE, CD])
    H3.name = "HN2"
    H4 = protonateDirection([NH2, NE, CZ])
    H4.name = "HN3"
    H5 = protonateDirection([NH2, NE, H1])
    H5.name = "HN4"

    return [H1,H2,H3,H4,H5]

def addHisHydrogen(residue):
    """
    Adds His hydrogen atoms to residues according to the 'old way'.
    """
    for atom in residue:
        if   atom.name == "CG":
            CG  = atom
        elif atom.name == "ND1":
            ND  = atom
        elif atom.name == "CD2":
            CD  = atom
        elif atom.name == "CE1":
            CE  = atom
        elif atom.name == "NE2":
            NE  = atom
    HD = protonateSP2([CG, ND, CE])
    HD.name = "HND"
    HE = protonateSP2([CD, NE, CE])
    HE.name = "HNE"
    return

def addTrpHydrogen(residue):
    """
    Adds Trp hydrogen atoms to residues according to the 'old way'.
    """
    CD = None
    NE = None
    DE = None
    for atom in residue:
        if   atom.name == "CD1":
            CD  = atom
        elif atom.name == "NE1":
            NE  = atom
        elif atom.name == "CE2":
            CE  = atom
    if CD == None or NE == None or CE == None:
        str = "Did not find all atoms in %s%4d - in %s" % (self.resName, self.resNumb, "addTrpHydrogen()")
        info(str)
        sys.exit(0)

    HE = protonateSP2([CD, NE, CE])
    HE.name = "HNE"

    return

def addAmdHydrogen(residue):
    """
    Adds Gln & Asn hydrogen atoms to residues according to the 'old way'.
    """
    C = None
    O = None
    N = None
    for atom in residue:
        if   (atom.resName == "GLN" and atom.name == "CD")  or (atom.resName == "ASN" and atom.name == "CG"):
            C = atom
        elif (atom.resName == "GLN" and atom.name == "OE1") or (atom.resName == "ASN" and atom.name == "OD1"):
            O = atom
        elif (atom.resName == "GLN" and atom.name == "NE2") or (atom.resName == "ASN" and atom.name == "ND2"):
            N = atom

    if C == None or O == None or N == None:
        str = "Did not find N, C and/or O in %s%4d - in %s" % (atom.resName, atom.resNumb, "addAmdHydrogen()")
        info(str)
        sys.exit(0)

    H1 = protonateDirection([N, O, C])
    H1.name = "HN1"
    H2 = protonateAverageDirection([N, C, O])
    H2.name = "HN2"

    return

def addBackBoneHydrogen(residue, O, C):
    """
    Adds hydrogen backbone atoms to residues according to the old way; dR is wrong for the N-terminus
    (i.e. first residue) but it doesn't affect anything at the moment. Could be improved, but works
    for now.
    """

    new_C     = None
    new_O     = None
    N     = None


    for atom in residue:
        if atom.name == "N":
            N = atom
        if atom.name == "C":
            new_C = atom
        if atom.name == "O":
            new_O = atom




    if None in [C,O,N]:
        return [new_O,new_C]


    if N.resName == "PRO":
        """ PRO doesn't have an H-atom; do nothing """
    else:
        H = protonateDirection([N, O, C])
        H.name = "H"

    return [new_O,new_C]





def protonateDirection(list):
    """
    Protonates an atom, X1, given a direction (X2 -> X3) [X1, X2, X3]
    """
    X1 = list[0]
    X2 = list[1]
    X3 = list[2]

    dX = (X3.x - X2.x)
    dY = (X3.y - X2.y)
    dZ = (X3.z - X2.z)
    length = math.sqrt( dX*dX + dY*dY + dZ*dZ )
    x = X1.x + dX/length
    y = X1.y + dY/length
    z = X1.z + dZ/length


    H  = make_new_H(X1,x,y,z)
    H.name = "H"


    return H


def protonateAverageDirection(list):
    """
    Protonates an atom, X1, given a direction (X1/X2 -> X3) [X1, X2, X3]
    Note, this one uses the average of X1 & X2 (N & O) unlike the previous
    N - C = O
    """
    X1 = list[0]
    X2 = list[1]
    X3 = list[2]

    dX = (X3.x + X1.x)*0.5 - X2.x
    dY = (X3.y + X1.y)*0.5 - X2.y
    dZ = (X3.z + X1.z)*0.5 - X2.z

    length = math.sqrt( dX*dX + dY*dY + dZ*dZ )
    x = X1.x + dX/length
    y = X1.y + dY/length
    z = X1.z + dZ/length

    H  = make_new_H(X1,x,y,z)
    H.name = "H"



    return H


def protonateSP2(list):
    """
    Protonates a SP2 atom, H2, given a list of [X1, X2, X3]
    X1-X2-X3
    """
    X1 = list[0]
    X2 = list[1]
    X3 = list[2]

    dX = (X1.x + X3.x)*0.5 - X2.x
    dY = (X1.y + X3.y)*0.5 - X2.y
    dZ = (X1.z + X3.z)*0.5 - X2.z
    length = math.sqrt( dX*dX + dY*dY + dZ*dZ )
    x = X2.x - dX/length
    y = X2.y - dY/length
    z = X2.z - dZ/length

    H  = make_new_H(X2,x,y,z)
    H.name = "H"

    return H


def make_new_H(atom, x,y,z):

    new_H = propka.atom.Atom()
    new_H.setProperty(numb    = None,
                      name    = 'H%s'%atom.name[1:],
                      resName = atom.resName,
                      chainID = atom.chainID,
                      resNumb = atom.resNumb,
                      x       = x,
                      y       = y,
                      z       = z,
                      occ     = None,
                      beta    = None)
    new_H.element = 'H'


    new_H.bonded_atoms = [atom]
    new_H.charge = 0
    new_H.steric_number = 0
    new_H.number_of_lone_pairs = 0
    new_H.number_of_protons_to_add = 0
    new_H.number_of_pi_electrons_in_double_and_triple_bonds = 0

    atom.bonded_atoms.append(new_H)
    atom.conformation_container.add_atom(new_H)

    return new_H

######## 3.0 style protonation methods end







#
# Desolvation methods
#


def radial_volume_desolvation(parameters, group):
    all_atoms = group.atom.conformation_container.get_non_hydrogen_atoms()
    volume = 0.0
    group.Nmass = 0
    min_distance_4th = 57.1914 # pow(2.75, 4)

    for atom in all_atoms:
        # ignore atoms in the same residue
        if atom.resNumb == group.atom.resNumb and atom.chainID == group.atom.chainID:
            continue

        sq_dist = squared_distance(group, atom)

        # desolvation
        if sq_dist < parameters.desolv_cutoff_squared:
            # use a default relative volume of 1.0 if the volume of the element is not found in parameters
            # insert check for methyl groups
            if atom.element == 'C' and atom.name not in ['CA','C']:
                dv = parameters.VanDerWaalsVolume['C4']
            else:
                dv = parameters.VanDerWaalsVolume.get(atom.element, 1.0)

            dv_inc = dv/max(min_distance_4th, sq_dist*sq_dist)
#            dv_inc = dv/(sq_dist*sq_dist) - dv/(parameters.desolv_cutoff_squared*parameters.desolv_cutoff_squared)
            volume += dv_inc
        # buried
        if sq_dist < parameters.buried_cutoff_squared:
            group.Nmass += 1

    group.buried = calculate_weight(parameters, group.Nmass)
    scale_factor = calculate_scale_factor(parameters, group.buried)
    volume_after_allowance = max(0.00, volume-parameters.desolvationAllowance)

    group.Emass = group.charge * parameters.desolvationPrefactor * volume_after_allowance * scale_factor
    # Emass, Nmass
    # Elocl, Nlocl -> reorganisation energy (count backbone hydorgen bond acceptors, C=O)



    #info('%s %5.2f %5.2f %4d'%(group, group.buried, group.Emass, group.Nmass))
    return



def contactDesolvation(parameters, group):
    """
    calculates the desolvation according to the Contact Model, the old default
    """

    local_radius = {'ASP': 4.5,
                    'GLU': 4.5,
                    'HIS': 4.5,
                    'CYS': 3.5,
                    'TYR': 3.5,
                    'LYS': 4.5,
                    'ARG': 5.0,
                    'C-': 4.5,
                    'N+': 4.5}

    all_atoms = group.atom.conformation_container.get_non_hydrogen_atoms()
    if residue.resName in version.desolvationRadii:
        local_cutoff = version.desolvationRadii[residue.resName]
    else:
        local_cutoff = 0.00
    residue.Nmass = 0
    residue.Nlocl = 0

    for atom in all_atoms:
        if atom.resNumb != group.atom.resNumb or atom.chainID != group.atom.chainID:
            dX = atom.x - residue.x
            dY = atom.y - residue.y
            dZ = atom.z - residue.z
            distance = math.sqrt(dX*dX + dY*dY + dZ*dZ)
            if distance < local_cutoff:
                group.Nlocl += 1
            if distance < parameters.buried_cutoff:
                group.Nmass += 1
        if residue.Nmass > 400:
            group.location = "BURIED "
        else:
            group.location = "SURFACE"
        group.Emass = group.charge * parameters.desolvationPrefactor * max(0.00, group.Nmass-parameters.desolvationAllowance)
        group.Elocl = group.charge * parameters.desolvationLocal * group.Nlocl
        # Buried ratio - new feature in propka3.0
        # Note, there will be an unforseen problem: e.g. if one residue has Nmass > Nmax and
        # the other Nmass < Nmax, the Npair will not be Nmass1 + Nmass2!
        residue.buried = calculateWeight(residue.Nmass)

        return 0.00, 0.00, 0.00, 0.00







def calculate_scale_factor(parameters, weight):
    scale_factor = 1.0 - (1.0 - parameters.desolvationSurfaceScalingFactor)*(1.0 - weight)
    return scale_factor


def calculate_weight(parameters, Nmass):
    """
    calculating the atom-based desolvation weight
    """
    weight = float(Nmass - parameters.Nmin)/float(parameters.Nmax - parameters.Nmin)
    weight = min(1.0, weight)
    weight = max(0.0, weight)

    return weight



def squared_distance(atom1, atom2):
#    if atom1 in atom2.squared_distances:
#        return atom2.squared_distances[atom1]

    dx = atom2.x - atom1.x
    dy = atom2.y - atom1.y
    dz = atom2.z - atom1.z

    res = dx*dx+dy*dy+dz*dz

    return res


def distance(atom1, atom2):
    return math.sqrt(squared_distance(atom1,atom2))



def get_smallest_distance(atoms1, atoms2):
    res_distance =1e6
    res_atom1 = None
    res_atom2 = None

    for atom1 in atoms1:
        for atom2 in atoms2:
            dist = squared_distance(atom1, atom2)
            if dist < res_distance:
                res_distance = dist
                res_atom1 = atom1
                res_atom2 = atom2

    return [res_atom1, math.sqrt(res_distance), res_atom2]


def calculatePairWeight(parameters, Nmass1, Nmass2):
    """
    calculating the atom-pair based desolvation weight
    """
    Nmass = Nmass1 + Nmass2
    Nmin  = 2*parameters.Nmin
    Nmax  = 2*parameters.Nmax
    weight = float(Nmass - Nmin)/float(Nmax - Nmin)
    weight = min(1.0, weight)
    weight = max(0.0, weight)

    return weight


def HydrogenBondEnergy(distance, dpka_max, cutoff, f_angle=1.0):
    """
    returns a hydrogen-bond interaction pKa shift
    """
    if   distance < cutoff[0]:
        value = 1.00
    elif distance > cutoff[1]:
        value = 0.00
    else:
        value = 1.0-(distance-cutoff[0])/(cutoff[1]-cutoff[0])

    dpKa  = dpka_max*value*f_angle

    return abs(dpKa)



def AngleFactorX(atom1=None, atom2=None, atom3=None, center=None):
    """
    Calculates the distance and angle-factor from three atoms for back-bone interactions,
    IMPORTANT: you need to use atom1 to be the e.g. ASP atom if distance is reset at return: [O1 -- H2-N3]
    Also generalized to be able to be used for residue 'centers' for C=O COO interactions.
    """

    dX_32 = atom2.x - atom3.x
    dY_32 = atom2.y - atom3.y
    dZ_32 = atom2.z - atom3.z

    distance_23 = math.sqrt( dX_32*dX_32 + dY_32*dY_32 + dZ_32*dZ_32 )

    dX_32 = dX_32/distance_23
    dY_32 = dY_32/distance_23
    dZ_32 = dZ_32/distance_23

    if atom1 == None:
        dX_21 = center[0] - atom2.x
        dY_21 = center[1] - atom2.y
        dZ_21 = center[2] - atom2.z
    else:
        dX_21 = atom1.x - atom2.x
        dY_21 = atom1.y - atom2.y
        dZ_21 = atom1.z - atom2.z

    distance_12 = math.sqrt( dX_21*dX_21 + dY_21*dY_21 + dZ_21*dZ_21 )

    dX_21 = dX_21/distance_12
    dY_21 = dY_21/distance_12
    dZ_21 = dZ_21/distance_12

    f_angle = dX_21*dX_32 + dY_21*dY_32 + dZ_21*dZ_32


    return distance_12, f_angle, distance_23



def hydrogen_bond_interaction(group1, group2, version):

    # find the smallest distance between interacting atoms
    atoms1 = group1.get_interaction_atoms(group2)
    atoms2 = group2.get_interaction_atoms(group1)
    [closest_atom1, distance, closest_atom2] = propka.calculations.get_smallest_distance(atoms1, atoms2)

    if None in [closest_atom1, closest_atom2]:
        warning('Side chain interaction failed for %s and %s' % (group1.label, group2.label))
        return None

    # get the parameters
    [dpka_max, cutoff] = version.get_hydrogen_bond_parameters(closest_atom1,closest_atom2)

    if dpka_max==None or None in cutoff:
        return None

    # check that the closest atoms are close enough
    if distance >= cutoff[1]:
        return None

    # check that bond distance criteria is met
    bond_distance_too_short = group1.atom.is_atom_within_bond_distance(group2.atom,
                                                                       version.parameters.min_bond_distance_for_hydrogen_bonds,1)
    if bond_distance_too_short:
        return None

    # set the angle factor
    #
    #  ---closest_atom1/2
    #                     .
    #                      .
    #                       the_hydrogen---closest_atom2/1---
    f_angle = 1.0
    if group2.type in version.parameters.angular_dependent_sidechain_interactions:
        if closest_atom2.element == 'H':
            heavy_atom = closest_atom2.bonded_atoms[0]
            hydrogen   = closest_atom2
            distance, f_angle, nada = propka.calculations.AngleFactorX(closest_atom1, hydrogen, heavy_atom)
        else:
            # Either the structure is corrupt (no hydrogen), or the heavy atom is closer to
            # the titratable atom than the hydrogen. In either case we set the angle factor
            # to 0
            f_angle = 0.0

    elif group1.type in version.parameters.angular_dependent_sidechain_interactions:
        if closest_atom1.element == 'H':
            heavy_atom = closest_atom1.bonded_atoms[0]
            hydrogen   = closest_atom1
            distance, f_angle, nada = propka.calculations.AngleFactorX(closest_atom2, hydrogen, heavy_atom)
        else:
            # Either the structure is corrupt (no hydrogen), or the heavy atom is closer to
            # the titratable atom than the hydrogen. In either case we set the angle factor
            # to 0
            f_angle = 0.0

    weight = version.calculatePairWeight(group1.Nmass, group2.Nmass)

    exception, value = version.checkExceptions(group1, group2)
    #exception = False # circumventing exception
    if exception == True:
        """ do nothing, value should have been assigned """
            #info(" exception for %s %s %6.2lf" % (group1.label, group2.label, value))
    else:
        value = version.calculateSideChainEnergy(distance, dpka_max, cutoff, weight, f_angle)

    # info('distance',distance)
    # info('dpka_max',dpka_max)
    # info('cutoff',cutoff)
    # info('f_angle',f_angle)
    # info('weight',weight)
    # info('value',value)
    # info('===============================================')

    return value



def HydrogenBondEnergy(distance, dpka_max, cutoff, f_angle=1.0):
    """
    returns a hydrogen-bond interaction pKa shift
    """
    if   distance < cutoff[0]:
        value = 1.00
    elif distance > cutoff[1]:
        value = 0.00
    else:
        value = 1.0-(distance-cutoff[0])/(cutoff[1]-cutoff[0])

    dpKa  = dpka_max*value*f_angle

    return abs(dpKa)




def electrostatic_interaction(group1, group2, distance, version):

    # check if we should do coulomb interaction at all
    if not version.checkCoulombPair(group1, group2, distance):
        return None

    weight = version.calculatePairWeight(group1.Nmass, group2.Nmass)
    value  = version.calculateCoulombEnergy(distance, weight)

    return value


def checkCoulombPair(parameters, group1, group2, distance):
    """
    Checks if this Coulomb interaction should be done - a propka2.0 hack
    """
    Npair = group1.Nmass + group2.Nmass
    do_coulomb = True

    # check if both groups are titratable (ions are taken care of in determinants::setIonDeterminants)
    if not (group1.titratable and group2.titratable):
        do_coulomb = False

    # check if the distance is not too big
    if distance > parameters.coulomb_cutoff2:
        do_coulomb = False

    # check that desolvation is ok
    if Npair < parameters.Nmin:
        do_coulomb = False

    return do_coulomb


def CoulombEnergy(distance, weight, parameters):
    """
    calculates the Coulomb interaction pKa shift based on Coulombs law
    eps = 60.0 for the moment; to be scaled with 'weight'
    """
    #diel = 80.0 -  60.0*weight

    diel = 160 - (160 -30)*weight
    R = max(distance, parameters.coulomb_cutoff1)
    scale = (R - parameters.coulomb_cutoff2)/(parameters.coulomb_cutoff1 -parameters.coulomb_cutoff2)
    scale = max(0.0, scale)
    scale = min(1.0, scale)

    dpka = 244.12/(diel*R) *scale

    return abs(dpka)



def BackBoneReorganization(parameters, conformation):
    """
    adding test stuff
    """
    titratable_groups = conformation.get_backbone_reorganisation_groups()
    BBC_groups = conformation.get_backbone_CO_groups()

    for titratable_group in titratable_groups:
        weight = titratable_group.buried
        dpKa = 0.00
        for BBC_group in BBC_groups:
            center = [titratable_group.x, titratable_group.y, titratable_group.z]
            distance, f_angle, nada = AngleFactorX(atom2=BBC_group.get_interaction_atoms(titratable_group)[0],
                                                   atom3=BBC_group.atom,
                                                   center=center)
            if distance <  6.0 and f_angle > 0.001:
                value = 1.0-(distance-3.0)/(6.0-3.0)
                dpKa += 0.80*min(1.0, value)

        titratable_group.Elocl = dpKa*weight
    return


#
# Exception methods
#

def checkExceptions(version, group1, group2):
    """
    checks for exceptions for this version - using defaults
    """
    resType1 = group1.type
    resType2 = group2.type

    if   (resType1 == "COO" and resType2 == "ARG"):
        exception, value = checkCooArgException(group1, group2, version)
    elif (resType1 == "ARG" and resType2 == "COO"):
        exception, value = checkCooArgException(group2, group1, version)
    elif (resType1 == "COO" and resType2 == "COO"):
        exception, value = checkCooCooException(group1, group2, version)
    elif (resType1 == "CYS" and resType2 == "CYS"):
        exception, value = checkCysCysException(group1, group2, version)
    elif (resType1 == "COO" and resType2 == "HIS") or \
            (resType1 == "HIS" and resType2 == "COO"):
        exception, value = checkCooHisException(group1, group2, version)
    elif (resType1 == "OCO" and resType2 == "HIS") or \
            (resType1 == "HIS" and resType2 == "OCO"):
        exception, value = checkOcoHisException(group1, group2, version)
    elif (resType1 == "CYS" and resType2 == "HIS") or \
            (resType1 == "HIS" and resType2 == "CYS"):
        exception, value = checkCysHisException(group1, group2, version)
    else:
        # do nothing, no exception for this pair
        exception = False; value = None

    return exception, value



def checkCooArgException(group_coo, group_arg, version):
    """
    checking Coo-Arg exception: uses the two shortes unique distances (involving 2+2 atoms)
    """

    str = "xxx"
    exception = True
    value_tot = 0.00

    #dpka_max = parameters.sidechain_interaction.get_value(group_coo.type,group_arg.type)
    #cutoff   = parameters.sidechain_cutoffs.get_value(group_coo.type,group_arg.type)

    # needs to be this way since you want to find shortest distance first
    #info("--- exception for %s %s ---" % (group_coo.label, group_arg.label))
    atoms_coo = []
    atoms_coo.extend(group_coo.get_interaction_atoms(group_arg))
    atoms_arg = []
    atoms_arg.extend(group_arg.get_interaction_atoms(group_coo))


    for iter in ["shortest", "runner-up"]:
        # find the closest interaction pair
        [closest_coo_atom, distance, closest_arg_atom] = get_smallest_distance(atoms_coo, atoms_arg)
        [dpka_max, cutoff] = version.get_hydrogen_bond_parameters(closest_coo_atom,closest_arg_atom)
     # calculate and sum up interaction energy
        f_angle = 1.00
        if group_arg.type in version.parameters.angular_dependent_sidechain_interactions:
            atom3 = closest_arg_atom.bonded_atoms[0]
            distance, f_angle, nada = AngleFactorX(closest_coo_atom, closest_arg_atom, atom3)

        value = HydrogenBondEnergy(distance, dpka_max, cutoff, f_angle)
        #info(iter, closest_coo_atom, closest_arg_atom,distance,value)
        value_tot += value
        # remove closest atoms before we attemp to find the runner-up pair
        atoms_coo.remove(closest_coo_atom)
        atoms_arg.remove(closest_arg_atom)


    return exception, value_tot


def checkCooCooException(group1, group2, version):
    """
    checking Coo-Coo hydrogen-bond exception
    """
    exception = True
    [closest_atom1, distance, closest_atom2] = get_smallest_distance(group1.get_interaction_atoms(group2),
                                                                     group2.get_interaction_atoms(group1))

    #dpka_max = parameters.sidechain_interaction.get_value(group1.type,group2.type)
    #cutoff   = parameters.sidechain_cutoffs.get_value(group1.type,group2.type)
    [dpka_max, cutoff] = version.get_hydrogen_bond_parameters(closest_atom1,closest_atom2)
    f_angle = 1.00
    value = HydrogenBondEnergy(distance, dpka_max, cutoff, f_angle)
    weight = calculatePairWeight(version.parameters, group1.Nmass, group2.Nmass)
    value = value * (1.0 + weight)

    return exception, value



def checkCooHisException(group1, group2, version):
    """
    checking Coo-His exception
    """
    exception = False
    if checkBuried(group1.Nmass, group2.Nmass):
        exception = True

    return exception, version.parameters.COO_HIS_exception


def checkOcoHisException(group1, group2, version):
    """
    checking Coo-His exception
    """
    exception = False
    if checkBuried(group1.Nmass, group2.Nmass):
        exception = True

    return exception, version.parameters.OCO_HIS_exception


def checkCysHisException(group1, group2, version):
    """
    checking Cys-His exception
    """
    exception = False
    if checkBuried(group1.Nmass, group2.Nmass):
        exception = True

    return exception, version.parameters.CYS_HIS_exception


def checkCysCysException(group1, group2, version):
    """
    checking Cys-Cys exception
    """
    exception = False
    if checkBuried(group1.Nmass, group2.Nmass):
        exception = True

    return exception, version.parameters.CYS_CYS_exception





def checkBuried(Nmass1, Nmass2):
    """
    returns True if an interaction is buried
    """

    if (Nmass1 + Nmass2 <= 900) and (Nmass1 <= 400 or Nmass2 <= 400):
        return False
    else:
        return True
