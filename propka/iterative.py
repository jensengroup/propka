
from __future__ import division
from __future__ import print_function

import math, time

import propka.lib as lib
from propka.determinant import Determinant
import propka.calculations
from propka.lib import info, warning, debug

# Some library functions for the interative pKa determinants


def addtoDeterminantList(group1, group2, distance, iterative_interactions, version):
    """
    Adds 'iterative determinants' to list ..., [[R1, R2], [side-chain, coulomb], [A1, A2]], ...
    Note, the sign is determined when the interaction is added to the iterative object!
    Note, distance < coulomb_cutoff here
    """

    hbond_value   = version.hydrogen_bond_interaction(group1, group2)
    coulomb_value = version.electrostatic_interaction(group1, group2, distance)

    # adding the interaction to 'iterative_interactions'
    if hbond_value or coulomb_value:
        pair         = [group1, group2]

        values       = [hbond_value, coulomb_value]
        while None in values:
            values[values.index(None)] = 0.0

        annihilation = [0., 0.]
        interaction  = [pair, values, annihilation]
        iterative_interactions.append(interaction)

    return


def addIterativeAcidPair(object1, object2, interaction):
    """
    Adding the Coulomb 'iterative' interaction (an acid pair):
    the higher pKa is raised  with QQ+HB
    the lower  pKa is lowered with HB
    """
    values       = interaction[1]
    annihilation = interaction[2]
    hbond_value   = values[0]
    coulomb_value = values[1]
    diff = coulomb_value + 2*hbond_value
    comp1 = object1.pKa_old + annihilation[0] + diff
    comp2 = object2.pKa_old + annihilation[1] + diff
    annihilation[0] = 0.
    annihilation[1] = 0.
    if comp1 > comp2:
        # side-chain
        determinant = [object2,  hbond_value]
        object1.determinants['sidechain'].append(determinant)
        determinant = [object1, -hbond_value]
        object2.determinants['sidechain'].append(determinant)
        # Coulomb
        determinant = [object2, coulomb_value]
        object1.determinants['coulomb'].append(determinant)
        annihilation[0] = -diff
    else:
        # side-chain
        determinant = [object1,  hbond_value]
        object2.determinants['sidechain'].append(determinant)
        determinant = [object2, -hbond_value]
        object1.determinants['sidechain'].append(determinant)
        # Coulomb
        determinant = [object1, coulomb_value]
        object2.determinants['coulomb'].append(determinant)
        annihilation[1] = -diff


def addIterativeBasePair(object1, object2, interaction):
    """
    Adding the Coulomb 'iterative' interaction (a base pair):
    the lower pKa is lowered
    """
    values       = interaction[1]
    annihilation = interaction[2]
    hbond_value   = values[0]
    coulomb_value = values[1]
    diff = coulomb_value + 2*hbond_value
    diff = -diff
    comp1 = object1.pKa_old + annihilation[0] + diff
    comp2 = object2.pKa_old + annihilation[1] + diff
    annihilation[0] = 0.
    annihilation[1] = 0.
    if comp1 < comp2:
        # side-chain
        determinant = [object2, -hbond_value]
        object1.determinants['sidechain'].append(determinant)
        determinant = [object1,  hbond_value]
        object2.determinants['sidechain'].append(determinant)
        # Coulomb
        determinant = [object2, -coulomb_value]
        object1.determinants['coulomb'].append(determinant)
        annihilation[0] = -diff
    else:
        # side-chain
        determinant = [object1, -hbond_value]
        object2.determinants['sidechain'].append(determinant)
        determinant = [object2,  hbond_value]
        object1.determinants['sidechain'].append(determinant)
        # Coulomb
        determinant = [object1, -coulomb_value]
        object2.determinants['coulomb'].append(determinant)
        annihilation[1] = -diff


def addIterativeIonPair(object1, object2, interaction, version):
    """
    Adding the Coulomb 'iterative' interaction (an acid-base pair):
    the pKa of the acid is lowered & the pKa of the base is raised
    """
    values       = interaction[1]
    annihilation = interaction[2]
    hbond_value   = values[0]
    coulomb_value = values[1]
    Q1 = object1.Q
    Q2 = object2.Q
    comp1 = object1.pKa_old + annihilation[0] + Q1*coulomb_value
    comp2 = object2.pKa_old + annihilation[1] + Q2*coulomb_value
    if object1.resName not in version.parameters.exclude_sidechain_interactions:
        comp1 += Q1*hbond_value
    if object2.resName not in version.parameters.exclude_sidechain_interactions:
        comp2 += Q2*hbond_value

    if   Q1 == -1.0 and comp1 < comp2:
      add_term = True  # pKa(acid) < pKa(base)
    elif Q1 ==  1.0 and comp1 > comp2:
      add_term = True  # pKa(base) > pKa(acid)
    else:
      add_term = False

    annihilation[0] = 0.00
    annihilation[1] = 0.00

    if add_term == True:

      # Coulomb
      if coulomb_value > 0.005:
        # residue1
        interaction = [object2, Q1*coulomb_value]
        annihilation[0] += -Q1*coulomb_value
        object1.determinants['coulomb'].append(interaction)
        # residue2
        interaction = [object1, Q2*coulomb_value]
        annihilation[1] += -Q2*coulomb_value
        object2.determinants['coulomb'].append(interaction)

      # Side-chain
      if hbond_value > 0.005:
        # residue1
        if object1.resName not in version.parameters.exclude_sidechain_interactions:
          interaction = [object2, Q1*hbond_value]
          annihilation[0] += -Q1*hbond_value
          object1.determinants['sidechain'].append(interaction)
        # residue2
        if object2.resName not in version.parameters.exclude_sidechain_interactions:
          interaction = [object1, Q2*hbond_value]
          annihilation[1] += -Q2*hbond_value
          object2.determinants['sidechain'].append(interaction)


def addDeterminants(iterative_interactions, version, options=None):
    """
    The iterative pKa scheme. Later it is all added in 'calculateTotalPKA'
    """
    # --- setup ---
    iteratives = []
    done_group = []

    # creating iterative objects with references to their real group counterparts
    for interaction in iterative_interactions:
        pair = interaction[0]
        for group in pair:
            if group in done_group:
                #print "done already"
                """ do nothing - already have an iterative object for this group """
            else:
                newIterative = Iterative(group)
                iteratives.append(newIterative)
                done_group.append(group)

    # Initialize iterative scheme
    debug("\n   --- pKa iterations (%d groups, %d interactions) ---" %
          (len(iteratives), len(iterative_interactions)))
    converged = False
    iteration = 0
    # set non-iterative pka values as first step
    for itres in iteratives:
      itres.pKa_iter.append(itres.pKa_NonIterative)


    # --- starting pKa iterations ---
    while converged == False:

      # initialize pKa_new
      iteration += 1
      for itres in iteratives:
        itres.determinants = {'sidechain':[],'backbone':[],'coulomb':[]}
        itres.pKa_new = itres.pKa_NonIterative


      # Adding interactions to temporary determinant container
      for interaction in iterative_interactions:
        pair   = interaction[0]
        values = interaction[1]
        annihilation = interaction[2]
        #print "len(interaction) = %d" % (len(interaction))
        object1, object2 = findIterative(pair, iteratives)
        Q1 = object1.Q
        Q2 = object2.Q
        if   Q1 < 0.0 and Q2 < 0.0:
            """ both are acids """
            addIterativeAcidPair(object1, object2, interaction)
        elif Q1 > 0.0 and Q2 > 0.0:
            """ both are bases """
            addIterativeBasePair(object1, object2, interaction)
        else:
            """ one of each """
            addIterativeIonPair(object1, object2, interaction, version)


      # Calculating pKa_new values
      for itres in iteratives:
        for type in ['sidechain','backbone','coulomb']:
          for determinant in itres.determinants[type]:
            itres.pKa_new += determinant[1]

      # Check convergence
      converged = True
      for itres in iteratives:
        if itres.pKa_new == itres.pKa_old:
          itres.converged = True
        else:
          itres.converged = False
          converged = False

      # reset pKa_old & storing pKa_new in pKa_iter
      for itres in iteratives:
        itres.pKa_old = itres.pKa_new
        itres.pKa_iter.append(itres.pKa_new)

      if iteration == 10:
          info("did not converge in %d iterations" % (iteration))
          break

    # --- Iterations finished ---

    # printing pKa iterations
    # formerly was conditioned on if options.verbosity >= 2 - now unnecessary
    str = "%12s" % (" ")
    for index in range(0, iteration+1 ):
        str += "%8d" % (index)
    debug(str)
    for itres in iteratives:
        str  = "%s   " % (itres.label)
        for pKa in itres.pKa_iter:
          str += "%8.2lf" % (pKa)
        if itres.converged == False:
          str += " *"
        debug(str)

    # creating real determinants and adding them to group object
    for itres in iteratives:
        for type in ['sidechain','backbone','coulomb']:
            for interaction in itres.determinants[type]:
                #info('done',itres.group.label,interaction[0],interaction[1])
                value = interaction[1]
                if value > 0.005 or value < -0.005:
                    g = interaction[0]
                    newDeterminant = Determinant(g, value)
                    itres.group.determinants[type].append(newDeterminant)



def findIterative(pair, iteratives):
    """
    Function to find the two 'iteratives' that corresponds to the groups in 'pair'
    """
    for iterative in iteratives:
        if iterative.group == pair[0]:
            iterative0 = iterative
        elif iterative.group == pair[1]:
            iterative1 = iterative

    return iterative0, iterative1



class Iterative:
    """
        Iterative class - pKa values and references of iterative groups
        Note, this class has a fake determinant list, true determinants are
              made after the iterations are finished.
    """

    def __init__(self, group):
        """
        Contructer of the iterative object
        """

        #print "creating 'iterative object' for %s" % (group.label)

        self.label    = group.label
        self.atom     = group.atom
        self.resName  = group.residue_type
        self.Q        = group.charge
        self.pKa_old  = None
        self.pKa_new  = None
        self.pKa_iter = []
        self.pKa_NonIterative = 0.00
        self.determinants = {'sidechain':[],'backbone':[],'coulomb':[]}
        self.group = group
        self.converged = True

        # Calculate the Non-Iterative part of pKa from the group object
        # Side chain
        side_chain = 0.00
        for determinant in group.determinants['sidechain']:
            value = determinant.value
            side_chain += value

        # Back bone
        back_bone  = 0.00
        for determinant in group.determinants['backbone']:
            value = determinant.value
            back_bone  += value

        # Coulomb
        coulomb    = 0.00
        for determinant in group.determinants['coulomb']:
            value = determinant.value
            coulomb    += value

        self.pKa_NonIterative  = group.model_pka
        self.pKa_NonIterative += group.Emass
        self.pKa_NonIterative += group.Elocl
        self.pKa_NonIterative += side_chain
        self.pKa_NonIterative += back_bone
        self.pKa_NonIterative += coulomb

        self.pKa_old = self.pKa_NonIterative


    def __eq__(self, other):
        """
        Check if two groups should be considered identical
        """
        if self.atom.type == 'atom':
            # In case of protein atoms we trust the labels
            return self.label==other.label
        else:
            # For heterogene atoms we also need to check the residue number
            return self.label==other.label and self.atom.resNumb == other.atom.resNumb

    def __hash__(self):
        """ Needed together with __eq__ - otherwise we can't make sets of groups """
        return id(self)
