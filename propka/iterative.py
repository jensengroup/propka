"""
Working with Determinants
=========================

Iterative functions for pKa calculations. These appear to mostly
involve :class:`propka.determinant.Determinant` instances.

"""
import logging
from typing import Dict, Iterable, List, Optional, Sequence, Tuple
from propka.determinant import Determinant
from propka.group import Group
from propka.version import Version


_LOGGER = logging.getLogger(__name__)


# TODO - these are undocumented constants
UNK_MIN_VALUE = 0.005

Interaction = list


def add_to_determinant_list(group1: Group, group2: Group, distance: float,
                            iterative_interactions: List[Interaction], version: Version):
    """Add iterative determinantes to the list.

    [[R1, R2], [side-chain, coulomb], [A1, A2]], ...

    NOTE - sign is determined when the interaction is added to the iterative
      object!
    NOTE - distance < coulomb_cutoff here

    Args:
        group1:  first group in pair
        group2:  second group in pair
        distance:  distance between groups
        iterative_interactions:  interaction list to modify
        version:  version object
    
    """
    hbond_value = version.hydrogen_bond_interaction(group1, group2)
    coulomb_value = version.electrostatic_interaction(group1, group2, distance)
    # adding the interaction to 'iterative_interactions'
    if hbond_value or coulomb_value:
        pair = [group1, group2]
        values = [hbond_value, coulomb_value]
        while None in values:
            values[values.index(None)] = 0.0
        annihilation = [0., 0.]
        interaction = [pair, values, annihilation]
        iterative_interactions.append(interaction)


def add_iterative_acid_pair(object1: "Iterative", object2: "Iterative",
                            interaction: Interaction):
    """Add the Coulomb 'iterative' interaction (an acid pair).

    The higher pKa is raised  with QQ+HB
    The lower pKa is lowered with HB

    Args:
        object1:  first object in pair
        object2:  second object in pair
        interaction:  list with [values, annihilation]
    """
    values = interaction[1]
    annihilation = interaction[2]
    hbond_value = values[0]
    coulomb_value = values[1]
    diff = coulomb_value + 2*hbond_value
    comp1 = object1.pka_old + annihilation[0] + diff
    comp2 = object2.pka_old + annihilation[1] + diff
    annihilation[0] = 0.0
    annihilation[1] = 0.0
    if comp1 > comp2:
        # side-chain
        determinant = [object2, hbond_value]
        object1.determinants['sidechain'].append(determinant)
        determinant = [object1, -hbond_value]
        object2.determinants['sidechain'].append(determinant)
        # Coulomb
        determinant = [object2, coulomb_value]
        object1.determinants['coulomb'].append(determinant)
        annihilation[0] = -diff
    else:
        # side-chain
        determinant = [object1, hbond_value]
        object2.determinants['sidechain'].append(determinant)
        determinant = [object2, -hbond_value]
        object1.determinants['sidechain'].append(determinant)
        # Coulomb
        determinant = [object1, coulomb_value]
        object2.determinants['coulomb'].append(determinant)
        annihilation[1] = -diff


def add_iterative_base_pair(object1: "Iterative", object2: "Iterative",
                            interaction: Interaction):
    """Add the Coulomb 'iterative' interaction (a base pair).

    The lower pKa is lowered

    Args:
        object1:  first object in pair
        object2:  second object in pair
        interaction:  list with [values, annihilation]
    """
    values = interaction[1]
    annihilation = interaction[2]
    hbond_value = values[0]
    coulomb_value = values[1]
    diff = coulomb_value + 2*hbond_value
    diff = -diff
    comp1 = object1.pka_old + annihilation[0] + diff
    comp2 = object2.pka_old + annihilation[1] + diff
    annihilation[0] = 0.0
    annihilation[1] = 0.0
    if comp1 < comp2:
        # side-chain
        determinant = [object2, -hbond_value]
        object1.determinants['sidechain'].append(determinant)
        determinant = [object1, hbond_value]
        object2.determinants['sidechain'].append(determinant)
        # Coulomb
        determinant = [object2, -coulomb_value]
        object1.determinants['coulomb'].append(determinant)
        annihilation[0] = -diff
    else:
        # side-chain
        determinant = [object1, -hbond_value]
        object2.determinants['sidechain'].append(determinant)
        determinant = [object2, hbond_value]
        object1.determinants['sidechain'].append(determinant)
        # Coulomb
        determinant = [object1, -coulomb_value]
        object2.determinants['coulomb'].append(determinant)
        annihilation[1] = -diff


def add_iterative_ion_pair(object1: "Iterative", object2: "Iterative",
                           interaction: Interaction, version: Version):
    """Add the Coulomb 'iterative' interaction (an acid-base pair)

    the pKa of the acid is lowered & the pKa of the base is raised

    Args:
        object1:  first object in pair
        object2:  second object in pair
        interaction:  list with [values, annihilation]
        version:  version object
    """
    values = interaction[1]
    annihilation = interaction[2]
    hbond_value = values[0]
    coulomb_value = values[1]
    q1 = object1.q
    q2 = object2.q
    comp1 = object1.pka_old + annihilation[0] + q1*coulomb_value
    comp2 = object2.pka_old + annihilation[1] + q2*coulomb_value
    if (object1.res_name
            not in version.parameters.exclude_sidechain_interactions):
        comp1 += q1*hbond_value
    if (object2.res_name
            not in version.parameters.exclude_sidechain_interactions):
        comp2 += q2*hbond_value
    if q1 == -1.0 and comp1 < comp2:
        # pKa(acid) < pKa(base)
        add_term = True
    elif q1 == 1.0 and comp1 > comp2:
        # pKa(base) > pKa(acid)
        add_term = True
    else:
        add_term = False
    annihilation[0] = 0.00
    annihilation[1] = 0.00
    if add_term:
        # Coulomb
        if coulomb_value > UNK_MIN_VALUE:
            # residue1
            interaction = [object2, q1*coulomb_value]
            annihilation[0] += -q1*coulomb_value
            object1.determinants['coulomb'].append(interaction)
            # residue2
            interaction = [object1, q2*coulomb_value]
            annihilation[1] += -q2*coulomb_value
            object2.determinants['coulomb'].append(interaction)
        # Side-chain
        if hbond_value > UNK_MIN_VALUE:
            # residue1
            if (object1.res_name
                    not in version.parameters.exclude_sidechain_interactions):
                interaction = [object2, q1*hbond_value]
                annihilation[0] += -q1*hbond_value
                object1.determinants['sidechain'].append(interaction)
            # residue2
            if (object2.res_name
                    not in version.parameters.exclude_sidechain_interactions):
                interaction = [object1, q2*hbond_value]
                annihilation[1] += -q2*hbond_value
                object2.determinants['sidechain'].append(interaction)


def add_determinants(iterative_interactions: List[Interaction], version: Version, _=None):
    """Add determinants iteratively.

    The iterative pKa scheme. Later it is all added in 'calculateTotalPKA'

    Args:
        iterative_interactions:  list of iterative interactions
        version:  version object
        _:  options object
    """
    # --- setup ---
    iteratives: List[Iterative] = []
    done_group = []
    # create iterative objects with references to their real group counterparts
    for interaction in iterative_interactions:
        pair = interaction[0]
        for group in pair:
            if group in done_group:
                # do nothing - already have an iterative object for this group
                pass
            else:
                new_iterative = Iterative(group)
                iteratives.append(new_iterative)
                done_group.append(group)
    # Initialize iterative scheme
    _LOGGER.debug(
        "\n   --- pKa iterations ({0:d} groups, {1:d} interactions) "
        "---".format(
            len(iteratives), len(iterative_interactions)
            )
    )
    converged = False
    iteration = 0
    # set non-iterative pka values as first step
    for iter_ in iteratives:
        iter_.pka_iter.append(iter_.pka_noniterative)
    # --- starting pKa iterations ---
    while not converged:
        # initialize pka_new
        iteration += 1
        for itres in iteratives:
            itres.determinants = {'sidechain': [], 'backbone': [],
                                  'coulomb': []}
            itres.pka_new = itres.pka_noniterative
        # Adding interactions to temporary determinant container
        for interaction in iterative_interactions:
            pair = interaction[0]
            object1, object2 = find_iterative(pair, iteratives)
            q1 = object1.q
            q2 = object2.q
            if q1 < 0.0 and q2 < 0.0:
                # both are acids
                add_iterative_acid_pair(object1, object2, interaction)
            elif q1 > 0.0 and q2 > 0.0:
                # both are bases
                add_iterative_base_pair(object1, object2, interaction)
            else:
                # one of each
                add_iterative_ion_pair(object1, object2, interaction, version)
        # Calculating pka_new values
        for itres in iteratives:
            for type_ in ['sidechain', 'backbone', 'coulomb']:
                for determinant in itres.determinants[type_]:
                    itres.pka_new += determinant[1]

        # Check convergence
        converged = True
        for itres in iteratives:
            if itres.pka_new == itres.pka_old:
                itres.converged = True
            else:
                itres.converged = False
                converged = False

        # reset pka_old & storing pka_new in pka_iter
        for itres in iteratives:
            assert itres.pka_new is not None
            itres.pka_old = itres.pka_new
            itres.pka_iter.append(itres.pka_new)

        if iteration == 10:
            _LOGGER.info(
                "did not converge in {0:d} iterations".format(iteration)
            )
            break
    # printing pKa iterations
    # formerly was conditioned on if options.verbosity >= 2 - now unnecessary
    str_ = '            '
    for index in range(iteration+1):
        str_ += "{0:>8d}".format(index)
    _LOGGER.debug(str_)
    for itres in iteratives:
        str_ = "{0:s}   ".format(itres.label)
        for pka in itres.pka_iter:
            str_ += "{0:>8.2f}".format(pka)
        if not itres.converged:
            str_ += " *"
        _LOGGER.debug(str_)
    # creating real determinants and adding them to group object
    for itres in iteratives:
        for type_ in ['sidechain', 'backbone', 'coulomb']:
            for interaction in itres.determinants[type_]:
                value: float = interaction[1]
                if value > UNK_MIN_VALUE or value < -UNK_MIN_VALUE:
                    group = interaction[0]
                    new_det = Determinant(group, value)
                    itres.group.determinants[type_].append(new_det)


def find_iterative(
    pair: Sequence[Group],
    iteratives: Iterable["Iterative"],
) -> Tuple["Iterative", "Iterative"]:
    """Find the 'iteratives' that correspond to the groups in 'pair'.

    Args:
        pair:  groups to match
        iteratives:  list of iteratives to search
    Returns:
        1. first matched iterative
        2. second matched iterative
    """
    iterative0: Optional[Iterative] = None
    iterative1: Optional[Iterative] = None
    for iterative in iteratives:
        if iterative.group == pair[0]:
            iterative0 = iterative
        elif iterative.group == pair[1]:
            iterative1 = iterative
    if iterative0 is None or iterative1 is None:
        raise LookupError("iteratives not found")
    return iterative0, iterative1


class Iterative:
    """Iterative class - pKa values and references of iterative groups.

    NOTE - this class has a fake determinant list, true determinants are made
    after the iterations are finished.
    """

    def __init__(self, group: Group):
        """Initialize object with group.

        Args:
            group:  group to use for initialization.
        """
        self.label = group.label
        self.atom = group.atom
        self.res_name = group.residue_type
        self.q = group.charge
        self.pka_old: Optional[float] = None
        self.pka_new: Optional[float] = None
        self.pka_iter: List[float] = []
        self.pka_noniterative = 0.00
        self.determinants: Dict[str, list] = {
            'sidechain': [],
            'backbone': [],
            'coulomb': []
        }
        self.group = group
        self.converged = True
        # Calculate the Non-Iterative part of pKa from the group object
        # Side chain
        side_chain = 0.00
        for determinant in group.determinants['sidechain']:
            value = determinant.value
            side_chain += value
        # Back bone
        back_bone = 0.00
        for determinant in group.determinants['backbone']:
            value = determinant.value
            back_bone += value
        # Coulomb
        coulomb = 0.00
        for determinant in group.determinants['coulomb']:
            value = determinant.value
            coulomb += value
        self.pka_noniterative = group.model_pka
        self.pka_noniterative += group.energy_volume
        self.pka_noniterative += group.energy_local
        self.pka_noniterative += side_chain
        self.pka_noniterative += back_bone
        self.pka_noniterative += coulomb
        self.pka_old = self.pka_noniterative

    def __eq__(self, other) -> bool:
        """Needed to use objects in sets."""
        assert isinstance(other, (Iterative, Group)), type(other)
        if self.atom.type == 'atom':
            # In case of protein atoms we trust the labels
            return self.label == other.label
        else:
            # For heterogene atoms we also need to check the residue number
            return (
                self.label == other.label
                and self.atom.res_num == other.atom.res_num)

    def __hash__(self):
        """Needed to use objects in sets."""
        return id(self)
