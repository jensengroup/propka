"""
Energy calculations
===================

Energy calculations.

"""
import math
import logging
from typing import TYPE_CHECKING, Optional, Sequence

from propka.atom import Atom
from propka.parameters import Parameters

if TYPE_CHECKING:
    from propka.conformation_container import ConformationContainer
    from propka.group import Group
    from propka.version import Version

from propka.calculations import squared_distance, get_smallest_distance


_LOGGER = logging.getLogger(__name__)


# TODO - I have no idea what these constants mean so I labeled them "UNK_"
UNK_MIN_DISTANCE = 2.75
MIN_DISTANCE_4TH = math.pow(UNK_MIN_DISTANCE, 4)
UNK_DIELECTRIC1 = 160
UNK_DIELECTRIC2 = 30
UNK_PKA_SCALING1 = 244.12
UNK_BACKBONE_DISTANCE1 = 6.0
UNK_BACKBONE_DISTANCE2 = 3.0
UNK_FANGLE_MIN = 0.001
UNK_PKA_SCALING2 = 0.8
COMBINED_NUM_BURIED_MAX = 900
SEPARATE_NUM_BURIED_MAX = 400


def radial_volume_desolvation(parameters, group: "Group") -> None:
    """Calculate desolvation terms for group.

    Args:
        parameters:  parameters for desolvation calculation
        group:  group of atoms for calculation
    """
    assert group.atom.conformation_container is not None
    all_atoms = group.atom.conformation_container.get_non_hydrogen_atoms()
    volume = 0.0
    group.num_volume = 0
    min_dist_4th = MIN_DISTANCE_4TH
    for atom in all_atoms:
        # ignore atoms in the same residue
        if (atom.res_num == group.atom.res_num
                and atom.chain_id == group.atom.chain_id):
            continue
        sq_dist = squared_distance(group, atom)
        # desolvation
        if sq_dist < parameters.desolv_cutoff_squared:
            # use a default relative volume of 1.0 if the volume of the element
            # is not found in parameters
            # insert check for methyl groups
            if atom.element == 'C' and atom.name not in ['CA', 'C']:
                dvol = parameters.VanDerWaalsVolume['C4']
            else:
                dvol = parameters.VanDerWaalsVolume.get(atom.element, 1.0)
            dv_inc = dvol/max(min_dist_4th, sq_dist*sq_dist)
            volume += dv_inc
        # buried
        if sq_dist < parameters.buried_cutoff_squared:
            group.num_volume += 1
    group.buried = calculate_weight(parameters, group.num_volume)
    scale_factor = calculate_scale_factor(parameters, group.buried)
    volume_after_allowance = max(0.00, volume-parameters.desolvationAllowance)
    group.energy_volume = (
        group.charge * parameters.desolvationPrefactor
        * volume_after_allowance * scale_factor)


def calculate_scale_factor(parameters, weight: float) -> float:
    """Calculate desolvation scaling factor.

    Args:
        parameters:  parameters for desolvation calculation
        weight:  weight for scaling factor
    Returns:
        scaling factor
    """
    scale_factor = (
        1.0 - (1.0 - parameters.desolvationSurfaceScalingFactor)
        * (1.0 - weight)
    )
    return scale_factor


def calculate_weight(parameters: Parameters, num_volume: float) -> float:
    """Calculate the atom-based desolvation weight.

    TODO - figure out why a similar function exists in version.py

    Args:
        parameters:  parameters for desolvation calculation
        num_volume:  number of heavy atoms within desolvation calculation
                     volume
    Returns:
        desolvation weight
    """
    weight = (
        float(num_volume - parameters.Nmin)
        / float(parameters.Nmax - parameters.Nmin))
    weight = min(1.0, weight)
    weight = max(0.0, weight)
    return weight


def calculate_pair_weight(parameters: Parameters, num_volume1: int, num_volume2: int) -> float:
    """Calculate the atom-pair based desolvation weight.

    Args:
        num_volume1:  number of heavy atoms within first desolvation volume
        num_volume2:  number of heavy atoms within second desolvation volume
    Returns:
        desolvation weight
    """
    num_volume = num_volume1 + num_volume2
    num_min = 2*parameters.Nmin
    num_max = 2*parameters.Nmax
    weight = float(num_volume - num_min)/float(num_max - num_min)
    weight = min(1.0, weight)
    weight = max(0.0, weight)
    return weight


def hydrogen_bond_energy(dist, dpka_max: float, cutoffs, f_angle=1.0) -> float:
    """Calculate hydrogen-bond interaction pKa shift.

    Args:
        dist:  distance for hydrogen bond
        dpka_max:  maximum pKa value shift
        cutoffs:  array with max and min distance values
        f_angle:  angle scaling factor
    Returns:
        pKa shift value
    """
    if dist < cutoffs[0]:
        value = 1.00
    elif dist > cutoffs[1]:
        value = 0.00
    else:
        value = 1.0 - (dist - cutoffs[0])/(cutoffs[1] - cutoffs[0])
    dpka = dpka_max*value*f_angle
    return abs(dpka)


def angle_distance_factors(
        atom1: Optional[Atom] = None,
        atom2: Atom = None,  # type: ignore[assignment]
        atom3: Atom = None,  # type: ignore[assignment]
        center: Optional[Sequence[float]] = None):
    """Calculate distance and angle factors for three atoms for backbone
    interactions.

    NOTE - you need to use atom1 to be the e.g. ASP atom if distance is reset
           at return: [O1 -- H2-N3].

    Also generalized to be able to be used for residue 'centers' for C=O COO
    interactions.

    Args:
        atom1:  first atom for calculation (could be None)
        atom2:  second atom for calculation
        atom3:  third atom for calculation
        center:  center point between atoms 1 and 2
    Returns
        [distance factor between atoms 1 and 2,
         angle factor,
         distance factor between atoms 2 and 3]
    """
    # The angle factor
    #
    #  ---closest_atom1/2
    #                     .
    #                      .
    #                       the_hydrogen---closest_atom2/1---
    dx_32 = atom2.x - atom3.x
    dy_32 = atom2.y - atom3.y
    dz_32 = atom2.z - atom3.z
    dist_23 = math.sqrt(dx_32 * dx_32 + dy_32 * dy_32 + dz_32 * dz_32)
    dx_32 = dx_32/dist_23
    dy_32 = dy_32/dist_23
    dz_32 = dz_32/dist_23
    if atom1 is None:
        assert center is not None
        dx_21 = center[0] - atom2.x
        dy_21 = center[1] - atom2.y
        dz_21 = center[2] - atom2.z
    else:
        dx_21 = atom1.x - atom2.x
        dy_21 = atom1.y - atom2.y
        dz_21 = atom1.z - atom2.z
    dist_12 = math.sqrt(dx_21 * dx_21 + dy_21 * dy_21 + dz_21 * dz_21)
    dx_21 = dx_21/dist_12
    dy_21 = dy_21/dist_12
    dz_21 = dz_21/dist_12
    f_angle = dx_21*dx_32 + dy_21*dy_32 + dz_21*dz_32
    return dist_12, f_angle, dist_23


def hydrogen_bond_interaction(group1: "Group", group2: "Group", version: "Version"):
    """Calculate energy for hydrogen bond interactions between two groups.

    Args:
        group1:  first interacting group
        group2:  second interacting group
        version:  an object that contains version-specific parameters
    Returns:
        hydrogen bond interaction energy
    """
    # find the smallest distance between interacting atoms
    atoms1 = group1.get_interaction_atoms(group2)
    atoms2 = group2.get_interaction_atoms(group1)
    [closest_atom1, dist, closest_atom2] = get_smallest_distance(
        atoms1, atoms2
    )
    if closest_atom1 is None or closest_atom2 is None:
        _LOGGER.warning(
            'Side chain interaction failed for {0:s} and {1:s}'.format(
                group1.label, group2.label))
        return None
    # get the parameters
    [dpka_max, cutoff] = version.get_hydrogen_bond_parameters(closest_atom1,
                                                              closest_atom2)
    if (dpka_max is None) or (None in cutoff):
        return None
    # check that the closest atoms are close enough
    if dist >= cutoff[1]:
        return None
    # check that bond distance criteria is met
    min_hbond_dist = version.parameters.min_bond_distance_for_hydrogen_bonds
    if group1.atom.is_atom_within_bond_distance(
            group2.atom, min_hbond_dist, 1
            ):
        return None
    # set angle factor
    f_angle = 1.0
    if (
            group2.type in
            version.parameters.angular_dependent_sidechain_interactions
            ):
        if closest_atom2.element == 'H':
            heavy_atom = closest_atom2.bonded_atoms[0]
            hydrogen = closest_atom2
            dist, f_angle, _ = angle_distance_factors(closest_atom1, hydrogen,
                                                      heavy_atom)
        else:
            # Either the structure is corrupt (no hydrogen), or the heavy atom
            # is closer to the titratable atom than the hydrogen. In either
            # case, we set the angle factor to 0
            f_angle = 0.0
    elif (
            group1.type in
            version.parameters.angular_dependent_sidechain_interactions
            ):
        if closest_atom1.element == 'H':
            heavy_atom = closest_atom1.bonded_atoms[0]
            hydrogen = closest_atom1
            dist, f_angle, _ = angle_distance_factors(closest_atom2, hydrogen,
                                                      heavy_atom)
        else:
            # Either the structure is corrupt (no hydrogen), or the heavy atom
            # is closer to the titratable atom than the hydrogen. In either
            # case, we set the angle factor to 0
            f_angle = 0.0
    weight = version.calculate_pair_weight(
        group1.num_volume, group2.num_volume
    )
    exception, value = version.check_exceptions(group1, group2)
    if exception:
        # Do nothing, value should have been assigned.
        pass
    else:
        value = version.calculate_side_chain_energy(
            dist, dpka_max, cutoff, weight, f_angle)
    return value


def electrostatic_interaction(group1, group2, dist, version):
    """Calculate electrostatic interaction betwee two groups.

    Args:
        group1:  first interacting group
        group2:  second interacting group
        dist:  distance between groups
        version:  version-specific object with parameters and functions
    Returns:
        electrostatic interaction energy or None (if no interaction is
        appropriate)
    """
    # check if we should do coulomb interaction at all
    if not version.check_coulomb_pair(group1, group2, dist):
        return None
    weight = version.calculate_pair_weight(
        group1.num_volume, group2.num_volume
    )
    value = version.calculate_coulomb_energy(dist, weight)
    return value


def check_coulomb_pair(parameters: Parameters, group1: "Group", group2: "Group", dist: float) -> bool:
    """Checks if this Coulomb interaction should be done.

    NOTE - this is a propka2.0 hack
    TODO - figure out why a similar function exists in version.py

    Args:
        parameters:  parameters for Coulomb calculations
        group1:  first interacting group
        group2:  second interacting group
        dist:  distance between groups
    Returns:
        Boolean
    """
    num_volume = group1.num_volume + group2.num_volume
    do_coulomb = True
    # check if both groups are titratable (ions are taken care of in
    # determinants::set_ion_determinants)
    if not (group1.titratable and group2.titratable):
        do_coulomb = False
    # check if the distance is not too big
    if dist > parameters.coulomb_cutoff2:
        do_coulomb = False
    # check that desolvation is ok
    if num_volume < parameters.Nmin:
        do_coulomb = False
    return do_coulomb


def coulomb_energy(dist: float, weight: float, parameters) -> float:
    """Calculates the Coulomb interaction pKa shift based on Coulomb's law.

    Args:
        dist:  distance for electrostatic interaction
        weight:  scaling of dielectric constant
        parameters:  parameter object for calculation
    Returns:
        pKa shift
    """
    diel = UNK_DIELECTRIC1 - (UNK_DIELECTRIC1 - UNK_DIELECTRIC2)*weight
    dist = max(dist, parameters.coulomb_cutoff1)
    scale = (
        (dist - parameters.coulomb_cutoff2)
        / (parameters.coulomb_cutoff1 - parameters.coulomb_cutoff2))
    scale = max(0.0, scale)
    scale = min(1.0, scale)
    dpka = UNK_PKA_SCALING1/(diel*dist)*scale
    return abs(dpka)


def backbone_reorganization(_, conformation: "ConformationContainer") -> None:
    """Perform calculations related to backbone reorganizations.

    NOTE - this was described in the code as "adding test stuff"
    NOTE - this function does not appear to be used
    TODO - figure out why a similar function exists in version.py

    Args:
        _:  not used
        conformation:  specific molecule conformation
    """
    titratable_groups = conformation.get_backbone_reorganisation_groups()
    bbc_groups = conformation.get_backbone_co_groups()

    for titratable_group in titratable_groups:
        weight = titratable_group.buried
        dpka = 0.00
        for bbc_group in bbc_groups:
            center = [
                titratable_group.x, titratable_group.y, titratable_group.z
            ]
            atom2 = bbc_group.get_interaction_atoms(titratable_group)[0]
            dist, f_angle, _ = angle_distance_factors(atom2=atom2,
                                                      atom3=bbc_group.atom,
                                                      center=center)
            if dist < UNK_BACKBONE_DISTANCE1 and f_angle > UNK_FANGLE_MIN:
                value = (
                    1.0 - (dist-UNK_BACKBONE_DISTANCE2)
                    / (UNK_BACKBONE_DISTANCE1-UNK_BACKBONE_DISTANCE2))
                dpka += UNK_PKA_SCALING2 * min(1.0, value)
        titratable_group.energy_local = dpka*weight


def check_exceptions(version, group1, group2):
    """Checks for atypical behavior in interactions between two groups.
    Checks are made based on group type.

    TODO - figure out why a similar function exists in version.py

    Args:
        version:  version object
        group1:  first group for check
        group2:  second group for check
    Returns:
        1. Boolean indicating atypical behavior,
        2. value associated with atypical interaction (None if Boolean is
           False)
    """
    res_type1 = group1.type
    res_type2 = group2.type
    if (res_type1 == "COO") and (res_type2 == "ARG"):
        exception, value = check_coo_arg_exception(group1, group2, version)
    elif (res_type1 == "ARG") and (res_type2 == "COO"):
        exception, value = check_coo_arg_exception(group2, group1, version)
    elif (res_type1 == "COO") and (res_type2 == "COO"):
        exception, value = check_coo_coo_exception(group1, group2, version)
    elif (res_type1 == "CYS") and (res_type2 == "CYS"):
        exception, value = check_cys_cys_exception(group1, group2, version)
    elif ((res_type1 == "COO") and (res_type2 == "HIS")
          or (res_type1 == "HIS") and (res_type2 == "COO")):
        exception, value = check_coo_his_exception(group1, group2, version)
    elif ((res_type1 == "OCO") and (res_type2 == "HIS")
          or (res_type1 == "HIS") and (res_type2 == "OCO")):
        exception, value = check_oco_his_exception(group1, group2, version)
    elif ((res_type1 == "CYS") and (res_type2 == "HIS")
          or (res_type1 == "HIS") and (res_type2 == "CYS")):
        exception, value = check_cys_his_exception(group1, group2, version)
    else:
        # do nothing, no exception for this pair
        exception = False
        value = None
    return exception, value


def check_coo_arg_exception(group_coo, group_arg, version):
    """Check for COO-ARG interaction atypical behavior.

    Uses the two shortest unique distances (involving 2+2 atoms)

    Args:
        group_coo:  COO group
        group_arg:  ARG group
        version:  version object
    Returns:
        1. Boolean indicating atypical behavior,
        2. value associated with atypical interaction (None if Boolean is
           False)
    """
    exception = True
    value_tot = 0.00
    # needs to be this way since you want to find shortest distance first
    atoms_coo = []
    atoms_coo.extend(group_coo.get_interaction_atoms(group_arg))
    atoms_arg = []
    atoms_arg.extend(group_arg.get_interaction_atoms(group_coo))
    for _ in ["shortest", "runner-up"]:
        # find the closest interaction pair
        [closest_coo_atom, dist, closest_arg_atom] = get_smallest_distance(
            atoms_coo, atoms_arg
        )
        if closest_coo_atom is None:
            _LOGGER.warning(f"COO interaction atoms missing for {group_coo}")
            continue
        if closest_arg_atom is None:
            _LOGGER.warning(f"ARG interaction atoms missing for {group_arg}")
            continue
        [dpka_max, cutoff] = version.get_hydrogen_bond_parameters(
            closest_coo_atom, closest_arg_atom
        )
        # calculate and sum up interaction energy
        f_angle = 1.00
        if (
                group_arg.type in
                version.parameters.angular_dependent_sidechain_interactions
                ):
            atom3 = closest_arg_atom.bonded_atoms[0]
            dist, f_angle, _ = angle_distance_factors(closest_coo_atom,
                                                      closest_arg_atom,
                                                      atom3)
        value = hydrogen_bond_energy(dist, dpka_max, cutoff, f_angle)
        value_tot += value
        # remove closest atoms before we attemp to find the runner-up pair
        atoms_coo.remove(closest_coo_atom)
        atoms_arg.remove(closest_arg_atom)
    return exception, value_tot


def check_coo_coo_exception(group1, group2, version):
    """Check for COO-COO hydrogen-bond atypical interaction behavior.

    Args:
        group1:  first group for check
        group2:  second group for check
        version:  version object
    Returns:
        1. Boolean indicating atypical behavior,
        2. value associated with atypical interaction (None if Boolean is
           False)
    """
    exception = True
    interact_groups12 = group1.get_interaction_atoms(group2)
    interact_groups21 = group2.get_interaction_atoms(group1)
    [closest_atom1, dist, closest_atom2] = get_smallest_distance(
        interact_groups12, interact_groups21
    )
    [dpka_max, cutoff] = version.get_hydrogen_bond_parameters(
        closest_atom1, closest_atom2
    )
    f_angle = 1.00
    value = hydrogen_bond_energy(dist, dpka_max, cutoff, f_angle)
    weight = calculate_pair_weight(
        version.parameters, group1.num_volume, group2.num_volume
    )
    value = value * (1.0 + weight)
    return exception, value


def check_coo_his_exception(group1, group2, version):
    """Check for COO-HIS atypical interaction behavior

    Args:
        group1:  first group for check
        group2:  second group for check
        version:  version object
    Returns:
        1. Boolean indicating atypical behavior,
        2. value associated with atypical interaction (None if Boolean is
           False)
    """
    exception = False
    if check_buried(group1.num_volume, group2.num_volume):
        exception = True
    return exception, version.parameters.COO_HIS_exception


def check_oco_his_exception(group1, group2, version):
    """Check for OCO-HIS atypical interaction behavior

    Args:
        group1:  first group for check
        group2:  second group for check
        version:  version object
    Returns:
        1. Boolean indicating atypical behavior,
        2. value associated with atypical interaction (None if Boolean is
           False)
    """
    exception = False
    if check_buried(group1.num_volume, group2.num_volume):
        exception = True
    return exception, version.parameters.OCO_HIS_exception


def check_cys_his_exception(group1, group2, version):
    """Check for CYS-HIS atypical interaction behavior

    Args:
        group1:  first group for check
        group2:  second group for check
        version:  version object
    Returns:
        1. Boolean indicating atypical behavior,
        2. value associated with atypical interaction (None if Boolean is
           False)
    """
    exception = False
    if check_buried(group1.num_volume, group2.num_volume):
        exception = True
    return exception, version.parameters.CYS_HIS_exception


def check_cys_cys_exception(group1, group2, version):
    """Check for CYS-CYS atypical interaction behavior

    Args:
        group1:  first group for check
        group2:  second group for check
        version:  version object
    Returns:
        1. Boolean indicating atypical behavior,
        2. value associated with atypical interaction (None if Boolean is
           False)
    """
    exception = False
    if check_buried(group1.num_volume, group2.num_volume):
        exception = True
    return exception, version.parameters.CYS_CYS_exception


def check_buried(num_volume1, num_volume2):
    """Check to see if an interaction is buried

    Args:
        num_volume1:  number of buried heavy atoms in volume 1
        num_volume2:  number of buried heavy atoms in volume 2
    Returns:
        True if interaction is buried, False otherwise
    """
    if ((num_volume1 + num_volume2 <= COMBINED_NUM_BURIED_MAX)
            and (num_volume1 <= SEPARATE_NUM_BURIED_MAX
                 or num_volume2 <= SEPARATE_NUM_BURIED_MAX)):
        return False
    return True
