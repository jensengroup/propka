#
# Container for molecular conformations
#

from __future__ import division
from __future__ import print_function

import propka.group, propka.determinants, propka.determinant, propka.ligand, propka.output, propka.coupled_groups, functools
from propka.lib import info, warning

class Conformation_container:
    def __init__(self, name='', parameters=None, molecular_container=None):
        self.molecular_container = molecular_container
        self.name=name
        self.parameters=parameters
        self.atoms = []
        self.groups = []
        self.chains = []
        self.current_iter_item = 0

        self.marvin_pkas_calculated = False
        self.non_covalently_coupled_groups = False

        return


    #
    # Group related methods
    #
    def extract_groups(self):
        """ Generates at list of molecular groups needed for calculating pKa values """
        for atom in self.get_non_hydrogen_atoms():
            # has this atom been checked for groups?
            if atom.groups_extracted == 0:
                group = propka.group.is_group(self.parameters, atom)
            else:
                group = atom.group
                # if the atom has been checked in a another conformation, check if it has a
                # group that should be used in this conformation as well
            if group:
                self.setup_and_add_group(group)

        return

    def additional_setup_when_reading_input_file(self):

        # if a group is coupled and we are reading a .propka_input file,
        # some more configuration might be needed

        # print coupling map
        map = propka.output.make_interaction_map('Covalent coupling map for %s'%self,
                                                 self.get_covalently_coupled_groups(),
                                                 lambda g1,g2: g1 in g2.covalently_coupled_groups)
        info(map)

        # check if we should set a common charge centre as well
        if self.parameters.common_charge_centre:
            self.set_common_charge_centres()

        return

    def set_common_charge_centres(self):
        for system in self.get_coupled_systems(self.get_covalently_coupled_groups(), propka.group.Group.get_covalently_coupled_groups):
            # make a list of the charge centre coordinates
            all_coordinates = list(map(lambda g: [g.x, g.y, g.z], system))
            # find the common charge center
            ccc = functools.reduce(lambda g1,g2: [g1[0]+g2[0], g1[1]+g2[1], g1[2]+g2[2]], all_coordinates)
            ccc = list(map(lambda c: c/len(system), ccc))
            # set the ccc for all coupled groups in this system
            for g in system:
                [g.x, g.y, g.z] = ccc
                g.common_charge_centre = True
        return




    def find_covalently_coupled_groups(self):
        """ Finds covalently coupled groups and sets common charge centres if needed """
        for group in self.get_titratable_groups():
            # Find covalently bonded groups
            bonded_groups = self.find_bonded_titratable_groups(group.atom, 1, group.atom)

            # couple groups
            for cg in bonded_groups:
                if cg in group.covalently_coupled_groups:
                    continue
                if cg.atom.sybyl_type == group.atom.sybyl_type:
                    group.couple_covalently(cg)

        # check if we should set a common charge centre as well
        if self.parameters.common_charge_centre:
            self.set_common_charge_centres()

        # print coupling map
        map = propka.output.make_interaction_map('Covalent coupling map for %s'%self,
                                                 #self.get_titratable_groups(),
                                                 self.get_covalently_coupled_groups(),
                                                 lambda g1,g2: g1 in g2.covalently_coupled_groups)
        info(map)


        return


    def find_non_covalently_coupled_groups(self, verbose=False):
        # check if non-covalent coupling has already been set up in an input file
        if len(list(filter(lambda g: len(g.non_covalently_coupled_groups)>0, self.get_titratable_groups())))>0:
            self.non_covalently_coupled_groups = True

        propka.coupled_groups.nccg.identify_non_covalently_coupled_groups(self,verbose=verbose)

        # re-do the check
        if len(list(filter(lambda g: len(g.non_covalently_coupled_groups)>0, self.get_titratable_groups())))>0:
            self.non_covalently_coupled_groups = True

        return


    def find_bonded_titratable_groups(self, atom, no_bonds, original_atom):
        res = set()
        for ba in atom.bonded_atoms:
            # skip the original atom
            if ba == original_atom:
                continue
            # check if this atom has a titratable group
            if ba.group and ba.group.titratable and no_bonds <= self.parameters.coupling_max_number_of_bonds:
                res.add(ba.group)
            # check for titratable groups bonded to this atom
            if no_bonds < self.parameters.coupling_max_number_of_bonds:
                res |= self.find_bonded_titratable_groups(ba,no_bonds+1, original_atom)

        return res


    def setup_and_add_group(self, group):
        """ Checks if we want to include this group in the calculations """

        # Is it recognized as a group at all?
        if not group:
            return

        # Other checks (include ligands, which chains etc.)

        # if all ok, accept the group
        self.init_group(group)
        self.groups.append(group)

    def init_group(self, group):
        """
        Initialize the given Group object.
        """
        # set up the group
        group.parameters=self.parameters
        group.setup()

        # If --titrate_only option is set, make non-specified residues un-titratable:
        titrate_only = self.molecular_container.options.titrate_only
        if titrate_only is not None:
            at = group.atom
            if not (at.chainID, at.resNumb, at.icode) in titrate_only:
                group.titratable = False
                if group.residue_type == 'CYS':
                    group.exclude_cys_from_results = True


    #
    # pka calculation methods
    #

    def calculate_pka(self, version, options):
        info('\nCalculating pKas for', self)

        # calculate desolvation
        for group in self.get_titratable_groups()+self.get_ions():
            version.calculate_desolvation(group)

        # calculate backbone interactions
        propka.determinants.setBackBoneDeterminants(self.get_titratable_groups(), self.get_backbone_groups(), version)

        # setting ion determinants
        propka.determinants.setIonDeterminants(self, version)

        # calculating the back-bone reorganization/desolvation term
        version.calculateBackBoneReorganization(self)

        # setting remaining non-iterative and iterative side-chain & Coulomb interaction determinants
        propka.determinants.setDeterminants(self.get_sidechain_groups(), version=version, options=options)

        # calculating the total pKa values
        for group in self.groups: group.calculate_total_pka()

        # take coupling effects into account
        penalised_labels = self.coupling_effects()

        if self.parameters.remove_penalised_group and len(penalised_labels)>0:
            info('Removing penalised groups!!!')

            for g in self.get_titratable_groups():
                g.remove_determinants(penalised_labels)

            # re-calculating the total pKa values
            for group in self.groups: group.calculate_total_pka()


        return


    def coupling_effects(self):
        #
        # Bases: The group with the highest pKa (the most stable one in the
        #        charged form) will be the first one to adopt a proton as pH
        #        is lowered and this group is allowed to titrate.
        #        The remaining groups are penalised
        #
        # Acids: The group with the highest pKa (the least stable one in the
        #        charged form) will be the last group to loose the proton as
        #        pH is raised and will be penalised.
        #        The remaining groups are allowed to titrate.
        #
        penalised_labels = []

        for all_groups in self.get_coupled_systems(self.get_covalently_coupled_groups(),
                                                   propka.group.Group.get_covalently_coupled_groups):

            # check if we should share determinants
            if self.parameters.shared_determinants:
                self.share_determinants(all_groups)

            # find the group that has the highest pKa value
            first_group = max(all_groups, key=lambda g:g.pka_value)

            # In case of acids
            if first_group.charge < 0:
                first_group.coupled_titrating_group = min(all_groups, key=lambda g:g.pka_value)
                penalised_labels.append(first_group.label) # group with the highest pKa is penalised

            # In case of bases
            else:
                for a in all_groups:
                    if a == first_group:
                        continue # group with the highest pKa is allowed to titrate...
                    a.coupled_titrating_group = first_group
                    penalised_labels.append(a.label) #... and the rest is penalised

        return penalised_labels


    def share_determinants(self, groups):

        # make a list of the determinants to share
        types = ['sidechain','backbone','coulomb']
        for type in types:
            # find maximum value for each determinant
            max_dets = {}
            for g in groups:
                for d in g.determinants[type]:
                    # update max dets
                    if d.group not in max_dets.keys():
                        max_dets[d.group] = d.value
                    else:
                        max_dets[d.group] = max(d.value, max_dets[d.group], key= lambda v: abs(v))

             # overwrite/add maximum value for each determinant
            for det_group in max_dets.keys():
                new_determinant = propka.determinant.Determinant(det_group, max_dets[det_group])
                for g in groups:
                    g.set_determinant(new_determinant,type)


        return


    def get_coupled_systems(self, groups, get_coupled_groups):
        """ This generator will yield one covalently coupled system at the time """
        groups = set(groups)
        while len(groups)>0:
            # extract a system of coupled groups ...
            system = set()
            self.get_a_coupled_system_of_groups(groups.pop(), system, get_coupled_groups)
            # ... and remove them from the list
            groups -= system

            yield system

        return


    def get_a_coupled_system_of_groups(self, new_group, coupled_groups, get_coupled_groups):
        coupled_groups.add(new_group)
        [self.get_a_coupled_system_of_groups(c, coupled_groups, get_coupled_groups) for c in get_coupled_groups(new_group) if c not in coupled_groups]
        return


    #
    # Energy/summary-related methods
    #
    def calculate_folding_energy(self, pH=None, reference=None):
        ddg = 0.0
        for group in self.groups:
            #info('Folding energy for %s at pH %f: %f'%(group,pH,group.calculate_folding_energy(self.parameters, pH=pH, reference=reference)))
            ddg += group.calculate_folding_energy(self.parameters, pH=pH, reference=reference)

        return ddg

    def calculate_charge(self, parmaeters, pH=None):
        unfolded = folded = 0.0
        for group in self.get_titratable_groups():
            unfolded += group.calculate_charge(parmaeters, pH=pH, state='unfolded')
            folded   += group.calculate_charge(parmaeters, pH=pH, state='folded')

        return unfolded,folded


    #
    # conformation/bookkeeping/atom methods
    #

    def get_backbone_groups(self):
        """ returns all backbone groups needed for the pKa calculations """
        return [group for group in self.groups if 'BB' in group.type]

    def get_sidechain_groups(self):
        """ returns all sidechain groups needed for the pKa calculations """
        return [group for group in self.groups if ('BB' not in group.type\
                    and not group.atom.cysteine_bridge)]

    def get_covalently_coupled_groups(self):
        return [g for g in self.groups if len(g.covalently_coupled_groups)>0]

    def get_non_covalently_coupled_groups(self):
        return [g for g in self.groups if len(g.non_covalently_coupled_groups)>0]

    def get_backbone_NH_groups(self):
        """ returns all NH backbone groups needed for the pKa calculations """
        return [group for group in self.groups if group.type == 'BBN']

    def get_backbone_CO_groups(self):
        """ returns all CO backbone groups needed for the pKa calculations """
        return [group for group in self.groups if group.type == 'BBC']

    def get_groups_in_residue(self, residue):
        return [group for group in self.groups if group.residue_type == residue]

    def get_titratable_groups(self):
        return [group for group in self.groups if group.titratable]

    def get_groups_for_calculations(self):
        """
        Returns a list of groups that should be included in results report.
        If --titrate_only option is specified, only residues that are titratable
        and are in that list are included; otherwise all titratable residues
        and CYS residues are included.
        """
        return [group for group in self.groups if group.use_in_calculations()]

    def get_acids(self):
        return [group for group in self.groups if (group.residue_type in self.parameters.acid_list
                                                   and not group.atom.cysteine_bridge)]

    def get_backbone_reorganisation_groups(self):
        return [group for group in self.groups if (group.residue_type in self.parameters.backbone_reorganisation_list
                                                   and not group.atom.cysteine_bridge)]

    def get_ions(self):
        return [group for group in self.groups if group.residue_type in self.parameters.ions.keys()]

    def get_group_names(self, list):
        return [group for group in self.groups if group.type in list]


    def get_ligand_atoms(self):
        return [atom for atom in self.atoms if atom.type=='hetatm']

    def get_heavy_ligand_atoms(self):
        return [atom for atom in self.atoms if atom.type=='hetatm' and atom.element != 'H']

    def get_chain(self,chain):
        return [atom for atom in self.atoms if atom.chainID != chain]


    def add_atom(self, atom):
        #info(self,'adding',atom)
        self.atoms.append(atom)
        if not atom.conformation_container:
            atom.conformation_container = self
        if not atom.molecular_container:
            atom.molecular_container = self.molecular_container

        # store chain id for bookkeeping
        if not atom.chainID in self.chains:
            self.chains.append(atom.chainID)

        return

    def copy_atom(self, atom):
        new_atom  = atom.makeCopy()
        self.atoms.append(new_atom)
        new_atom.conformation_container = self

        return

    def get_non_hydrogen_atoms(self):
        return [atom for atom in self.atoms if atom.element!='H']


    def top_up(self, other):
        """ Tops up self with all atoms found in other but not in self """
        my_residue_labels = { a.residue_label for a in self.atoms }
        for atom in other.atoms:
            if not atom.residue_label in my_residue_labels:
                self.copy_atom(atom)
        return

    def find_group(self, group):
        for g in self.groups:
            if g.atom.residue_label == group.atom.residue_label:
                if g.type == group.type:
                    return g
        return False


    def set_ligand_atom_names(self):
        for atom in self.get_ligand_atoms():
            propka.ligand.assign_sybyl_type(atom)
        return



    def __str__(self):
        return'Conformation container %s with %d atoms and %d groups'%(self.name,len(self),len(self.groups))

    def __len__(self):
        return len(self.atoms)

    def sort_atoms(self):
        # sort the atoms ...
        self.atoms.sort(key=self.sort_atoms_key)
        # ... and re-number them
        for i in range(len(self.atoms)):
            self.atoms[i].numb = i+1

        return

    def sort_atoms_key(self, atom):
        key = ord(atom.chainID)*1e7
        key += atom.resNumb*1000
        if len(atom.name) > len(atom.element):
            key += ord(atom.name[len(atom.element)])
            #info(atom,ord(atom.name[len(atom.element)]), '|%s||%s|'%(atom.name,atom.element))
        return key
