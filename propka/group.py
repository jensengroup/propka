#
# Class for storing groups important for propka calculations
#

from __future__ import division
from __future__ import print_function

import propka.ligand, propka.determinant, propka.ligand_pka_values, math, propka.protonate
from propka.lib import info, warning

my_protonator = propka.protonate.Protonate(verbose=False)

expected_atoms_acid_interactions = {
    'COO':{'O':2},
    'HIS':{'H':2, 'N':2},
    'CYS':{'S':1},
    'TYR':{'O':1},
    'LYS':{'N':1},
    'ARG':{'H':5, 'N':3},
    'ROH':{'O':1},
    'AMD':{'H':2, 'N':1},
    'TRP':{'H':1, 'N':1},
    'N+': {'N':1},
    'C-': {'O':2},
    'BBN':{'H':1, 'N':1,},
    'BBC':{'O':1},
    'NAR':{'H':1, 'N':1},
    'NAM':{'H':1, 'N':1},
    'F':  {'F':1},
    'Cl': {'Cl':1},
    'OH': {'H':1, 'O':1},
    'OP': {'O':1},
    'O3': {'O':1},
    'O2': {'O':1},
    'SH': {'S':1},
    'CG': {'H':5, 'N':3},
    'C2N':{'H':4, 'N':2},
    'OCO':{'O':2},
    'N30':{'H':4, 'N':1},
    'N31':{'H':3, 'N':1},
    'N32':{'H':2, 'N':1},
    'N33':{'H':1, 'N':1},
    'NP1':{'H':2, 'N':1},
    'N1' :{'N':1}
}

expected_atoms_base_interactions = {
    'COO':{'O':2},
    'HIS':{'N':2},
    'CYS':{'S':1},
    'TYR':{'O':1},
    'LYS':{'N':1},
    'ARG':{'N':3},
    'ROH':{'O':1},
    'AMD':{'O':1},
    'TRP':{'N':1},
    'N+': {'N':1},
    'C-': {'O':2},
    'BBN':{'H':1, 'N':1,},
    'BBC':{'O':1},
    'NAR':{'H':1, 'N':1},
    'NAM':{'H':1, 'N':1},
    'F':  {'F':1},
    'Cl': {'Cl':1},
    'OH': {'H':1, 'O':1},
    'OP': {'O':1},
    'O3': {'O':1},
    'O2': {'O':1},
    'SH': {'S':1},
    'CG': {'N':3},
    'C2N':{'N':2},
    'OCO':{'O':2},
    'N30':{'N':1},
    'N31':{'N':1},
    'N32':{'N':1},
    'N33':{'N':1},
    'NP1':{'N':1},
    'N1' :{'N':1}
}


class Group:
    def __init__(self, atom):
        #info('Made new %s group from %s'%(type,atom))
        self.atom = atom
        self.type = ''
        atom.group = self

        # set up data structures
        self.determinants = {'sidechain':[],'backbone':[],'coulomb':[]}
        self.pka_value = 0.0
        self.model_pka = 0.0

        self.Emass = 0.0
        self.Nmass = 0.0
        self.Elocl = 0.0
        self.Nlocl = 0.0
        self.buried = 0.0
        self.x = 0.0
        self.y = 0.0
        self.z = 0.0
        self.charge = 0
        self.interaction_atoms_for_acids = []
        self.interaction_atoms_for_bases = []
        self.model_pka_set = False

        # information on covalent and non-covalent coupling
        self.non_covalently_coupled_groups = []
        self.covalently_coupled_groups = []
        self.coupled_titrating_group = None
        self.common_charge_centre = False


        self.residue_type = self.atom.resName
        if self.atom.terminal:
            self.residue_type = self.atom.terminal

        if self.atom.type=='atom':
            self.label = '%-3s%4d%2s'%(self.residue_type, atom.resNumb, atom.chainID)
        elif self.atom.resName in ['DA ','DC ','DG ','DT ']:
            self.label = '%1s%1s%1s%4d%2s'%(self.residue_type[1],
                                            atom.element,
                                            atom.name.replace('\'','')[-1],
                                            atom.resNumb,
                                            atom.chainID)

#            self.label = '%1s%1s%1s%4d%2s'%(self.residue_type[1], atom.element,atom.name[-1], atom.resNumb, atom.chainID)
        else:
            self.label = '%-3s%4s%2s'%(self.residue_type, atom.name, atom.chainID)


        # container for squared distances
        self.squared_distances = {}

        return


    #
    # Coupling-related methods
    #
    def couple_covalently(self, other):
        """ Couple this group with another group """
        # do the coupling
        if not other in self.covalently_coupled_groups:
            self.covalently_coupled_groups.append(other)

        if not self in other.covalently_coupled_groups:
            other.covalently_coupled_groups.append(self)

        return

    def couple_non_covalently(self, other):
        """ Couple this group with another group """
        # do the coupling
        if not other in self.non_covalently_coupled_groups:
            self.non_covalently_coupled_groups.append(other)

        if not self in other.non_covalently_coupled_groups:
            other.non_covalently_coupled_groups.append(self)

        return

    def get_covalently_coupled_groups(self): return self.covalently_coupled_groups
    def get_non_covalently_coupled_groups(self): return self.non_covalently_coupled_groups


    def share_determinants(self, others):

        # for each determinant type
        for other in others:
            if other == self:
                continue

            for type in ['sidechain','backbone','coulomb']:
                for g in other.determinants[type]: self.share_determinant(g,type)

        # recalculate pka values
        self.calculate_total_pka()
        other.calculate_total_pka()

        return


    def share_determinant(self, new_determinant, type):
        added = False
        # first check if we already have a determinant with this label
        for own_determinant in self.determinants[type]:
            if own_determinant.group == new_determinant.group:
                # if so, find the average value
                avr = (own_determinant.value + new_determinant.value)/2.0
                own_determinant.value = avr
                new_determinant.value = avr
                added = True

        # otherwise we just add the determinant to our list
        if not added:
            self.determinants[type].append(propka.determinant.Determinant(new_determinant.group,
                                                                          new_determinant.value))

        return

    def make_covalently_coupled_line(self):

        # first check if there are any coupled groups at all
        if len(self.covalently_coupled_groups) == 0:
            return ''

        line = 'CCOUPL%5d'%self.atom.numb

        # extract and sort numbers of coupled groups
        coupled = []
        for g in self.covalently_coupled_groups:
            coupled.append(g.atom.numb)
        coupled.sort()

        # write 'em out
        for b in coupled:
            line += '%5d'%b
        line += '\n'

        return line

    def make_non_covalently_coupled_line(self):
        # first check if there are any coupled groups at all
        if len(self.non_covalently_coupled_groups) == 0:
            return ''

        line = 'NCOUPL%5d'%self.atom.numb

        # extract and sort numbers of coupled groups
        coupled = []
        for g in self.non_covalently_coupled_groups:
            coupled.append(g.atom.numb)
        coupled.sort()

        # write 'em out
        for b in coupled:
            line += '%5d'%b
        line += '\n'


        return line

    #
    # Bookkeeping methods
    #
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

    def __iadd__(self, other):
        if self.type != other.type:
            raise Exception('Cannot add groups of different types (%s and %s)'%(self.type,other.type))

        # add all values
        self.pka_value += other.pka_value
        self.Nmass += other.Nmass
        self.Emass += other.Emass
        self.Nlocl += other.Nlocl
        self.Elocl += other.Elocl
        self.buried += other.buried
        # and add all determinants
        for type in ['sidechain','backbone','coulomb']:
            for determinant in other.determinants[type]:
                self.add_determinant(determinant, type)
        return self


    def add_determinant(self, new_determinant, type):
        """ Adds to current and creates non-present determinants"""
        # first check if we already have a determinant with this label
        for own_determinant in self.determinants[type]:
            if own_determinant.group == new_determinant.group:
                # if so, add the value
                own_determinant.value += new_determinant.value
                return

        # otherwise we just add the determinant to our list
        self.determinants[type].append(propka.determinant.Determinant(new_determinant.group,
                                                                      new_determinant.value))

        return

    def set_determinant(self, new_determinant, type):
        """ Overwrites current and creates non-present determinants"""
        # first check if we already have a determinant with this label
        for own_determinant in self.determinants[type]:
            if own_determinant.group == new_determinant.group:
                # if so, overwrite the value
                own_determinant.value = new_determinant.value
                return

        # otherwise we just add the determinant to our list
        self.determinants[type].append(propka.determinant.Determinant(new_determinant.group,
                                                                      new_determinant.value))

        return

    def remove_determinants(self, labels):
        """ removes all determinants with label in labels """
        for type in ['sidechain','backbone','coulomb']:
            matches = list(filter(lambda d: d.label in labels, [d for d in self.determinants[type]]))
            for m in matches: self.determinants[type].remove(m)

        return

    def __truediv__(self, value):
        value = float(value)
        # divide all values
        self.pka_value /= value
        self.Nmass /= value
        self.Emass /= value
        self.Nlocl /= value
        self.Elocl /= value
        self.buried /= value
        # and all determinants
        for type in ['sidechain','backbone','coulomb']:
            for determinant in self.determinants[type]:
                determinant.value /= value
        return self

    def clone(self):
        res = Group(self.atom)
        res.type = self.type
        res.residue_type = self.residue_type
        res.model_pka = self.model_pka
        res.coupled_titrating_group = self.coupled_titrating_group
        res.covalently_coupled_groups = self.covalently_coupled_groups
        res.non_covalently_coupled_groups = self.non_covalently_coupled_groups
        res.titratable = self.titratable
        res.exclude_cys_from_results = self.exclude_cys_from_results
        res.charge = self.charge
        return res

    def setup(self):
        # set the charges
        if self.type in self.parameters.charge.keys():
            self.charge = self.parameters.charge[self.type]
        if self.residue_type in self.parameters.ions.keys():
            self.charge = self.parameters.ions[self.residue_type]

        #info('ION setup',self,self.residue_type, self.charge)

        # find the center and the interaction atoms
        self.setup_atoms()

        # set the model pka value
        self.titratable = False
        if self.residue_type in self.parameters.model_pkas.keys():
            if not self.model_pka_set:
                self.model_pka = self.parameters.model_pkas[self.residue_type]
                # check if we should apply a custom model pka
                key = '%s-%s'%(self.atom.resName.strip(), self.atom.name.strip())
                if key in self.parameters.custom_model_pkas.keys():
                    self.model_pka = self.parameters.custom_model_pkas[key]

                self.model_pka_set = True

        if self.model_pka_set and not self.atom.cysteine_bridge:
            self.titratable = True
        self.exclude_cys_from_results = False

        return

    def setup_atoms(self):
        # This method is overwritten for some types of groups
        # set the center at the position of the main atom
        self.set_center([self.atom])
        # set the main atom as interaction atom
        self.set_interaction_atoms([self.atom], [self.atom])
        return

    def set_interaction_atoms(self, interaction_atoms_for_acids, interaction_atoms_for_bases):
        [a.set_group_type(self.type) for a in interaction_atoms_for_acids+interaction_atoms_for_bases]

        self.interaction_atoms_for_acids = interaction_atoms_for_acids
        self.interaction_atoms_for_bases = interaction_atoms_for_bases

        # check if all atoms have been identified
        ok = True
        for [expect, found, t] in [[expected_atoms_acid_interactions, self.interaction_atoms_for_acids, 'acid'],
                                [expected_atoms_base_interactions, self.interaction_atoms_for_bases, 'base']]:
            if self.type in expect.keys():
                for e in expect[self.type].keys():
                    if len([a for a in found if a.element==e]) != expect[self.type][e]:
                        ok = False

        if not ok:
            warning('Missing atoms or failed protonation for %s (%s) -- please check the structure' % (self.label, self.type))
            warning('%s' % self)
            Na = sum([expected_atoms_acid_interactions[self.type][e] for e in expected_atoms_acid_interactions[self.type].keys()])
            Nb = sum([expected_atoms_base_interactions[self.type][e] for e in expected_atoms_base_interactions[self.type].keys()])

            warning('Expected %d interaction atoms for acids, found:' % Na)
            for i in range(len(self.interaction_atoms_for_acids)):
                 warning('             %s' % self.interaction_atoms_for_acids[i])

            warning('Expected %d interaction atoms for bases, found:' % Nb)
            for i in range(len(self.interaction_atoms_for_bases)):
                 warning('             %s' % self.interaction_atoms_for_bases[i])


                    #return

        return

    def get_interaction_atoms(self, interacting_group):
        if interacting_group.residue_type in self.parameters.base_list:
            return self.interaction_atoms_for_bases
        else:
            return self.interaction_atoms_for_acids #default is acid interaction atoms - cf. 3.0

    def set_center(self, atoms):
        if not atoms:
            raise ValueError("At least one atom must be specified")

        # reset center
        self.x = 0.0; self.y = 0.0; self.z = 0.0

        # find the average positon of atoms
        for atom in atoms:
            self.x += atom.x; self.y += atom.y; self.z += atom.z

        self.x /= float(len(atoms))
        self.y /= float(len(atoms))
        self.z /= float(len(atoms))

        return

    def getDeterminantString(self, remove_penalised_group=False):
        if self.coupled_titrating_group and remove_penalised_group:
            return ''

        number_of_sidechain = len(self.determinants['sidechain'])
        number_of_backbone  = len(self.determinants['backbone'])
        number_of_coulomb   = len(self.determinants['coulomb'])

        number_of_lines     = max(1, number_of_sidechain, number_of_backbone, number_of_coulomb)
        str  = ""
        for line_number in range(number_of_lines):
            str += "%s" % (self.label)
            if line_number == 0:
                str += " %6.2lf" %(self.pka_value)
                if len(self.non_covalently_coupled_groups)>0:
                    str+='*'
                else:
                    str+=' '

                # if self.atom.cysteine_bridge:
                #     str += " BONDED "
                # else:
                str += " %4d%2s " % (int(100.0*self.buried), "%")

                str += " %6.2lf %4d" % (self.Emass, self.Nmass)
                str += " %6.2lf %4d" % (self.Elocl, self.Nlocl)
            else:
                str += "%40s" % (" ")

            # add the determinants
            for type in ['sidechain','backbone','coulomb']:
                str += self.get_determinant_for_string(type,line_number)

            # adding end-of-line
            str += "\n"

        str += "\n"

        return str

    def get_determinant_for_string(self, type, number):
        if number >= len(self.determinants[type]):
            empty_determinant = "%s%4d%2s" % ("XXX", 0, "X")
            return "%8.2lf %s" % (0.0, empty_determinant)
        else:
            determinant = self.determinants[type][number]
            return "%8.2lf %s" % (determinant.value, determinant.label)


    def calculate_total_pka(self):
        # if this is a cysteine involved in a di-sulphide bond
        if self.atom.cysteine_bridge:
            self.pka_value = 99.99
            return


        self.pka_value = self.model_pka + self.Emass + self.Elocl

        for determinant_type in ['sidechain', 'backbone', 'coulomb']:
            for determinant in self.determinants[determinant_type]:
                self.pka_value += determinant.value

        return




    def calculate_intrinsic_pka(self):
        """
        Calculates the intrinsic pKa values from the desolvation determinants, back-bone hydrogen bonds,
        and side-chain hydrogen bond to non-titratable residues
        """
        back_bone  = 0.0
        for determinant in self.determinants['backbone']:
            value = determinant.value
            back_bone += value

        side_chain = 0.0
        for determinant in self.determinants['sidechain']:
            if determinant.label[0:3] not in ['ASP','GLU','LYS','ARG','HIS','CYS','TYR','C- ','N+ ']:
                value = determinant.value
                side_chain += value

        self.intrinsic_pKa = self.model_pka + self.Emass + self.Elocl + back_bone + side_chain

        return





    def getSummaryString(self, remove_penalised_group=False):
        if self.coupled_titrating_group and remove_penalised_group:
            return ''

        ligand_type = ''
        if self.atom.type == 'hetatm':
            ligand_type = self.type

        penalty = ''
        if self.coupled_titrating_group:
            penalty = ' NB: Discarded due to coupling with %s'%self.coupled_titrating_group.label

        str = "   %9s %8.2lf %10.2lf %18s   %s\n" % (self.label,
                                                     self.pka_value,
                                                     self.model_pka,ligand_type,
                                                     penalty)

        return str

    def __str__(self):
        return 'Group (%s) for %s'%(self.type,self.atom)



    #
    # Energy-related methods
    #

    def calculate_folding_energy(self, parameters, pH=None, reference=None):
        """
        returning the electrostatic energy of this residue at pH 'pH'
        """
        if pH == None:
            pH = parameters.pH
        if reference == None:
            reference = parameters.reference

        # If not titratable, the contribution is zero

        if not self.titratable:
            return 0.00

        # calculating the ddG(neutral --> low-pH) contribution
        ddG_neutral = 0.00
        if reference == 'neutral' and self.charge > 0.00:
            pka_prime = self.pka_value
            for determinant in self.determinants['coulomb']:
                if determinant.value > 0.00:
                    pka_prime -= determinant.value
            ddG_neutral = -1.36*(pka_prime - self.model_pka)

        # calculating the ddG(low-pH --> pH) contribution
        # folded
        x =  pH - self.pka_value
        y = 10**x
        Q_pro = math.log10(1+y)

        # unfolded
        x =  pH - self.model_pka
        y = 10**x
        Q_mod = math.log10(1+y)

        ddG_low = -1.36*(Q_pro - Q_mod)
        ddG = ddG_neutral + ddG_low

        return ddG

    def calculate_charge(self, parmaeters, pH=7.0, state='folded'):

        if state == "unfolded":
            x =  self.charge * (self.model_pka - pH)
        else:
            x =  self.charge * (self.pka_value - pH)

        y = 10**x
        charge = self.charge*(y/(1.0+y))

        return charge

    def use_in_calculations(self):
        """
        Whether this group should be included in the results report. If
        --titrate_only option is specified, only residues that are titratable
        and are in that list are included; otherwise all titratable residues
        and CYS residues are included.
        """
        return self.titratable or (self.residue_type == 'CYS' and \
                                   not self.exclude_cys_from_results)


class COO_group(Group):
    def __init__(self, atom):
        Group.__init__(self,atom)
        self.type = 'COO'

    def setup_atoms(self):
        # Identify the two caroxyl oxygen atoms
        the_oxygens = self.atom.get_bonded_elements('O')

        # set the center using the two oxygen carboxyl atoms (if present)
        if the_oxygens:
            self.set_center(the_oxygens)
        else:
            self.set_center([self.atom])
            # FIXME perhaps it would be better to ignore this group completely
            # if the oxygen is missing from this residue?

        self.set_interaction_atoms(the_oxygens, the_oxygens)
        return


class HIS_group(Group):
    def __init__(self, atom):
        Group.__init__(self,atom)
        self.type = 'HIS'

    def setup_atoms(self):
        # Find the atoms in the histidine ring
        ring_atoms = propka.ligand.is_ring_member(self.atom)
        if len(ring_atoms) != 5:
            warning('His group does not seem to contain a ring', self)

        # protonate ring
        for r in ring_atoms:
            my_protonator.protonate_atom(r)

        # set the center using the ring atoms
        if ring_atoms:
            self.set_center(ring_atoms)
        else:
            # Missing side-chain atoms
            self.set_center([self.atom])
            # FIXME perhaps it would be better to ignore this group completely?

        # find the hydrogens on the ring-nitrogens
        hydrogens = []
        nitrogens = [ra for ra in ring_atoms if ra.element == 'N']

        for nitrogen in nitrogens:
            hydrogens.extend(nitrogen.get_bonded_elements('H'))

        self.set_interaction_atoms(hydrogens+nitrogens, nitrogens)

        return


class CYS_group(Group):
    def __init__(self, atom):
        Group.__init__(self,atom)
        self.type = 'CYS'


class TYR_group(Group):
    def __init__(self, atom):
        Group.__init__(self,atom)
        self.type = 'TYR'


class LYS_group(Group):
    def __init__(self, atom):
        Group.__init__(self,atom)
        self.type = 'LYS'


class ARG_group(Group):
    def __init__(self, atom):
        Group.__init__(self,atom)
        self.type = 'ARG'

    def setup_atoms(self):
        # set the center at the position of the main atom
        self.set_center([self.atom])

        # find all the hydrogens on the nitrogen atoms
        nitrogens = self.atom.get_bonded_elements('N')
        for n in nitrogens:
            my_protonator.protonate_atom(n)

        hydrogens = []
        for nitrogen in nitrogens:
            hydrogens.extend(nitrogen.get_bonded_elements('H'))
        self.set_interaction_atoms(nitrogens+hydrogens, nitrogens)

        return


class ROH_group(Group):
    def __init__(self, atom):
        Group.__init__(self,atom)
        self.type = 'ROH'

class SER_group(Group):
    def __init__(self, atom):
        Group.__init__(self,atom)
        self.type = 'SER'

class AMD_group(Group):
    def __init__(self, atom):
        Group.__init__(self,atom)
        self.type = 'AMD'

    def setup_atoms(self):
        # Identify the oxygen and nitrogen amide atoms
        the_oxygen = self.atom.get_bonded_elements('O')
        the_nitrogen = self.atom.get_bonded_elements('N')

        # add protons to the nitrogen
        my_protonator.protonate_atom(the_nitrogen[0])
        the_hydrogens = the_nitrogen[0].get_bonded_elements('H')

        # set the center using the oxygen and nitrogen amide atoms
        self.set_center(the_oxygen+the_nitrogen)

        self.set_interaction_atoms(the_nitrogen+the_hydrogens,the_oxygen)

        return


class TRP_group(Group):
    def __init__(self, atom):
        Group.__init__(self,atom)
        self.type = 'TRP'

    def setup_atoms(self):
        # set the center at the position of the main atom
        self.set_center([self.atom])

        # find the hydrogen on the nitrogen atom
        my_protonator.protonate_atom(self.atom)
        the_hydrogen = self.atom.get_bonded_elements('H')
        self.set_interaction_atoms(the_hydrogen+[self.atom], [self.atom])
        return


class Nterm_group(Group):
    def __init__(self, atom):
        Group.__init__(self,atom)
        self.type = 'N+'


class Cterm_group(Group):
    def __init__(self, atom):
        Group.__init__(self,atom)
        self.type = 'COO'  # this is to deal with the COO-C- parameter unification.

    def setup_atoms(self):
        # Identify the carbon and other oxygen carboxyl atoms
        the_carbons = self.atom.get_bonded_elements('C')
        if not the_carbons:
            self.set_center([self.atom])
            # FIXME perhaps it would be better to ignore this group completely
            # if the carbon is missing from this residue?
        else:
            the_other_oxygen = the_carbons[0].get_bonded_elements('O')
            the_other_oxygen.remove(self.atom)

            # set the center and interaction atoms
            the_oxygens = [self.atom]+ the_other_oxygen
            self.set_center(the_oxygens)
            self.set_interaction_atoms(the_oxygens, the_oxygens)


        return


class BBN_group(Group):
    def __init__(self, atom):
        Group.__init__(self,atom)
        self.type = 'BBN'
        self.residue_type = 'BBN'


    def setup_atoms(self):
        # Identify the hydrogen
        my_protonator.protonate_atom(self.atom)
        the_hydrogen = self.atom.get_bonded_elements('H')

        # set the center using the nitrogen
        self.set_center([self.atom])
        self.set_interaction_atoms(the_hydrogen+[self.atom], the_hydrogen+[self.atom])
        return


class BBC_group(Group):
    def __init__(self, atom):
        Group.__init__(self,atom)
        self.type = 'BBC'
        self.residue_type = 'BBC'

    def setup_atoms(self):
        # Identify the oxygen
        the_oxygen = self.atom.get_bonded_elements('O')

        # set the center using the nitrogen
        self.set_center([self.atom])
        self.set_interaction_atoms(the_oxygen, the_oxygen)
        return


class NAR_group(Group):
    def __init__(self, atom):
        Group.__init__(self,atom)
        self.type = 'NAR'
        self.residue_type = 'NAR'
        info('Found NAR group:', atom)
        return


    def setup_atoms(self):
        # Identify the hydrogen
        my_protonator.protonate_atom(self.atom)
        the_hydrogen = self.atom.get_bonded_elements('H')

        # set the center using the nitrogen
        self.set_center([self.atom])
        self.set_interaction_atoms(the_hydrogen+[self.atom], the_hydrogen+[self.atom])
        return





class NAM_group(Group):
    def __init__(self, atom):
        Group.__init__(self,atom)
        self.type = 'NAM'
        self.residue_type = 'NAM'
        info('Found NAM group:', atom)
        return


    def setup_atoms(self):
        # Identify the hydrogen
        my_protonator.protonate_atom(self.atom)
        the_hydrogen = self.atom.get_bonded_elements('H')

        # set the center using the nitrogen
        self.set_center([self.atom])
        self.set_interaction_atoms(the_hydrogen+[self.atom], the_hydrogen+[self.atom])
        return



class F_group(Group):
    def __init__(self, atom):
        Group.__init__(self,atom)
        self.type = 'F'
        self.residue_type = 'F'
        info('Found F   group:', atom)
        return

class Cl_group(Group):
    def __init__(self, atom):
        Group.__init__(self,atom)
        self.type = 'Cl'
        self.residue_type = 'Cl'
        info('Found Cl   group:', atom)
        return

class OH_group(Group):
    def __init__(self, atom):
        Group.__init__(self,atom)
        self.type = 'OH'
        self.residue_type = 'OH'
        info('Found OH  group:', atom)
        return


    def setup_atoms(self):
        # Identify the hydrogen
        my_protonator.protonate_atom(self.atom)
        the_hydrogen = self.atom.get_bonded_elements('H')

        # set the center using the nitrogen
        self.set_center([self.atom])
        self.set_interaction_atoms(the_hydrogen+[self.atom], the_hydrogen+[self.atom])
        return

class OP_group(Group):
    def __init__(self, atom):
        Group.__init__(self,atom)
        self.type = 'OP'
        self.residue_type = 'OP'
        info('Found OP  group:', atom)
        return


    def setup_atoms(self):
        # Identify the hydrogen
        my_protonator.protonate_atom(self.atom)
        #the_hydrogen = self.atom.get_bonded_elements('H')

        # set the center using the oxygen
        self.set_center([self.atom])
        #self.set_interaction_atoms(the_hydrogen+[self.atom], the_hydrogen+[self.atom])
        self.set_interaction_atoms([self.atom], [self.atom])
        return


class O3_group(Group):
    def __init__(self, atom):
        Group.__init__(self,atom)
        self.type = 'O3'
        self.residue_type = 'O3'
        info('Found O3  group:', atom)
        return


class O2_group(Group):
    def __init__(self, atom):
        Group.__init__(self,atom)
        self.type = 'O2'
        self.residue_type = 'O2'
        info('Found O2  group:', atom)
        return

class SH_group(Group):
    def __init__(self, atom):
        Group.__init__(self,atom)
        self.type = 'SH'
        self.residue_type = 'SH'
        info('Found SH  group:', atom)
        return


class CG_group(Group):
    """Guadinium group"""
    def __init__(self, atom):
        Group.__init__(self,atom)
        self.type = 'CG'
        self.residue_type = 'CG'
        info('Found CG  group:', atom)
        return

    def setup_atoms(self):
        # Identify the nitrogens
        the_nitrogens = self.atom.get_bonded_elements('N')

        # set the center using the nitrogen
        self.set_center([self.atom])

        the_hydrogens = []
        for n in the_nitrogens:
            my_protonator.protonate_atom(n)
            the_hydrogens += n.get_bonded_elements('H')
        self.set_interaction_atoms(the_hydrogens+the_nitrogens, the_nitrogens)

        return

class C2N_group(Group):
    """Amidinium group"""
    def __init__(self, atom):
        Group.__init__(self,atom)
        self.type = 'C2N'
        self.residue_type = 'C2N'
        info('Found C2N group:', atom)
        return

    def setup_atoms(self):
        # Identify the nitrogens
        the_nitrogens = self.atom.get_bonded_elements('N')
        the_nitrogens = [n for n in the_nitrogens if len(n.get_bonded_heavy_atoms())==1]

        # set the center using the nitrogen
        self.set_center([self.atom])

        the_hydrogens = []
        for n in the_nitrogens:
            my_protonator.protonate_atom(n)
            the_hydrogens += n.get_bonded_elements('H')

        self.set_interaction_atoms(the_hydrogens+the_nitrogens, the_nitrogens)
        return

class OCO_group(Group):
    def __init__(self, atom):
        Group.__init__(self,atom)
        self.type = 'OCO'
        self.residue_type = 'OCO'
        info('Found OCO group:', atom)
        return

    def setup_atoms(self):
        # Identify the two caroxyl oxygen atoms
        the_oxygens = self.atom.get_bonded_elements('O')

        # set the center using the two oxygen carboxyl atoms
        self.set_center(the_oxygens)
        self.set_interaction_atoms(the_oxygens, the_oxygens)
        return



class N30_group(Group):
    def __init__(self, atom):
        Group.__init__(self,atom)
        self.type = 'N30'
        self.residue_type = 'N30'
        info('Found N30 group:', atom)
        return

    def setup_atoms(self):
        # Identify the nitrogens
        my_protonator.protonate_atom(self.atom)
        the_hydrogens = self.atom.get_bonded_elements('H')
        # set the center using the nitrogen
        self.set_center([self.atom])
        self.set_interaction_atoms(the_hydrogens+[self.atom], [self.atom])
        return

class N31_group(Group):
    def __init__(self, atom):
        Group.__init__(self,atom)
        self.type = 'N31'
        self.residue_type = 'N31'
        info('Found N31 group:', atom)
        return

    def setup_atoms(self):
        # Identify the nitrogens
        my_protonator.protonate_atom(self.atom)
        the_hydrogens = self.atom.get_bonded_elements('H')
        # set the center using the nitrogen
        self.set_center([self.atom])
        self.set_interaction_atoms(the_hydrogens+[self.atom], [self.atom])
        return

class N32_group(Group):
    def __init__(self, atom):
        Group.__init__(self,atom)
        self.type = 'N32'
        self.residue_type = 'N32'
        info('Found N32 group:', atom)
        return

    def setup_atoms(self):
        # Identify the nitrogens
        my_protonator.protonate_atom(self.atom)
        the_hydrogens = self.atom.get_bonded_elements('H')
        # set the center using the nitrogen
        self.set_center([self.atom])
        self.set_interaction_atoms(the_hydrogens+[self.atom], [self.atom])
        return

class N33_group(Group):
    def __init__(self, atom):
        Group.__init__(self,atom)
        self.type = 'N33'
        self.residue_type = 'N33'
        info('Found N33 group:', atom)
        return

    def setup_atoms(self):
        # Identify the nitrogens
        my_protonator.protonate_atom(self.atom)
        the_hydrogens = self.atom.get_bonded_elements('H')
        # set the center using the nitrogen
        self.set_center([self.atom])
        self.set_interaction_atoms(the_hydrogens+[self.atom], [self.atom])
        return

class NP1_group(Group):
    def __init__(self, atom):
        Group.__init__(self,atom)
        self.type = 'NP1'
        self.residue_type = 'NP1'
        info('Found NP1 group:', atom)
        return


    def setup_atoms(self):
        # Identify the nitrogens
        my_protonator.protonate_atom(self.atom)
        the_hydrogens = self.atom.get_bonded_elements('H')
        # set the center using the nitrogen
        self.set_center([self.atom])
        self.set_interaction_atoms(the_hydrogens+[self.atom], [self.atom])
        return

class N1_group(Group):
    def __init__(self, atom):
        Group.__init__(self,atom)
        self.type = 'N1'
        self.residue_type = 'N1'
        info('Found N1 group:', atom)
        return



class Ion_group(Group):
    def __init__(self, atom):
        Group.__init__(self,atom)
        self.type = 'ION'
        self.residue_type = atom.resName.strip()
        info('Found ion group:', atom)
        return


class non_titratable_ligand_group(Group):
    def __init__(self, atom):
        Group.__init__(self, atom)
        self.type = 'LG'
        self.residue_type = 'LG'
#        info('Non-titratable ligand group',atom)
        return

class titratable_ligand_group(Group):
    def __init__(self, atom):
        Group.__init__(self, atom)
        # set the charge and determine type (acid or base)
        self.charge = atom.charge
        if self.charge <0:
            self.type = 'ALG'
            self.residue_type = 'ALG'
        elif self.charge > 0:
            self.type = 'BLG'
            self.residue_type = 'BLG'
        else:
            raise Exception('Unable to determine type of ligand group - charge not set?')


        # check if marvin model pka has been calculated
        # this is not true if we are reading an input file
        if atom.marvin_pka:
            self.model_pka = atom.marvin_pka
            info('Titratable ligand group    ', atom, self.model_pka, self.charge)
        self.model_pka_set = True

        return


def is_group(parameters, atom):
    atom.groups_extracted = 1

    # check if this atom belongs to a protein group
    protein_group = is_protein_group(parameters, atom)
    if protein_group: return protein_group

    # check if this atom belongs to a ion group
    ion_group = is_ion_group(parameters, atom)
    if ion_group: return ion_group

    # check if this atom belongs to a ligand group
    if parameters.ligand_typing == 'marvin':
        ligand_group = is_ligand_group_by_marvin_pkas(parameters, atom)
    elif parameters.ligand_typing == 'sybyl':
        ligand_group = is_ligand_group_by_sybyl_types(parameters, atom)
    elif parameters.ligand_typing == 'groups':
        ligand_group = is_ligand_group_by_groups(parameters, atom)
    else:
        raise Exception('Unknown ligand typing method \'%s\''%parameters.ligand_typing)

    if ligand_group: return ligand_group



    return None


def is_protein_group(parameters,atom):
    if atom.type != 'atom':
        return None

    ### Check for termial groups
    if atom.terminal == 'N+':
        return Nterm_group(atom)
    elif atom.terminal == 'C-':
        return Cterm_group(atom)

    ### Backbone
    if atom.type == 'atom' and atom.name == 'N':
        # ignore proline backbone nitrogens
        if atom.resName != 'PRO':
            return BBN_group(atom)
    if atom.type == 'atom' and atom.name == 'C':
        # ignore C- carboxyl
        if atom.count_bonded_elements('O') == 1:
            return BBC_group(atom)

    ### Filters for side chains based on PDB protein atom names
    key = '%s-%s'%(atom.resName, atom.name)

    if key in parameters.protein_group_mapping.keys():
        return eval('%s_group(atom)'%parameters.protein_group_mapping[key])

    return None

def is_ligand_group_by_sybyl_types(parameters, atom):


    return None


def is_ligand_group_by_groups(parameters, atom):
    ### Ligand group filters
    if atom.type != 'hetatm':
        return None

    my_protonator.protonate_atom(atom)

    if atom.sybyl_type == 'N.ar':
        if len(atom.get_bonded_heavy_atoms())==2:
            return NAR_group(atom)

    if atom.sybyl_type == 'N.am':
        return NAM_group(atom)

    if atom.sybyl_type in ['N.3', 'N.4']:
        heavy_bonded = atom.get_bonded_heavy_atoms()
        if len(heavy_bonded) == 0:
            return N30_group(atom)
        elif len(heavy_bonded) == 1:
            return N31_group(atom)
        elif len(heavy_bonded) == 2:
            return N32_group(atom)
        elif len(heavy_bonded) == 3:
            return N33_group(atom)

    if atom.sybyl_type == 'N.1':
        return N1_group(atom)

    if atom.sybyl_type == 'N.pl3':
        # make sure that this atom is not part of a guadinium or amidinium group
        bonded_carbons = atom.get_bonded_elements('C')
        if len(bonded_carbons) == 1:
            bonded_nitrogens = bonded_carbons[0].get_bonded_elements('N')
            if len(bonded_nitrogens) == 1:
                return NP1_group(atom)


    if atom.sybyl_type == 'C.2':
        # Guadinium and amidinium groups
        bonded_nitrogens = atom.get_bonded_elements('N')
        npls = [n for n in bonded_nitrogens if (n.sybyl_type == 'N.pl3' and len(n.get_bonded_heavy_atoms())==1)]
        if len(npls) == 2:
            n_with_max_two_heavy_atom_bonds = [n for n in bonded_nitrogens if len(n.get_bonded_heavy_atoms())<3]
            if len(n_with_max_two_heavy_atom_bonds) == 2:
                return C2N_group(atom)
            if len(bonded_nitrogens) == 3:
                return CG_group(atom)
        # carboxyl group
        bonded_oxygens = atom.get_bonded_elements('O')
        bonded_oxygens = [b for b in bonded_oxygens if 'O.co2' in b.sybyl_type]
        if len(bonded_oxygens) == 2:
            return OCO_group(atom)


    if atom.sybyl_type == 'F':
        return F_group(atom)

    if atom.sybyl_type == 'Cl':
        return Cl_group(atom)

    if atom.sybyl_type == 'O.3':
        if len(atom.get_bonded_heavy_atoms()) == 1:
            # phosphate group
            if atom.count_bonded_elements('P') == 1:
                return OP_group(atom)
            # hydroxyl group
            else:
                return OH_group(atom)
        # other SP3 oxygen
        else:
            return O3_group(atom)

    if atom.sybyl_type == 'O.2':
        return O2_group(atom)


    if atom.sybyl_type == 'S.3':
        # sulphide group
        if len(atom.get_bonded_heavy_atoms()) == 1:
            return SH_group(atom)
            # other SP3 oxygen
        #else:
        #    return S3_group(atom)


    return None


def is_ligand_group_by_marvin_pkas(parameters, atom):
    if atom.type != 'hetatm':
        return None

    # calculate Marvin ligand pkas for this conformation container
    # if not already done
    if not atom.conformation_container.marvin_pkas_calculated:
        lpka = propka.ligand_pka_values.ligand_pka_values(parameters)
        lpka.get_marvin_pkas_for_molecular_container(atom.molecular_container,
                                                     min_pH=parameters.min_ligand_model_pka,
                                                     max_pH=parameters.max_ligand_model_pka)


    if atom.marvin_pka:
        return titratable_ligand_group(atom)

    # Special case of oxygen in carboxyl group not assigned a pka value by marvin
    if atom.sybyl_type == 'O.co2':
        atom.charge = -1.0
        other_oxygen = [o for o in atom.get_bonded_elements('C')[0].get_bonded_elements('O') if not o==atom][0]
        atom.marvin_pka = other_oxygen.marvin_pka
        return titratable_ligand_group(atom)


    if atom.element in parameters.hydrogen_bonds.elements:
        return non_titratable_ligand_group(atom)

    return None


def is_ion_group(parameters, atom):

    if atom.resName.strip() in parameters.ions.keys():
        return Ion_group(atom)

    return None
