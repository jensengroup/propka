
from __future__ import division
from __future__ import print_function

import math
import propka.lib as lib
import sys, os
from propka.lib import info, warning

import pkg_resources

# names and types of all key words in configuration file
matrices =            ['interaction_matrix']

pair_wise_matrices =  ['sidechain_cutoffs']

number_dictionaries = ['VanDerWaalsVolume','charge','model_pkas','ions',
                      'valence_electrons','custom_model_pkas']

list_dictionaries =   ['backbone_NH_hydrogen_bond','backbone_CO_hydrogen_bond']

string_dictionaries = ['protein_group_mapping']

string_lists =        ['ignore_residues','angular_dependent_sidechain_interactions',
                       'acid_list','base_list','exclude_sidechain_interactions',
                       'backbone_reorganisation_list','write_out_order']

distances =           ['desolv_cutoff','buried_cutoff','coulomb_cutoff1','coulomb_cutoff2']

parameters =          ['Nmin','Nmax','desolvationSurfaceScalingFactor','desolvationPrefactor',
                       'desolvationAllowance','coulomb_diel','COO_HIS_exception','OCO_HIS_exception',
                       'CYS_HIS_exception','CYS_CYS_exception','min_ligand_model_pka','max_ligand_model_pka',
                       'include_H_in_interactions','coupling_max_number_of_bonds',
                       'min_bond_distance_for_hydrogen_bonds','coupling_penalty',
                       'shared_determinants','common_charge_centre','hide_penalised_group', 'remove_penalised_group',
                       'max_intrinsic_pKa_diff','min_interaction_energy','max_free_energy_diff','min_swap_pka_shift',
                       'min_pka','max_pka','sidechain_interaction']

strings =             ['version','output_file_tag','ligand_typing','pH','reference']




class Parameters:
    def __init__(self, parameter_file):

        self.set_up_data_structures()
        self.read_parameters(parameter_file)

        #self.print_interaction_parameters()
        #self.print_interaction_parameters_latex()
        #####self.print_interactions_latex()
        #sys.exit(0)


        return


    def read_parameters(self, file):
        # try to locate the parameters file
        try:
            ifile = pkg_resources.resource_filename(__name__, file)
            input = lib.open_file_for_reading(ifile)
        except:
            input = lib.open_file_for_reading(file)

        for line in input:
            self.parse_line(line)

        return


    def parse_line(self, line):
        # first, remove comments
        comment_pos = line.find('#')
        if comment_pos != -1:
            line = line[:comment_pos]

        # split the line into words
        words = line.split()
        if len(words) == 0:
            return

        # parse the words
        if len(words)==3 and words[0] in number_dictionaries:
            self.parse_to_number_dictionary(words)
        elif len(words)==2 and words[0] in string_lists:
            self.parse_to_string_list(words)
        elif len(words)==2 and words[0] in distances:
            self.parse_distance(words)
        elif len(words)==2 and words[0] in parameters:
            self.parse_parameter(words)
        elif len(words)==2 and words[0] in strings:
            self.parse_string(words)
        elif len(words)>2 and words[0] in list_dictionaries:
            self.parse_to_list_dictionary(words)
        elif words[0] in matrices+pair_wise_matrices:
            self.parse_to_matrix(words)
        elif len(words)==3 and words[0] in string_dictionaries:
            self.parse_to_string_dictionary(words)


        #info(words)

        return


    def parse_to_number_dictionary(self, words):
        exec('self.%s[\'%s\'] = %s'%tuple(words))
        return

    def parse_to_string_dictionary(self, words):
        exec('self.%s[\'%s\'] = \'%s\''%tuple(words))
        return

    def parse_to_list_dictionary(self, words):
        exec('if not \'%s\' in self.%s.keys(): self.%s[\'%s\'] = []'%(words[1],words[0],words[0],words[1]))
        for word in words[2:]:
            exec('self.%s[\'%s\'].append(%s)'%(words[0],words[1],word))

        return

    def parse_to_string_list(self, words):
        exec('self.%s.append(\'%s\')'%tuple(words))
        return

    def parse_to_matrix(self, words):
        exec('self.%s.add(%s)'%(words[0],tuple(words[1:])))
        return

    def parse_distance(self, words):
        # float check needed
        exec('self.%s = %s'%tuple(words))
        exec('self.%s_squared = pow(%s,2)'%tuple(words))
        return

    def parse_parameter(self, words):
        exec('self.%s = %s'%tuple(words))
        return

    def parse_string(self, words):
        #info('self.%s = \'%s\''%tuple(words))
        exec('self.%s = \'%s\''%tuple(words))
        return


    def set_up_data_structures(self):
        for key_word in number_dictionaries+list_dictionaries+string_dictionaries:
            exec('self.%s = {}'%key_word)
        for key_word in string_lists:
            exec('self.%s = []'%key_word)
        for key_word in strings:
            exec('self.%s = \'\''%key_word)
        for key_word in matrices:
            exec('self.%s = Interaction_matrix(\'%s\')'%(key_word,key_word))
        for key_word in pair_wise_matrices:
            exec('self.%s =Pair_wise_matrix(\'%s\')'%(key_word,key_word))

        return

    def print_interaction_parameters(self):
        info('--------------- Model pKa values ----------------------')
        for k in self.model_pkas.keys():
            info('%3s %8.2f' % (k, self.model_pkas[k]))

        info('')
        info('--------------- Interactions --------------------------')
        agroups = ['COO', 'HIS', 'CYS', 'TYR', 'SER', 'N+', 'LYS', 'AMD', 'ARG', 'TRP', 'ROH', 'CG', 'C2N', 'N30', 'N31', 'N32', 'N33', 'NAR', 'OCO', 'NP1', 'OH', 'O3', 'CL', 'F', 'NAM', 'N1', 'O2', 'OP', 'SH']
        lgroups = ['CG', 'C2N', 'N30', 'N31', 'N32', 'N33', 'NAR', 'OCO', 'NP1', 'OH', 'O3', 'CL', 'F', 'NAM', 'N1', 'O2', 'OP', 'SH']

        map = {
            'CG' :['ARG'],
            'C2N':['ARG'],
            'N30':['N+','LYS'],
            'N31':['N+','LYS'],
            'N32':['N+','LYS'],
            'N33':['N+','LYS'] ,
            'NAR':['HIS'],
            'OCO':['COO'],
            'OP' :[],#['TYR','SER'],
            'SH' :['CYS'] ,
            'NP1':[],
            'OH' :['ROH'],
            'O3' :[] ,
            'CL' :[],
            'F'  :[],
            'NAM':['AMD'],
            'N1' :[],
            'O2' :[]}


        for g1 in agroups:
            for g2 in lgroups:

                interaction = '%3s %3s %1s %4s %4s'%(g1,g2,
                                                    self.interaction_matrix[g1][g2],
                                                    self.sidechain_cutoffs.get_value(g1,g2)[0],
                                                    self.sidechain_cutoffs.get_value(g1,g2)[1])

                map_interaction = ''
                if g2 in map:
                    for m in map[g2]:
                        map_interaction += '|%3s %3s %1s %4s %4s'%(g1,m,
                                                                  self.interaction_matrix[g1][m],
                                                                  self.sidechain_cutoffs.get_value(g1,m)[0],
                                                                  self.sidechain_cutoffs.get_value(g1,m)[1])
                        if self.interaction_matrix[g1][m] != self.interaction_matrix[g1][g2]:
                            map_interaction += '* '
                        if self.sidechain_cutoffs.get_value(g1,m)[0] != self.sidechain_cutoffs.get_value(g1,g2)[0] or \
                                self.sidechain_cutoffs.get_value(g1,m)[1] != self.sidechain_cutoffs.get_value(g1,g2)[1]:
                            map_interaction += '! '
                        else:
                            map_interaction += '  '
                    if len(map[g2])==0 and (self.sidechain_cutoffs.get_value(g1,g2)[0] !=3 or self.sidechain_cutoffs.get_value(g1,g2)[1] != 4):
                        map_interaction += '?  '

                info(interaction, map_interaction)

                if g1==g2:
                    break
            info('-')

        info('--------------- Exceptions ----------------------------')
        info('COO-HIS', self.COO_HIS_exception)
        info('OCO-HIS', self.OCO_HIS_exception)
        info('CYS-HIS', self.CYS_HIS_exception)
        info('CYS-CYS', self.CYS_CYS_exception)


        info('--------------- Mapping -------------------------------')
        info("""
Titratable:
CG  ARG
C2N ARG
N30 N+/LYS
N31 N+/LYS
N32 N+/LYS
N33 N+/LYS
NAR HIS
OCO COO
OP  TYR/SER?
SH  CYS

Non-titratable:
NP1 AMD?
OH  ROH
O3  ?
CL
F
NAM
N1
O2
""")
        return






    def print_interaction_parameters_latex(self):
#         info('--------------- Model pKa values ----------------------')
#         for k in self.model_pkas.keys():
#             info('%3s %8.2f'%(k,self.model_pkas[k]))

#         info('')
#         info('--------------- Interactions --------------------------')
        agroups = ['COO', 'HIS', 'CYS', 'TYR', 'SER', 'N+', 'LYS', 'AMD', 'ARG', 'TRP', 'ROH', 'CG', 'C2N', 'N30', 'N31', 'N32', 'N33', 'NAR', 'OCO', 'NP1', 'OH', 'O3', 'CL', 'F', 'NAM', 'N1', 'O2', 'OP', 'SH']
        lgroups = ['CG', 'C2N', 'N30', 'N31', 'N32', 'N33', 'NAR', 'OCO', 'NP1', 'OH', 'O3', 'CL', 'F', 'NAM', 'N1', 'O2', 'OP', 'SH']

        map = {
            'CG' :['ARG'],
            'C2N':['ARG'],
            'N30':['N+','LYS'],
            'N31':['N+','LYS'],
            'N32':['N+','LYS'],
            'N33':['N+','LYS'] ,
            'NAR':['HIS'],
            'OCO':['COO'],
            'OP' :[],#['TYR','SER'],
            'SH' :['CYS'] ,
            'NP1':['AMD'],
            'OH' :['ROH'],
            'O3' :[] ,
            'CL' :[],
            'F'  :[],
            'NAM':[],
            'N1' :[],
            'O2' :[]}


        s = """
\\begin{longtable}{lllll}
\\caption{Ligand interaction parameters. For interactions not listed, the default value of %s is applied.}
\\label{tab:ligand_interaction_parameters}\\\\

\\toprule
Group1 & Group2 & Interaction & c1 &c2 \\\\
\\midrule
\\endfirsthead

\\multicolumn{5}{l}{\\emph{continued from the previous page}}\\\\
\\toprule
Group1 & Group2 & Interaction & c1 &c2 \\\\
\\midrule
\\endhead

\\midrule
\\multicolumn{5}{r}{\\emph{continued on the next page}}\\\\
\\endfoot

\\bottomrule
\\endlastfoot

"""%(self.sidechain_cutoffs.default)
        for g1 in agroups:
            for g2 in lgroups:
                if self.interaction_matrix[g1][g2]=='-':
                    continue
                if self.sidechain_cutoffs.get_value(g1,g2)==self.sidechain_cutoffs.default:
                    continue


                s+= '%3s & %3s & %1s & %4s & %4s\\\\ \n'%(g1,g2,
                                                       self.interaction_matrix[g1][g2],
                                                       self.sidechain_cutoffs.get_value(g1,g2)[0],
                                                       self.sidechain_cutoffs.get_value(g1,g2)[1])

                if g1==g2:
                    break

        s += '  \\end{longtable}\n'
        info(s)
        return

    def print_interactions_latex(self):
        agroups = ['COO', 'HIS', 'CYS', 'TYR', 'SER', 'N+', 'LYS', 'AMD', 'ARG', 'TRP', 'ROH', 'CG', 'C2N', 'N30', 'N31', 'N32', 'N33', 'NAR', 'OCO', 'NP1', 'OH', 'O3', 'CL', 'F', 'NAM', 'N1', 'O2', 'OP', 'SH']
        lgroups = ['CG', 'C2N', 'N30', 'N31', 'N32', 'N33', 'NAR', 'OCO', 'NP1', 'OH', 'O3', 'CL', 'F', 'NAM', 'N1', 'O2', 'OP', 'SH']


        s = """
\\begin{longtable}{%s}
\\caption{Ligand interaction parameters. For interactions not listed, the default value of %s is applied.}
\\label{tab:ligand_interaction_parameters}\\\\

\\toprule
Group1 & Group2 & Interaction & c1 &c2 \\\\
\\midrule
\\endfirsthead

\\multicolumn{5}{l}{\\emph{continued from the previous page}}\\\\
\\toprule
Group1 & Group2 & Interaction & c1 &c2 \\\\
\\midrule
\\endhead

\\midrule
\\multicolumn{5}{r}{\\emph{continued on the next page}}\\\\
\\endfoot

\\bottomrule
\\endlastfoot

"""%('l'*len(agroups),self.sidechain_cutoffs.default)
        for g1 in agroups:
            for g2 in agroups:

                s+= '%3s & %3s & %1s & %4s & %4s\\\\ \n'%(g1,g2,
                                                       self.interaction_matrix[g1][g2],
                                                       self.sidechain_cutoffs.get_value(g1,g2)[0],
                                                       self.sidechain_cutoffs.get_value(g1,g2)[1])

                if g1==g2:
                    break

        s += '  \\end{longtable}\n'
        info(s)
        return




class Interaction_matrix:
    def __init__(self, name):
        self.name = name
        self.ordered_keys = []
        self.dictionary = {}
        return

    def add(self,words):
        new_group = words[0]
        self.ordered_keys.append(new_group)

        if not new_group in self.dictionary.keys():
            self.dictionary[new_group] = {}

        for i in range(len(self.ordered_keys)):
            group = self.ordered_keys[i]
            if len(words)>i+1:
                try:
                    exec('self.value = %s'%words[i+1])
                except:
                    self.value = words[i+1]
                self.dictionary[group][new_group] = self.value
                self.dictionary[new_group][group] = self.value


        return

    def get_value(self, item1, item2):
        try:
            return self.dictionary[item1][item2]
        except:
            return None

    def __getitem__(self, group):
        if group not in self.dictionary.keys():
            raise Exception('%s not found in interaction matrix %s'%(group,self.name))
        return self.dictionary[group]


    def keys(self):
        return self.dictionary.keys()

    def __str__(self):
        s = '      '
        for k1 in self.ordered_keys:
            s+='%3s '%k1
        s+='\n'
        for k1 in self.ordered_keys:
            s+='%3s '%k1
            for k2 in self.ordered_keys:
                s+='%3s '%self[k1][k2]
            s+='\n'

        return s
#         ks = ['COO', 'SER', 'ARG', 'LYS', 'HIS', 'AMD', 'CYS', 'TRP','ROH','TYR','N+','CG', 'C2N', 'N30', 'N31', 'N32', 'N33', 'NAR', 'OCO', 'NP1', 'OH', 'O3', 'CL', 'F', 'NAM', 'N1', 'O2', 'OP', 'SH']

#         p = ''
#         n=0
#         for i in range(len(ks)):
#             for j in range(i,len(ks)):
#                 if not [0.0,0.0]==self[ks[i]][ks[j]]:
#                     if not [3.0,4.0]==self[ks[i]][ks[j]]:
#                         p+='sidechain_cutoff %3s %3s %s\n'%(ks[i],ks[j],self[ks[i]][ks[j]])
#                         n+=1

#         info('total',n,len(ks))
#         return p



class Pair_wise_matrix:
    def __init__(self, name):
        self.name = name
        self.dictionary = {}
        self.default = [0.0, 0.0]
        return

    def add(self,words):
        # assign the default value
        if len(words)==3 and words[0]=='default':
            self.default = [float(words[1]), float(words[2])]
            return

        # assign non-default values
        g1 = words[0]
        g2 = words[1]
        v = [float(words[2]), float(words[3])]

        self.insert(g1,g2,v)
        self.insert(g2,g1,v)

        return

    def insert(self, k1,k2,v):

        if k1 in self.dictionary.keys() and k2 in self.dictionary[k1].keys():
            if k1!=k2:
                warning('Parameter value for %s, %s defined more than once' % (k1, k2))

        if not k1 in self.dictionary:
            self.dictionary[k1] = {}

        self.dictionary[k1][k2] =v

        return

    def get_value(self, item1, item2):

        try:
            return self.dictionary[item1][item2]
        except:
            return self.default

    def __getitem__(self, group):
        if group not in self.dictionary.keys():
            raise Exception('%s not found in interaction matrix %s'%(group,self.name))
        return self.dictionary[group]


    def keys(self):
        return self.dictionary.keys()

    def __str__(self):
        s=''
        for k1 in self.keys():
            for k2 in self[k1].keys():
                s += '%s %s %s\n'%(k1,k2,self[k1][k2])

        return s























