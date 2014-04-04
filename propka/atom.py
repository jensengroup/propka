
from __future__ import division
from __future__ import print_function

import string, propka.lib, propka.group

from . import hybrid36

class Atom:
    """
      Atom class - contains all atom information found in the pdbfile
    """

    def __init__(self, line=None, verbose=False):

        self.set_properties(line)

        self.residue_label = "%-3s%4d%2s" % (self.name,self.resNumb, self.chainID)

        self.groups_extracted = 0
        self.group = None
        self.group_type = None
        self.number_of_bonded_elements = {}

        self.cysteine_bridge = False

        self.bonded_atoms = []
        self.residue = None
        self.conformation_container = None
        self.molecular_container = None
        self.is_protonated = False
        self.steric_number_and_lone_pairs_set = False

        self.terminal = None
        self.charge = 0
        self.charge_set = False
        self.steric_number = 0
        self.number_of_lone_pairs = 0
        self.number_of_protons_to_add = 0
        self.number_of_pi_electrons_in_double_and_triple_bonds = 0
        self.number_of_pi_electrons_in_conjugate_double_and_triple_bonds = 0

        # ligand atom types
        self.sybyl_type = ''
        self.sybyl_assigned = False
        self.marvin_pka = False


        return




    def set_properties(self, line):

        self.name = ''
        self.numb = 0
        self.x = 0.0
        self.y = 0.0
        self.z = 0.0
        self.resNumb = 0
        self.resName = ''
        self.chainID = 'A'
        self.type = ''
        self.occ = '1.0'
        self.beta = '0.0'
        self.element = ''
        self.icode = ''

        if line:
            self.name = line[12:16].strip()
            self.numb = int( hybrid36.decode(line[ 6:11]) )
            self.x = float( line[30:38].strip() )
            self.y = float( line[38:46].strip() )
            self.z = float( line[46:54].strip() )
            self.resNumb = int( line[22:26].strip() )
            self.resName = "%-3s" % (line[17:20].strip())
            self.chainID = line[21]
            # Set chain id to "_" if it is just white space.
            if not self.chainID.strip():
                self.chainID = '_'
            self.type = line[:6].strip().lower()

            if self.resName in ['DA ','DC ','DG ','DT ']:
                self.type = 'hetatm'

            self.occ = line[55:60].strip()
            self.beta = line[60:66].strip()
            self.icode = line[26:27]

            # Set the element using the position of the name in the pdb file
            self.element = line[12:14].strip().strip(string.digits)
            if len(self.name) == 4:
                self.element = self.element[0]
            if len(self.element)==2:
                self.element = '%1s%1s'%(self.element[0], self.element[1].lower())

        return

    def set_group_type(self, type):
        self.group_type = type
        return

    def count_bonded_elements(self, element):
        res = 0
        for ba in self.bonded_atoms:
            if element == ba.element:
                res +=1
        return res


    def get_bonded_elements(self, element):
        res = []
        for ba in self.bonded_atoms:
                if ba.element == element:
                    res.append(ba)
        return res

    def get_bonded_heavy_atoms(self):
        return [ba for ba in self.bonded_atoms if ba.element!='H']


    def is_atom_within_bond_distance(self, other_atom, max_bonds, cur_bond):
        """ check if <other_atom> is found within <max_bonds> bonds of self """
        for ba in self.bonded_atoms:
            if ba == other_atom:
                return True
            if max_bonds > cur_bond:
                if ba.is_atom_within_bond_distance(other_atom, max_bonds, cur_bond+1):
                    return True

        return False



    def setProperty(self,
                    numb    = None,
                    name    = None,
                    resName = None,
                    chainID = None,
                    resNumb = None,
                    x       = None,
                    y       = None,
                    z       = None,
                    occ     = None,
                    beta    = None):
        """
        sets properties of the atom object
        """

        if numb    != None: self.numb    = numb
        if name    != None: self.name    = name
        if resName != None: self.resName = resName
        if chainID != None: self.chainID = chainID
        if resNumb != None: self.resNumb = resNumb
        if x       != None: self.x       = x
        if y       != None: self.y       = y
        if z       != None: self.z       = z
        if occ     != None: self.occ     = occ
        if beta    != None: self.beta    = beta



    def makeCopy(self):
        """
        making a copy of this atom
        """

        newAtom = Atom()
        newAtom.type = self.type
        newAtom.numb = self.numb
        newAtom.name = self.name
        newAtom.element = self.element
        newAtom.resName = self.resName
        newAtom.resNumb = self.resNumb
        newAtom.chainID = self.chainID
        newAtom.x = self.x
        newAtom.y = self.y
        newAtom.z = self.z
        newAtom.occ = self.occ
        newAtom.beta = self.beta
        newAtom.terminal = self.terminal
        newAtom.residue_label = self.residue_label
        newAtom.icode = self.icode



        return newAtom


    def make_input_line(self):
        # making string
        group = '-'
        model_pka = '-'
        if self.group:
            group = self.group.type
            if self.terminal == 'C-':
                group = 'C-' ## circumventing C-/COO parameter unification

            if self.group.titratable:
                model_pka = '%6.2f'%self.group.model_pka


        str  = "%-6s%5d %s %s%2s%4d%12.3lf%8.3lf%8.3lf%6s%6s \n" % (self.type.upper(),
                                                                    self.numb,
                                                                    propka.lib.makeTidyAtomLabel(self.name,
                                                                                                 self.element),
                                                                    self.resName,
                                                                    self.chainID,
                                                                    self.resNumb,
                                                                    self.x,
                                                                    self.y,
                                                                    self.z,
                                                                    group,
                                                                    model_pka)




        return str

    def make_conect_line(self):
        res = 'CONECT%5d'%self.numb

        # extract and sort numbers of bonded residues
        bonded = []
        for atom in self.bonded_atoms:
            bonded.append(atom.numb)
        bonded.sort()

        # write 'em out
        for b in bonded:
            res += '%5d'%b
        res += '\n'
        return res


    def get_input_parameters(self):
        """ Method for getting the input parameters stored in the
        occupancy and b-factor fields in input files"""

        # Set the group type
        if self.occ != '-':
            # make sure to set the terminal
            if self.occ in ['N+','C-']:
                self.terminal = self.occ

            # save the ligand group charge
            if self.occ =='BLG':
                self.charge=+1
            elif self.occ =='ALG':
                self.charge=-1

            # generic ions
            if self.occ in ['1P','2P','1N','2N']:
                self.resName=self.occ
                self.occ='Ion'

            # correct the group type
            self.occ = self.occ.replace('N+','Nterm')
            self.occ = self.occ.replace('C-','Cterm')
            self.occ = self.occ.replace('ION','Ion')
            self.occ = self.occ.replace('ALG','titratable_ligand')
            self.occ = self.occ.replace('BLG','titratable_ligand')
            self.occ = self.occ.replace('LG','non_titratable_ligand')

            # try to initialise the group
            try:
                exec('self.group = propka.group.%s_group(self)'%self.occ)
            except:
                raise Exception('%s in input_file is not recognized as a group'%self.occ)

        # set the model pKa value
        if self.beta != '-':
            self.group.model_pka = float(self.beta)
            self.group.model_pka_set = True

        # set occ and beta to standard values
        self.occ = '1.00'
        self.beta = '0.00'



        return



    def make_pdb_line(self):

        # making string
        str  = "%-6s%5d %s %s%2s%4d%12.3lf%8.3lf%8.3lf%6s%6s\n" % (self.type.upper(),
                                                                   self.numb,
                                                                   propka.lib.makeTidyAtomLabel(self.name,
                                                                                                self.element),
                                                                   self.resName,
                                                                   self.chainID,
                                                                   self.resNumb,
                                                                   self.x,
                                                                   self.y,
                                                                   self.z,
                                                                   self.occ,
                                                                   self.beta)


        return str



    def make_mol2_line(self,id):
#1      S1     3.6147     2.0531     1.4795     S.3     1       noname  -0.1785
        # making string
        str  = "%-4d %-4s %10.4f %10.4f %10.4f %6s %6d %10s %10.4f\n" % (id,
                                                                         propka.lib.makeTidyAtomLabel(self.name,
                                                                                                      self.element),
                                                                         self.x,
                                                                         self.y,
                                                                         self.z,
                                                                         self.sybyl_type.replace('-',''),
                                                                         self.resNumb,
                                                                         self.resName,
                                                                         0.0)#self.charge)


        return str



    def makePDBLine(self,
                    numb    = None,
                    name    = None,
                    resName = None,
                    chainID = None,
                    resNumb = None,
                    x       = None,
                    y       = None,
                    z       = None,
                    occ     = None,
                    beta    = None):
        """
        returns a pdb ATOM-line for various purposes;
        specifying arguments over-writes.
        """
        if numb    == None: numb    = self.numb
        if name    == None: name    = self.name
        if resName == None: resName = self.resName
        if chainID == None: chainID = self.chainID
        if resNumb == None: resNumb = self.resNumb
        if x       == None: x       = self.x
        if y       == None: y       = self.y
        if z       == None: z       = self.z
        if occ     == None: occ     = self.occ
        if beta    == None: beta    = self.beta

        # making string
        str  = "ATOM "
        str += "%6d" % (numb)
        str += " %s" % (propka.lib.makeTidyAtomLabel(name,self.element))
        str += " %s" % (resName)
        str += "%2s" % (chainID)
        str += "%4d" % (resNumb)
        str += "%12.3lf" % (x)
        str += "%8.3lf" % (y)
        str += "%8.3lf" % (z)
        str += "%6.2lf" % (occ)
        str += "%6.2lf" % (beta)
        str += '\n'
        return str


    def getTidyLabel(self):
        """
        Returns a 'tidier' atom label for printing the new pdbfile
        """
        return propka.lib.makeTidyAtomLabel(self.name,self.element)


    def __str__(self):
        return '%5d-%4s %5d-%3s (%1s) [%8.3f %8.3f %8.3f] %s' %(self.numb, self.name, self.resNumb, self.resName, self.chainID, self.x, self.y, self.z,self.element)


#    def get_element(self):
#        """ try to extract element if not already done"""
#        if self.element == '':
#            self.element = self.name.strip(string.digits)[0]
#        return self.element


    def set_residue(self, residue):
        """ Makes a references to the parent residue"""
        if self.residue == None:
            self.residue = residue




