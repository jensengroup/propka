
from __future__ import division
from __future__ import print_function

import pickle,sys,os,math,propka.calculations

import pkg_resources
from propka.lib import info, warning

class bondmaker:
    def __init__(self):

        # predefined bonding distances
        self.distances = {
            'S-S':2.5,
            'F-F':1.7}

        self.distances_squared = {}
        for k in self.distances.keys():
            self.distances_squared[k]=self.distances[k]*self.distances[k]

        self.H_dist = 1.5;
        self.default_dist = 2.0;

        self.H_dist_squared  = self.H_dist * self.H_dist
        self.default_dist_squared = self.default_dist * self.default_dist

        self.max_sq_distance = max(list(self.distances_squared.values())+[self.default_dist_squared])

        # protein bonding data
        self.data_file_name = pkg_resources.resource_filename(__name__, 'protein_bonds.dat')

        data = open(self.data_file_name,'rb')
        self.protein_bonds = pickle.load(data)
        data.close()


        self.intra_residue_backbone_bonds = {'N': ['CA'],
                                             'CA':['N','C'],
                                             'C': ['CA','O'],
                                             'O': ['C']}

        self.number_of_pi_electrons_in_bonds_in_backbone = {'C':1,
                                                            'O':1}

        self.number_of_pi_electrons_in_conjugate_bonds_in_backbone = {'N':1}

        self.number_of_pi_electrons_in_bonds_in_sidechains = {'ARG-CZ' :1,
                                                              'ARG-NH1':1,
                                                              'ASN-OD1':1,
                                                              'ASN-CG' :1,
                                                              'ASP-OD1':1,
                                                              'ASP-CG' :1,
                                                              'GLU-OE1':1,
                                                              'GLU-CD' :1,
                                                              'GLN-OE1':1,
                                                              'GLN-CD' :1,
                                                              'HIS-CG' :1,
                                                              'HIS-CD2':1,
                                                              'HIS-ND1':1,
                                                              'HIS-CE1':1,
                                                              'PHE-CG' :1,
                                                              'PHE-CD1':1,
                                                              'PHE-CE1':1,
                                                              'PHE-CZ' :1,
                                                              'PHE-CE2':1,
                                                              'PHE-CD2':1,
                                                              'TRP-CG' :1,
                                                              'TRP-CD1':1,
                                                              'TRP-CE2':1,
                                                              'TRP-CD2':1,
                                                              'TRP-CE3':1,
                                                              'TRP-CZ3':1,
                                                              'TRP-CH2':1,
                                                              'TRP-CZ2':1,
                                                              'TYR-CG' :1,
                                                              'TYR-CD1':1,
                                                              'TYR-CE1':1,
                                                              'TYR-CZ' :1,
                                                              'TYR-CE2':1,
                                                              'TYR-CD2':1}


        self.number_of_pi_electrons_in_conjugate_bonds_in_sidechains = {'ARG-NE' :1,
                                                                        'ARG-NH2':1,
                                                                        'ASN-ND2':1,
                                                                        'GLN-NE2':1,
                                                                        'HIS-NE2':1,
                                                                        'TRP-NE1':1}

        self.number_of_pi_electrons_in_bonds_ligands = {'C.ar':1,
                                                        'N.pl3':0,
                                                        'C.2':1,
                                                        'O.2':1,
                                                        'O.co2':1,
                                                        'N.ar':1,
                                                        'C.1':2,
                                                        'N.1':2}

        self.number_of_pi_electrons_in_conjugate_bonds_in_ligands = {'N.am':1,'N.pl3':1}


        self.backbone_atoms = list(self.intra_residue_backbone_bonds.keys())

        self.terminal_oxygen_names = ['OXT','O\'\'']


        return


    def find_bonds_for_protein(self, protein):
        """ Finds bonds proteins based on the way atoms
        normally bond in proteins"""

        info('++++ Side chains ++++')
        # side chains
        for chain in protein.chains:
            for residue in chain.residues:
                if residue.resName.replace(' ','') not in ['N+','C-']:
                    self.find_bonds_for_side_chain(residue.atoms)

        info('++++ Backbones ++++')
        # backbone
        last_residues = []
        for chain in protein.chains:
            for i in range(1,len(chain.residues)):
                if chain.residues[i-1].resName.replace(' ','') not in ['N+','C-']:
                    if chain.residues[i].resName.replace(' ','') not in ['N+','C-']:
                        self.connect_backbone(chain.residues[i-1], chain.residues[i])
                        last_residues.append(chain.residues[i])

        info('++++ terminal oxygen ++++')
        # terminal OXT
        for last_residue in last_residues:
            self.find_bonds_for_terminal_oxygen(last_residue)

        info('++++ cysteines ++++')
        # Cysteines
        for chain in protein.chains:
            for i in range(0,len(chain.residues)):
                if chain.residues[i].resName == 'CYS':
                    for j in range(0,len(chain.residues)):
                        if chain.residues[j].resName == 'CYS' and j != i:
                            self.check_for_cysteine_bonds(chain.residues[i],
                                                          chain.residues[j])
        return

    def check_for_cysteine_bonds(self, cys1, cys2):
        for atom1 in cys1.atoms:
            if atom1.name == 'SG':
                for atom2 in cys2.atoms:
                    if atom2.name == 'SG':
                        if propka.calculations.squared_distance(atom1,atom2) < self.SS_dist_squared:
                            self.make_bond(atom1, atom2)


        return

    def find_bonds_for_terminal_oxygen(self, residue):
        for atom1 in residue.atoms:
            if atom1.name in self.terminal_oxygen_names:
                for atom2 in residue.atoms:
                    if atom2.name == 'C':
                        self.make_bond(atom1, atom2)

        return


    def connect_backbone(self, residue1, residue2):
        """ Sets up bonds in the backbone """
        # residue 1
        self.find_bonds_for_residue_backbone(residue1)

        # residue 2
        self.find_bonds_for_residue_backbone(residue2)

        # inter-residue bond
        for atom1 in residue1.atoms:
            if atom1.name == 'C':
                for atom2 in residue2.atoms:
                    if atom2.name == 'N':
                        if propka.calculations.squared_distance(atom1,atom2) < self.default_dist_squared:
                            self.make_bond(atom1, atom2)

        return

    def find_bonds_for_residue_backbone(self, residue):
        for atom1 in residue.atoms:
            if atom1.name in list(self.number_of_pi_electrons_in_bonds_in_backbone.keys()):
                atom1.number_of_pi_electrons_in_double_and_triple_bonds = self.number_of_pi_electrons_in_bonds_in_backbone[atom1.name]
            if atom1.name in list(self.number_of_pi_electrons_in_conjugate_bonds_in_backbone.keys()) and len(atom1.bonded_atoms)>1: # last part to avoid including N-term
                atom1.number_of_pi_electrons_in_conjugate_double_and_triple_bonds = self.number_of_pi_electrons_in_conjugate_bonds_in_backbone[atom1.name]

            if atom1.name in self.backbone_atoms:
                for atom2 in residue.atoms:
                    if atom2.name in self.intra_residue_backbone_bonds[atom1.name]:
                        self.make_bond(atom1, atom2)

        return


    def find_bonds_for_side_chain(self, atoms):
        """ Finds bonds for a side chain """
        for atom1 in atoms:

            key = '%s-%s'%(atom1.resName,atom1.name)
            if key in list(self.number_of_pi_electrons_in_bonds_in_sidechains.keys()):
                atom1.number_of_pi_electrons_in_double_and_triple_bonds = self.number_of_pi_electrons_in_bonds_in_sidechains[key]
            if key in list(self.number_of_pi_electrons_in_conjugate_bonds_in_sidechains.keys()):
                atom1.number_of_pi_electrons_in_conjugate_double_and_triple_bonds = self.number_of_pi_electrons_in_conjugate_bonds_in_sidechains[key]

            if not atom1.name in self.backbone_atoms:
                if not atom1.name in self.terminal_oxygen_names:
                    for atom2 in atoms:
                        if atom2.name in self.protein_bonds[atom1.resName][atom1.name]:
                            self.make_bond(atom1,atom2)

        return


    def find_bonds_for_ligand(self, ligand):
        """ Finds bonds for all atoms in the molecule """
        # identify bonding atoms
        self.find_bonds_for_atoms(ligand.atoms)
        self.add_pi_electron_table_info(ligand.atoms)

        return

    def add_pi_electron_table_info(self, atoms):
        # apply table information on pi-electrons
        for atom in atoms:
            # for ligands
            if atom.type == 'hetatm':
                if atom.sybyl_type in self.number_of_pi_electrons_in_bonds_ligands.keys():
                    atom.number_of_pi_electrons_in_double_and_triple_bonds = self.number_of_pi_electrons_in_bonds_ligands[atom.sybyl_type]
                if atom.sybyl_type in self.number_of_pi_electrons_in_conjugate_bonds_in_ligands.keys():
                    atom.number_of_pi_electrons_in_conjugate_double_and_triple_bonds = self.number_of_pi_electrons_in_conjugate_bonds_in_ligands[atom.sybyl_type]

            # for protein
            if atom.type == 'atom':
                key = '%s-%s'%(atom.resName,atom.name)
                if key in list(self.number_of_pi_electrons_in_bonds_in_sidechains.keys()):
                    atom.number_of_pi_electrons_in_double_and_triple_bonds = self.number_of_pi_electrons_in_bonds_in_sidechains[key]
                if key in list(self.number_of_pi_electrons_in_conjugate_bonds_in_sidechains.keys()):
                    atom.number_of_pi_electrons_in_conjugate_double_and_triple_bonds = self.number_of_pi_electrons_in_conjugate_bonds_in_sidechains[key]

                if atom.name in list(self.number_of_pi_electrons_in_bonds_in_backbone.keys()):
                    atom.number_of_pi_electrons_in_double_and_triple_bonds = self.number_of_pi_electrons_in_bonds_in_backbone[atom.name]
                if atom.name in list(self.number_of_pi_electrons_in_conjugate_bonds_in_backbone.keys()) and len(atom.bonded_atoms)>1: # last part to avoid including N-term
                    atom.number_of_pi_electrons_in_conjugate_double_and_triple_bonds = self.number_of_pi_electrons_in_conjugate_bonds_in_backbone[atom.name]

        return


    def find_bonds_for_protein_by_distance(self, molecule):
        """ Finds bonds for all atoms in the molecule """
        #self.find_bonds_for_protein(molecule)
        atoms = []
        for chain in molecule.chains:
            for residue in chain.residues:
                if residue.resName.replace(' ','') not in ['N+','C-']:
                    for atom in residue.atoms:
                        atoms.append(atom)

        self.find_bonds_for_atoms_using_boxes(atoms) #####
        #self.find_bonds_for_atoms(atoms) #####
        return atoms

    def find_bonds_for_atoms(self, atoms):
        """ Finds all bonds for a list of atoms"""
        no_atoms = len(atoms)

        for i in range(no_atoms):
            for j in range(i+1,no_atoms):

                if atoms[i] in atoms[j].bonded_atoms:
                    continue

                if self.check_distance(atoms[i], atoms[j]):
                    self.make_bond(atoms[i],atoms[j])
                    # di-sulphide bonds
                    if atoms[i].element == 'S' and atoms[j].element == 'S':
                        atoms[i].cysteine_bridge = True
                        atoms[j].cysteine_bridge = True


        return


    def check_distance(self, atom1, atom2):
        sq_dist = propka.calculations.squared_distance(atom1, atom2)

        if sq_dist > self.max_sq_distance:
            return False

        key = '%s-%s'%(atom1.element,atom2.element)
        h_count = key.count('H')

        if sq_dist <  self.H_dist_squared and h_count==1:
            return True
        if sq_dist < self.default_dist_squared and h_count==0:
            return True
        if key in self.distances_squared.keys():
            if sq_dist < self.distances_squared[key]:
                return True



        return False

    def find_bonds_for_molecules_using_boxes(self, molecules):
        """ Finds all bonds for a molecular container"""
        for name in molecules.conformation_names:
            self.find_bonds_for_atoms_using_boxes(molecules.conformations[name].atoms)
            #for chain in molecules.conformations[name].chains:
            #    self.find_bonds_for_atoms_using_boxes(molecules.conformations[name].get_chain(chain))

        return

    def add_pi_electron_information(self, molecules):
        for name in molecules.conformation_names:
            self.add_pi_electron_table_info(molecules.conformations[name].atoms)
        return


    def find_bonds_for_atoms_using_boxes(self, atoms):
        """ Finds all bonds for a list of atoms"""

        box_size = 2.5

        # find min and max coordinates
        xmin = 1e6; xmax = -1e6
        ymin = 1e6; ymax = -1e6
        zmin = 1e6; zmax = -1e6
        for atom in atoms:
            if atom.x > xmax:
                xmax = atom.x
            if atom.y > ymax:
                ymax = atom.y
            if atom.z > zmax:
                zmax = atom.z
            if atom.x < xmin:
                xmin = atom.x
            if atom.y < ymin:
                ymin = atom.y
            if atom.z < zmin:
                zmin = atom.z

        xlen = xmax-xmin
        ylen = ymax-ymin
        zlen = zmax-zmin

        #info('x range: [%6.2f;%6.2f] %6.2f'%(xmin,xmax,xlen))
        #info('y range: [%6.2f;%6.2f] %6.2f'%(ymin,ymax,ylen))
        #info('z range: [%6.2f;%6.2f] %6.2f'%(zmin,zmax,zlen))

        # how many boxes do we need in each dimension?
        # NOTE: math.ceil() returns an int in python3 and a float in
        # python2, so we need to convert it to an int for range() to work in
        # both versions. See PEP 3141.
        self.no_box_x = max(1, int(math.ceil(xlen/box_size)))
        self.no_box_y = max(1, int(math.ceil(ylen/box_size)))
        self.no_box_z = max(1, int(math.ceil(zlen/box_size)))

        #info('No. box x: %6.2f'%self.no_box_x)
        #info('No. box y: %6.2f'%self.no_box_y)
        #info('No. box z: %6.2f'%self.no_box_z)

        # initialize boxes
        self.boxes = {}
        for x in range(self.no_box_x):
            for y in range(self.no_box_y):
                for z in range(self.no_box_z):
                    self.boxes[(x,y,z)] = []

        # put atoms into boxes
        for atom in atoms:
            x = math.floor((atom.x-xmin)/box_size)
            y = math.floor((atom.y-ymin)/box_size)
            z = math.floor((atom.z-zmin)/box_size)
            self.put_atom_in_box(x,y,z,atom)

        # assign bonds
        for key, value in self.boxes.items():
            self.find_bonds_for_atoms(value)


        return

    def put_atom_in_box(self,x,y,z,atom):
        # atom in the x,y,z box and the up to 7 neighboring boxes on
        # one side of the x,y,z box in each dimension

        for bx in [x,x+1]:
            for by in [y,y+1]:
                for bz in [z,z+1]:
                    key = (bx,by,bz)
                    try:
                        self.boxes[key].append(atom)
                    except KeyError:
                        # No box exists for this coordinate
                        pass

                        #info(atom,'->',key,':',len(self.boxes[key]))

        return

    def has_bond(self, atom1, atom2):
        if atom1 in atom2.bonded_atoms or atom2 in atom1.bonded_atoms:
            return True
        return False


    def make_bond(self, atom1, atom2):
        """ Makes a bond between atom1 and atom2 """
        if atom1 == atom2:
            return

        #info('making bond for',atom1,atom2)
        if not atom1 in atom2.bonded_atoms:
            atom2.bonded_atoms.append(atom1)

        if not atom2 in atom1.bonded_atoms:
            atom1.bonded_atoms.append(atom2)

        return


    def generate_protein_bond_dictionary(self, atoms):

        for atom in atoms:
            for bonded_atom in atom.bonded_atoms:
                resi_i = atom.resName
                name_i = atom.name
                resi_j = bonded_atom.resName
                name_j = bonded_atom.name

                if not name_i in self.backbone_atoms or\
                        not name_j in self.backbone_atoms:
                    if not name_i in self.terminal_oxygen_names and\
                            not name_j in self.terminal_oxygen_names:

                        if not resi_i in list(self.protein_bonds.keys()):
                            self.protein_bonds[resi_i] = {}
                        if not name_i in self.protein_bonds[resi_i]:
                            self.protein_bonds[resi_i][name_i] = []

                        if not name_j in self.protein_bonds[resi_i][name_i]:
                            self.protein_bonds[resi_i][name_i].append(name_j)


                        if not resi_j in list(self.protein_bonds.keys()):
                            self.protein_bonds[resi_j] = {}
                        if not name_j in self.protein_bonds[resi_j]:
                            self.protein_bonds[resi_j][name_j] = []

                        if not name_i in self.protein_bonds[resi_j][name_j]:
                            self.protein_bonds[resi_j][name_j].append(name_i)


        return


if __name__ == '__main__':
    # If called directly, set up protein bond dictionary
    import protein, pdb, sys,os
    arguments = sys.argv
    if len(arguments) != 2:
        info('Usage: bonds.py <pdb_file>')
        sys.exit(0)

    filename = arguments[1]
    if not os.path.isfile(filename):
        info('Error: Could not find \"%s\"' % filename)
        sys.exit(1)

    pdblist = pdb.readPDB(filename)
    my_protein = protein.Protein(pdblist,'test.pdb')

    for chain in my_protein.chains:
        for residue in chain.residues:
            residue.atoms = [atom for atom in residue.atoms if atom.element != 'H']

    b = bondmaker()
    #b.protein_bonds = {}
    atoms = b.find_bonds_for_protein_by_distance(my_protein)
    #    b.generate_protein_bond_dictionary(atoms)

    #file = open(b.data_file_name,'wb')
    #pickle.dump(b.protein_bonds, file)
    #file.close()
