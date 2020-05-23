"""PROPKA representation of bonds."""
import math
import json
import pkg_resources
import propka.calculations
# TODO - replace the info/warning imports with logging functionality
from propka.lib import info


# TODO - should these constants be defined higher up in the module?
# TODO - I don't know what some of these constants mean
DISULFIDE_DISTANCE = 2.5
FLUORIDE_DISTANCE = 1.7
HYDROGEN_DISTANCE = 1.5
DEFAULT_DISTANCE = 2.0
BOX_SIZE = 2.5
POS_MAX = 1e6


class BondMaker:
    """Makes bonds?

    TODO - the documentation for this class needs to be improved.
    """
    def __init__(self):
        # predefined bonding distances
        self.distances = {'S-S' : DISULFIDE_DISTANCE, 'F-F' : FLUORIDE_DISTANCE}
        self.distances_squared = {}
        for key in self.distances:
            self.distances_squared[key] = self.distances[key] * self.distances[key]
        self.h_dist = HYDROGEN_DISTANCE
        self.default_dist = DEFAULT_DISTANCE
        self.h_dist_squared = self.h_dist * self.h_dist
        self.default_dist_squared = self.default_dist * self.default_dist
        distances = list(self.distances_squared.values()) + [self.default_dist_squared]
        self.max_sq_distance = max(distances)
        # protein bonding data
        self.data_file_name = pkg_resources.resource_filename(__name__, 'protein_bonds.json')
        with open(self.data_file_name, 'rt') as json_file:
            self.protein_bonds = json.load(json_file)
        self.intra_residue_backbone_bonds = {'N': ['CA'], 'CA': ['N', 'C'],
                                             'C': ['CA', 'O'], 'O': ['C']}
        self.num_pi_elec_bonds_backbone = {'C': 1, 'O': 1}
        self.num_pi_elec_conj_bonds_backbone = {'N': 1}
        self.num_pi_elec_bonds_sidechains = {'ARG-CZ' : 1, 'ARG-NH1': 1,
                                             'ASN-OD1': 1, 'ASN-CG' : 1,
                                             'ASP-OD1': 1, 'ASP-CG' : 1,
                                             'GLU-OE1': 1, 'GLU-CD' : 1,
                                             'GLN-OE1': 1, 'GLN-CD' : 1,
                                             'HIS-CG' : 1, 'HIS-CD2': 1,
                                             'HIS-ND1': 1, 'HIS-CE1': 1,
                                             'PHE-CG' : 1, 'PHE-CD1': 1,
                                             'PHE-CE1': 1, 'PHE-CZ' : 1,
                                             'PHE-CE2': 1, 'PHE-CD2': 1,
                                             'TRP-CG' : 1, 'TRP-CD1': 1,
                                             'TRP-CE2': 1, 'TRP-CD2': 1,
                                             'TRP-CE3': 1, 'TRP-CZ3': 1,
                                             'TRP-CH2': 1, 'TRP-CZ2': 1,
                                             'TYR-CG' : 1, 'TYR-CD1': 1,
                                             'TYR-CE1': 1, 'TYR-CZ' : 1,
                                             'TYR-CE2': 1, 'TYR-CD2': 1}
        self.num_pi_elec_conj_bonds_sidechains = {'ARG-NE': 1, 'ARG-NH2': 1,
                                                  'ASN-ND2': 1, 'GLN-NE2': 1,
                                                  'HIS-NE2': 1, 'TRP-NE1': 1}
        self.num_pi_elec_bonds_ligands = {'C.ar': 1, 'N.pl3': 0, 'C.2': 1,
                                          'O.2': 1, 'O.co2': 1, 'N.ar': 1,
                                          'C.1': 2, 'N.1': 2}
        self.num_pi_elec_conj_bonds_ligands = {'N.am': 1, 'N.pl3': 1}
        self.backbone_atoms = list(self.intra_residue_backbone_bonds.keys())
        self.terminal_oxygen_names = ['OXT', 'O\'\'']
        self.boxes = {}
        self.num_box_x = None
        self.num_box_y = None
        self.num_box_z = None

    def find_bonds_for_protein(self, protein):
        """Finds bonds proteins based on the way atoms normally bond in proteins.

        Args:
            protein:  the protein to search for bonds
        """
        info('++++ Side chains ++++')
        # side chains
        for chain in protein.chains:
            for residue in chain.residues:
                if residue.res_name.replace(' ', '') not in ['N+', 'C-']:
                    self.find_bonds_for_side_chain(residue.atoms)
        info('++++ Backbones ++++')
        # backbone
        last_residues = []
        for chain in protein.chains:
            for i in range(1, len(chain.residues)):
                if chain.residues[i-1].res_name.replace(' ', '') not in ['N+', 'C-']:
                    if chain.residues[i].res_name.replace(' ', '') not in ['N+', 'C-']:
                        self.connect_backbone(chain.residues[i-1], chain.residues[i])
                        last_residues.append(chain.residues[i])
        info('++++ terminal oxygen ++++')
        # terminal OXT
        for last_residue in last_residues:
            self.find_bonds_for_terminal_oxygen(last_residue)
        info('++++ cysteines ++++')
        # Cysteines
        for chain in protein.chains:
            for i in range(0, len(chain.residues)):
                if chain.residues[i].res_name == 'CYS':
                    for j in range(0, len(chain.residues)):
                        if chain.residues[j].res_name == 'CYS' and j != i:
                            self.check_for_cysteine_bonds(chain.residues[i],
                                                          chain.residues[j])

    def check_for_cysteine_bonds(self, cys1, cys2):
        """Looks for potential bonds between two cysteines.

        Args:
            cys1:  one of the cysteines to check
            cys1:  one of the cysteines to check
        """
        for atom1 in cys1.atoms:
            if atom1.name == 'SG':
                for atom2 in cys2.atoms:
                    if atom2.name == 'SG':
                        dist = propka.calculations.squared_distance(atom1, atom2)
                        # TODO - is SS_dist_squared an attribute of this class?
                        if dist < self.SS_dist_squared:
                            self.make_bond(atom1, atom2)

    def find_bonds_for_terminal_oxygen(self, residue):
        """Look for bonds for terminal oxygen.

        Args:
            residue - test residue
        """
        for atom1 in residue.atoms:
            if atom1.name in self.terminal_oxygen_names:
                for atom2 in residue.atoms:
                    if atom2.name == 'C':
                        self.make_bond(atom1, atom2)

    # TODO - stopped here.
    def connect_backbone(self, residue1, residue2):
        """Sets up bonds in the backbone

        Args:
            residue1:  first residue to connect
            residue2:  second residue to connect
        """
        self.find_bonds_for_residue_backbone(residue1)
        self.find_bonds_for_residue_backbone(residue2)
        for atom1 in residue1.atoms:
            if atom1.name == 'C':
                for atom2 in residue2.atoms:
                    if atom2.name == 'N':
                        if propka.calculations.squared_distance(atom1, atom2) \
                            < self.default_dist_squared:
                            self.make_bond(atom1, atom2)

    def find_bonds_for_residue_backbone(self, residue):
        """Find bonds for this residue's backbone.

        Args:
            residue:  reside to search for backbone bonds.
        """
        for atom1 in residue.atoms:
            if atom1.name in list(self.num_pi_elec_bonds_backbone.keys()):
                atom1.num_pi_elec_2_3_bonds \
                    = self.num_pi_elec_bonds_backbone[atom1.name]
            if atom1.name in \
                list(self.num_pi_elec_conj_bonds_backbone.keys()) \
                    and len(atom1.bonded_atoms) > 1: # last part to avoid including N-term
                atom1.num_pi_elec_conj_2_3_bonds \
                    = self.num_pi_elec_conj_bonds_backbone[atom1.name]

            if atom1.name in self.backbone_atoms:
                for atom2 in residue.atoms:
                    if atom2.name in self.intra_residue_backbone_bonds[atom1.name]:
                        self.make_bond(atom1, atom2)

    def find_bonds_for_side_chain(self, atoms):
        """Finds bonds for a side chain.

        Args:
            atoms:  list of atoms to check for bonds
        """
        for atom1 in atoms:
            key = '%s-%s' % (atom1.res_name, atom1.name)
            if key in list(self.num_pi_elec_bonds_sidechains.keys()):
                atom1.num_pi_elec_2_3_bonds \
                    = self.num_pi_elec_bonds_sidechains[key]
            if key in list(self.num_pi_elec_conj_bonds_sidechains.keys()):
                atom1.num_pi_elec_conj_2_3_bonds \
                    = self.num_pi_elec_conj_bonds_sidechains[key]

            if not atom1.name in self.backbone_atoms:
                if not atom1.name in self.terminal_oxygen_names:
                    for atom2 in atoms:
                        if atom2.name in self.protein_bonds[atom1.res_name][atom1.name]:
                            self.make_bond(atom1, atom2)

    def find_bonds_for_ligand(self, ligand):
        """Finds bonds for all atoms in the ligand molecule

        Args:
            ligand:  ligand molecule to search for bonds
        """
        # identify bonding atoms
        self.find_bonds_for_atoms(ligand.atoms)
        self.add_pi_electron_table_info(ligand.atoms)

    def add_pi_electron_table_info(self, atoms):
        """Add table information on pi electrons

        Args:
            atoms:  list of atoms for pi electron table information checking
        """
        # apply table information on pi-electrons
        for atom in atoms:
            # for ligands
            if atom.type == 'hetatm':
                if atom.sybyl_type in self.num_pi_elec_bonds_ligands.keys():
                    atom.num_pi_elec_2_3_bonds = self.num_pi_elec_bonds_ligands[atom.sybyl_type]
                if atom.sybyl_type in self.num_pi_elec_conj_bonds_ligands.keys():
                    atom.num_pi_elec_conj_2_3_bonds \
                        = self.num_pi_elec_conj_bonds_ligands[atom.sybyl_type]
            # for protein
            if atom.type == 'atom':
                key = '%s-%s' % (atom.res_name, atom.name)
                if key in list(self.num_pi_elec_bonds_sidechains.keys()):
                    atom.num_pi_elec_2_3_bonds = self.num_pi_elec_bonds_sidechains[key]
                if key in list(self.num_pi_elec_conj_bonds_sidechains.keys()):
                    atom.num_pi_elec_conj_2_3_bonds = self.num_pi_elec_conj_bonds_sidechains[key]

                if atom.name in list(self.num_pi_elec_bonds_backbone.keys()):
                    atom.num_pi_elec_2_3_bonds = self.num_pi_elec_bonds_backbone[atom.name]
                if atom.name in list(self.num_pi_elec_conj_bonds_backbone.keys()) \
                    and len(atom.bonded_atoms) > 1:
                    # last part to avoid including N-term
                    atom.num_pi_elec_conj_2_3_bonds \
                        = self.num_pi_elec_conj_bonds_backbone[atom.name]

    def find_bonds_for_protein_by_distance(self, molecule):
        """Finds bonds for all atoms in the molecule.

        Args:
            molecule:  molecule in which to find bonds.
        Returns:
            list of atoms
        """
        #self.find_bonds_for_protein(molecule)
        atoms = []
        for chain in molecule.chains:
            for residue in chain.residues:
                if residue.res_name.replace(' ', '') not in ['N+', 'C-']:
                    for atom in residue.atoms:
                        atoms.append(atom)
        self.find_bonds_for_atoms_using_boxes(atoms)
        return atoms

    def find_bonds_for_atoms(self, atoms):
        """Finds all bonds for a list of atoms

        Args:
            atoms:  list of atoms in which to find bonds.
        """
        no_atoms = len(atoms)
        for i in range(no_atoms):
            for j in range(i+1, no_atoms):
                if atoms[i] in atoms[j].bonded_atoms:
                    continue
                if self.check_distance(atoms[i], atoms[j]):
                    self.make_bond(atoms[i], atoms[j])
                    # di-sulphide bonds
                    if atoms[i].element == 'S' and atoms[j].element == 'S':
                        atoms[i].cysteine_bridge = True
                        atoms[j].cysteine_bridge = True

    def check_distance(self, atom1, atom2):
        """Check distance between two atoms

        Args:
            atom1:  first atom for distance check
            atom2:  second atom for distance check
        Returns:
            True if within distance, False otherwise
        """
        sq_dist = propka.calculations.squared_distance(atom1, atom2)
        if sq_dist > self.max_sq_distance:
            return False
        key = '%s-%s' % (atom1.element, atom2.element)
        h_count = key.count('H')
        if sq_dist < self.h_dist_squared and h_count == 1:
            return True
        if sq_dist < self.default_dist_squared and h_count == 0:
            return True
        if key in self.distances_squared.keys():
            if sq_dist < self.distances_squared[key]:
                return True
        return False

    def find_bonds_for_molecules_using_boxes(self, molecules):
        """ Finds all bonds for a molecular container.

        Args:
            molecules:  list of molecules for finding bonds.
        """
        for name in molecules.conformation_names:
            self.find_bonds_for_atoms_using_boxes(molecules.conformations[name].atoms)

    def add_pi_electron_information(self, molecules):
        """Add pi electron information to a molecule.

        Args:
            molecules:  list of molecules for adding pi electron information.
        """
        for name in molecules.conformation_names:
            self.add_pi_electron_table_info(molecules.conformations[name].atoms)

    def find_bonds_for_atoms_using_boxes(self, atoms):
        """Finds all bonds for a list of atoms.

        Args:
            atoms:  list of atoms for finding bonds
        """
        box_size = BOX_SIZE
        xmin = POS_MAX
        xmax = -1.0 * POS_MAX
        ymin = POS_MAX
        ymax = -1.0 * POS_MAX
        zmin = POS_MAX
        zmax = -1.0 * POS_MAX
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
        xlen = xmax - xmin
        ylen = ymax - ymin
        zlen = zmax - zmin
        self.num_box_x = max(1, int(math.ceil(xlen/box_size)))
        self.num_box_y = max(1, int(math.ceil(ylen/box_size)))
        self.num_box_z = max(1, int(math.ceil(zlen/box_size)))
        self.boxes = {}
        for x in range(self.num_box_x):
            for y in range(self.num_box_y):
                for z in range(self.num_box_z):
                    self.boxes[(x, y, z)] = []
        for atom in atoms:
            x = math.floor((atom.x - xmin)/box_size)
            y = math.floor((atom.y - ymin)/box_size)
            z = math.floor((atom.z - zmin)/box_size)
            self.put_atom_in_box(x, y, z, atom)
        for value in self.boxes.values():
            self.find_bonds_for_atoms(value)

    def put_atom_in_box(self, x, y, z, atom):
        """Put an atom in a box.

        Args:
            x:  box x-coordinates
            y:  box y-coordinates
            z:  box z-coordinates
            atom:  the atom to place in a box
        """
        for box_x in [x, x+1]:
            for box_y in [y, y+1]:
                for box_z in [z, z+1]:
                    key = (box_x, box_y, box_z)
                    try:
                        self.boxes[key].append(atom)
                    except KeyError:
                        pass

    @staticmethod
    def has_bond(atom1, atom2):
        """Look for bond between two atoms.

        Args:
            atom1:  first atom to check
            atom2:  second atom to check
        Returns:
            True if there is a bond between atoms
        """
        if atom1 in atom2.bonded_atoms or atom2 in atom1.bonded_atoms:
            return True
        return False

    @staticmethod
    def make_bond(atom1, atom2):
        """Makes a bond between atom1 and atom2

        Args:
            atom1:  first atom to bond
            atom2:  second atom to bond
        """
        if atom1 == atom2:
            return
        if not atom1 in atom2.bonded_atoms:
            atom2.bonded_atoms.append(atom1)
        if not atom2 in atom1.bonded_atoms:
            atom1.bonded_atoms.append(atom2)

    def generate_protein_bond_dictionary(self, atoms):
        """Generate dictionary of protein bonds.

        Args:
            atoms:  list of atoms for bonding
        """
        for atom in atoms:
            for bonded_atom in atom.bonded_atoms:
                resi_i = atom.res_name
                name_i = atom.name
                resi_j = bonded_atom.res_name
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
