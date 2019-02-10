# -*- coding: utf-8 -*-
"""
force field Code for Masters project.
Author: Craig Melton
28/11/2017
"""



import periodic_table as pt
import numpy as np
import ringalgothem_final_2 as ring
import math
import re
#import simple_ringfinding_algothem  as ring
# Molecule and atom classes
import collections


class atom:
    """atom types contian the atomic type and the xyz coordinates, polariseably, 
    electronegativity"""
    def __init__(self, symbol, mass, position, groups,pdb):
        self.symbol = symbol
        self.mass = mass
        self.position = position
        if self.symbol != 'Q':
            try:
                self.max_valence = pt.max_valence[self.symbol]
                
            except KeyError:
                self.max_valence = 6
                
        elif self.symbol == 'Q':
            self.max_valence = 1
            
        else:
            self.max_valence = 0
        self.connectivity = []
        self.connectivity2 = []
        self.raduces = 0
        self.charge = 0
        self.ExtendedConnectivity = 0
        self.non_H_connections = 0
        self.Valence = 0
        self.ringfinder = 0
        self.ringfinder2 = 0
        self.ringfinder_list = []
        
        self.groups = groups
        self.pdb = pdb
        
        try:
            self.ID = int(pdb[6:11])
        except ValueError:
            self.ID = 0
        self.angels = []
        
        self.SpareValence = 0
        self.options = 0
        if len(self.pdb) > 0: 
            try:
                self.resSeq = int(self.pdb[22:26])
            except ValueError:
                self.resSeq = -1
        else:
            self.resSeq = -1
        
        
        
    def __str__(self):
        """what the atom class returns when printing"""
        return '({}, {}, {})'.format(self.symbol, self.mass, self.position)
    
    def __repr__(self):
        """evaluate atom class returns string"""
        return '({}, {}, {})'.format(self.symbol, self.mass, self.position)
    def __iter__(self):
        """iterable for atom in sets"""
        return '({}, {}, {})'.format(self.symbol, self.mass, self.position)
    def __key(self):
        """ """
        return (self.position[0],self.position[1],self.position[2])
    def __eq__(self, other):
        """ """
        if  isinstance(other, self.__class__) and other.__key() == self.__key():
            return True
        else:
            return False
    def __hash__(self):
        """hash number for atom"""
        return hash(self.__key())
   
class Bond:
    """defines bonds in a a molecule or system"""
    def __init__(self, length, atoms, bondorder=1):
        self.length = length
        self.atoms = atoms
        self.bondorder  = bondorder
        self.connections = []
        self.SP_list = []
        
    def __str__(self):
        """what the bond class returns when printing"""
        atms = []
        for atom in self.atoms:
            atms.append(atom.ID)
        atms.sort()  
        return '({}, {}, {})'.format(self.length, self.bondorder, atms)
    
    def __repr__(self):
        """evaluate bond class returns string"""
        atms = []
        for atom in self.atoms:
            atms.append(atom.ID)
        atms.sort()  
        return '({}, {}, {})'.format(self.length, self.bondorder, atms)#self.atoms)
    
    def __key(self):
        """ """
        (self.atoms).sort(key= lambda x: x.__hash__())
        atms = tuple(self.atoms)
          
        return (self.length, atms)
    def __eq__(self, other):
        """ """
        if  isinstance(other, self.__class__) and other.__key() == self.__key():
            return True
        else:
            return False
    def __hash__(self):
        """hash number for atom"""
        return hash(self.__key())

class Angel:
    """the angle class contains the atoms need to generate the angle and the angle
    found.""" 
    def __init__(self, angle, atoms):
        self.angle = angle
        self.atoms = atoms
    def __key(self):
        """ """
        (self.atoms).sort(key= lambda x: x.__hash__())
        atms = tuple(self.atoms)
        return (self.length, atms)
    def __eq__(self, other):
        """ """
        if  isinstance(other, self.__class__) and other.__key() == self.__key():
            return True
        else:
            return False
    def __hash__(self):
        """hash number for atom"""
        return hash(self.__key())
 
def bond_hack(molecule):
    """changes connectivity tuples into bonds of class Bond. for all atom 
    in the molecule"""
    for atom in molecule.atoms:
        connection2 = []
        for bond in atom.connectivity: 
           
            if atom.symbol != 'Q':
                connection2.append(Bond(bond[0], [atom, bond[1]]))
            else:
                connection2.append(Bond(bond[0], [atom, bond[1]]))
        atom.connectivity = connection2
#    for atom1 in molecule.atoms:
#        
#        k = 0
#        for bond in atom1.connectivity:
#            atom2 = bond.atoms[1]
#            try:
#                if pt.max_valence[atom1.symbol] - atom1.max_valence - k != 0:
#                    try:
#                        n = pt.max_valence[atom2.symbol] - atom2.max_valence
#                        bond.bondorder += n
#                        k += n
#                    except KeyError:
#                        n = 6 - atom2.max_valence
#                        bond.bondorder += n
#                        k += n
#                        
#            except KeyError:
#                if 6 - atom1.max_valence - k != 0:
#                    try:
#                        n = pt.max_valence[atom2.symbol] - atom2.max_valence
#                        bond.bondorder += n
#                        k += n
#                    except KeyError:
#                        n = 6 - atom2.max_valence
#                        bond.bondorder += n
#                        k += n
#    for atom in molecule.atoms:
#        if atom.symbol != 'Q':
#            connection2 = []
#            for bond in atom.connectivity: 
#                try:
#                    connection2.append(Bond(bond[0], [atom, bond[1]]))
#                except TypeError:
#                    continue

    bonds2 = []
    for atom in molecule.atoms:
        for bond2 in atom.connectivity:
            bonds2.append(bond2)
    
    molecule.bonds = list(set(bonds2))
    
    return None


        
     
class molecule:
    """molecule class containing data on the atoms in the molecule in the form 
    of the atom class, bonds, angles, dihedrals"""
    def __init__(self, atoms):
        self.atoms = atoms
        self.bonds = set()
        self.errors = []
        self.aromatic = []
        self.rings = []
        self.group_list = []
        self.angels = []
        self.ring_errors = []
        self.carbonyl_group = []
        self.intersting_groups = []
        self.bonds_from_file = []
        self.ring_cycles = 0
        
        self.amino_acids_dic = {}
#        self.ringnumbers = [] 
        
    def distance_matrix(self):
        """creates distance matrix"""
        coords = []
        for atom in self.atoms:
            coords.append(atom.position)
        coords = np.array(coords)
        
        x2 = np.sum(coords**2,axis=1)
        y2 = x2[:,np.newaxis]
        xy = np.dot(coords,coords.T)
        r2 = x2 + y2 - 2*xy
        zeros = np.zeros_like(r2)
        r2m = np.maximum(r2,zeros)
        distance_matrix = np.sqrt(r2m)
        return distance_matrix

    def connectivity(self, n):    
        """Creates the distance matrix with numpy functions and arrays"""
        atoms = self.atoms 
        bonds = self.bonds
        errors = self.errors
        aromatic = self.aromatic
        
#        for atom in atoms:    
#            if atom.max_valence < 5:
#                for atom10 in atoms:
#                    length = np.linalg.norm(atom.position - atom10.position)
#                    
#                    if length < 2.80:
#                        (atom.connectivity).append((length, atom10))
        
        distancematrix = self.distance_matrix() 
        i = 0
        for atom1 in self.atoms:
            
            j = 0
            for atom2 in self.atoms:
                atom1.connectivity2.append((distancematrix[i][j],atom2))
                j += 1
            
            i += 1
        for atom1 in self.atoms:
            if atom1.max_valence < 5:
                sorted(atom1.connectivity2, key = lambda x: x[0])
                if atom1.symbol == 'N':
                    k = 0
                    for num in range(2):
                        if sorted(atom1.connectivity2, key = lambda x: x[0])[num+1][1].symbol == 'H':
                            k +=1
#                            if atom1.ID == 244 or atom1.ID == 209:
#                                print('ID {} atom1.max_valence {},\natom1.connectivity {}'.format(atom1.ID,atom1.max_valence,sorted(atom1.connectivity2, key = lambda x: x[0])[:4]))
                    
                    if k >= 2:
                        atom1.max_valence += 1
                        atom1.Valence += 1
                
                atom1.connectivity = (sorted(atom1.connectivity2, key = lambda x: x[0]))[1:atom1.max_valence+1]
                
        for atom1 in self.atoms:
            for bond in atom1.connectivity:
                if bond[0] > 2.80:
                    atom1.connectivity.remove(bond)
        
                n = 0
                
                        
                
        #print(atoms[-1].connectivity)
#        try:  # check for equil lengths from atom, atoms 
#            for atom2 in atoms:
#                connections = (atom2.connectivity)
#                atom2.connectivity2 = connections
#                sorted(atom2.connectivity2, key = lambda x: x[0])
#                connections = (sorted(connections, key = lambda x: x[0]))[1:atom2.max_valence+1]
#                
#                atom2.connectivity = connections
#        except TypeError:
#            print('atom in error {}\n'.format(atom2))
#            print('connections\n')
#            for atom21 in atom2.connectivity:
#                print(atom21)
        
        for atom3 in self.atoms:
            self.bonds.update(atom3.connectivity)
            
        bonds_removed = [0]
        while len(bonds_removed) > 0:
            bonds_removed = []
            for atom5 in self.atoms:
                
                for bond1 in atom5.connectivity:
                    i = 0
                    for bond2 in bond1[1].connectivity:
                        
                        if atom5 == bond2[1]:
                           i = 1
                           
                          
                    if i == 0:
        #                    length = len((atom5.connectivity))
                        try:
                            (atom5.connectivity).remove(bond1)
                            
                            bonds_removed.append(bond1)
                        except ValueError:
                            continue
        
        
        for atom5 in self.atoms:
            
            for atom6 in atom5.connectivity:
                i = 0
                for atom7 in atom6[1].connectivity:
                    
                    if atom5 == atom7[1]:
                       i += 1
                if i == 0:
                    (atom5.connectivity).remove(atom6)
                    
        
        
#        for atom8 in atoms:
#            try:
#                reduced_valence = pt.max_valence[atom8.symbol] - atom8.max_valence
#            except KeyError:
#                reduced_valence = 6 - atom8.max_valence
#            if reduced_valence == 0:
#                continue
#            i = 0
#            for atom9 in atom8.connectivity:
#                try:
#                    i += pt.max_valence[atom9[1].symbol] -  atom9[1].max_valence
#                except KeyError:
#                    i += 6 -  atom9[1].max_valence
#              
#            if reduced_valence != i and atom8.max_valence % 1 == 0:
#                errors.append(atom8)
#            elif atom8.max_valence % 1 != 0:
#                aromatic.append(atom8)
                
                
        for atom3 in self.atoms:
            self.bonds.update(atom3.connectivity)    
        self.bonds = list(set(self.bonds))
        return self.bonds, errors
    
    def metal_connections(self):
        """goes back and add the connectivities for the metals atoms with valence of 5 
        and above"""
        for atom1 in self.atoms:
            if len(atom1.connectivity) >= 5 and atom1.symbol != 'Q':
                
                atom1.connectivity = (sorted(atom1.connectivity2, key = lambda x: x[0]))[1:atom1.max_valence+1]
    
#    def higher_valance(self):
#        """take the 5 or 6 max valance atoms and fills in there connectivity"""
#        for atom in self.atoms:
#            
#            if atom.max_valence > 4:
#                for atom10 in self.atoms:
#                    length = np.linalg.norm(atom.position - atom10.position)
#                    
#                    if length < 2.80:
#                        (atom.connectivity).append((length, atom10))
#        return None
    
    
    def amino_acids(self,amino_acids):
        """takes atoms from molecule and then returns a list directery with the number of each atom type present"""
        amino_acids_dic1 = {}

        for amino_acid in amino_acids:

            amino_acids_dic1[amino_acid.amino_acid[1]] = 1 + amino_acids_dic1.get(amino_acid.amino_acid[1], 0)

        self.amino_acids_dic = amino_acids_dic1

        return None 
    
    
#    def amino_acids(self):
#        """takes atoms from molecule and then returns a list directery with the number of each atom type present"""
#        amino_acids_dic1 = {}
#        resSeq = 0
#        self.atoms.sort(key=lambda atom:atom.resSeq)
#        for atom in self.atoms:
#            if atom.resSeq > -1:
##            print(atom.pdb[22:26])
#                if int(atom.pdb[22:26]) != resSeq:
#                    amino_acids_dic1[atom.pdb[17:20]] = 1 + amino_acids_dic1.get(atom.pdb[17:20], 0)
#                    resSeq = int(atom.pdb[22:26])
#        self.amino_acids_dic = amino_acids_dic1
#        return None 
    
#    def Hydrogen_duplite_remove(self):
#        """if to atoms are to close together the one with the longer bond is 
#        removed along with its bond"""
#        for atom in self.atoms:
#            if atom.symbol == 'H':
#                for atom2 in self.atoms:
#                    if atom2.symbol == 'H':
#                        if atom.connectivity[0][0] < atom2.connectivity[0][0]:
#                            atom2.connectivity[0][1].connectivity.remove((atom2.connectivity[0][0],atom2))
#                            self.atoms.remove(atom2)
#                        else:
#                            atom.connectivity[0][1].connectivity.remove((atom.connectivity[0][0],atom))
#                            self.atoms.remove(atom)
#        for atom_H in self.atoms:
#            if atom_H.symbol == 'H':
#                if atom_H.symbol == atom_H.connectivity[0][1].symbol and atom_H.connectivity[0][0] < 0.6:
#                    if atom_H.connectivity2[1][0] < atom_H.connectivity[0][1].connectivity[1][0]:
#                        atom_H.connectivity[0][1].connectivity[0][1].connectivity.remove((atom_H.connectivity[0][0],atom_H.connectivity[0][1]))
#                        self.atoms.remove(atom_H.connectivity[0][1])
#                    else:
#                        atom_H.connectivity[0][1].connectivity.remove((atom_H.connectivity[0][0],atom_H))
#                        self.atoms.remove(atom_H)
#        return None
                    
                    
                    
    
    def non_H_connections_num(self):
        """take connections data from the atoms and find the number of non H 
        atoms connected to each atoms then changes non_H_connections porpertie
        of the atom"""
        atoms =  self.atoms
        
        
        for atom in atoms:
            if atom.symbol == ('H' or  'Q'):
                continue
            
            for atom2 in atom.connectivity:
                if atom2[1].symbol != ('H' or 'Q'):
                    atom.non_H_connections += 1
        
                    
    
    def Extendedconnections(self):
        """this function carries out the Morgan Algorithm is the molecule and 
        the atoms it contians to calculate the Extended connections of the atoms"""
        atoms = self.atoms
        
        
        for atom in atoms:
            try:
                num_of_doublebonds = pt.max_valence[atom.symbol] - atom.max_valence
            except KeyError:
                num_of_doublebonds = 6 - atom.max_valence
            
            atom.ExtendedConnectivity += atom.non_H_connections + num_of_doublebonds
            atom.ExtendedConnectivity *= 10
            if atom.symbol == 'C':
                atom.ExtendedConnectivity += 2
            elif atom.symbol == 'N':
                atom.ExtendedConnectivity += 3
            elif atom.symbol == 'O':
                atom.ExtendedConnectivity += 4

            EC_list = []
            for atom in atoms:
                EC_list.append(atom.ExtendedConnectivity)
                
            NEC = len(set(EC_list))
            NTEC = NEC+1
            i = 0
            while NEC < NTEC:
                
                if i < 1:
                    NEC = NTEC
                for atom in atoms:
                    if atom.symbol == 'H' or atom.symbol == 'Q':
                        continue
                    atom.ExtendedConnectivity *= 5
                connections = atom.connectivity
                for atom2 in connections:
                    if atom2[1].symbol != 'H' or atom2[1].symbol != 'Q':
                        atom.ExtendedConnectivity += atom2[1].ExtendedConnectivity
                TEC_list = []
                for atom in atoms:
                    TEC_list.append(atom.ExtendedConnectivity)
                NTEC = len(set(EC_list))
                i += 1

    def find_angles(self):
        """finding all the angles in a molecule"""
        for atom1 in self.atoms:
            atoms2 = connected_atom_bond(atom1.connectivity)
            i = 1
            for bond1 in atoms2:
                for bond2 in atoms2[i:]:
                    atoms = [atom1, bond1[1],bond2[1]]
                    length = np.linalg.norm(bond1[1].position - bond2[1].position)
                    cos = (bond1[0]^2 + bond2[0]^2 -length^2)/(2*bond1[0]*bond2[0])
                    angel = math.degrees(math.acos(cos))
                    atom1.angles.append(Angel(angel,atoms))
                    
                i += 1
        all_angels = []
        for atom1 in self.atoms:
            for angel in atom1.angels:
                all_angels.append(angel)
        self.angels = list(set(all_angels))
        return None




 

 
 
#Notes:
#coords is (Natom x 3) numpy array (so basically your coordinates list converted into numpy array form)
#the “r2m” line is just to ensure numerical stability – there is a chance that r2 contains “negative zeros” due to round-off error (e.g. -1.e-17) which cause the sqrt operator to fail. This line just ensures that these are set to exactly zero.
#



# loading in molecules atomic position data

def connected_atom_bond(connectivity):
    """from bond class returns the list of exteral atoms to the connectivity of
     the with the connectivity list"""
    atoms = []
    for bond1 in connectivity:
        atoms.append((bond1.length,bond1.atoms[1]))
    return connectivity

def read_xyzfile(filename): # week 1 
        """read in the atomic postions and symbols to a list of list containing the 
        atomic symbols and typle of carteasin coordinates of the atom"""
        Molecule = open(filename)
        
        
        n = Molecule.readline()
        n = int(n)
        
        Molecule.readline()
        
        atom_range = range(n)
        
        atom_names = []
        for num in atom_range:
            name = "atom_{}".format(num)
            atom_names.append(name)
        i = 0
        atom_list = []
        while i < n:
            data = Molecule.readline()
            data = data.split()
            
            atomic_data = []
            for value in data:
                try:
                    value = float(value)
                    atomic_data.append(value)
                except ValueError:
                    value = value.strip()
                    atomic_data.append(value)
                    
            data = atomic_data
            
            array = np.array((data[1], data[2], data[3]))
            mass = pt.masses[data[0]]
            atom_names[i] = atom(data[0], mass, array, ring.group([],0))
            atom_list.append(atom_names[i])
            i += 1
            

        
        return n, atom_list


def read_in_molecule_data(filename):
    """read in the atoms and molecular data for a molecule"""
    n, atoms = read_xyzfile(filename)
    Molecule = molecule(atoms)
    return Molecule, n

def read_pdbfile(file): # week 1 
        """read in the atomic postions and symbols to a list of list containing the 
        atomic symbols and typle of carteasin coordinates of the atom"""
        n = len(file)
        atom_list = []
        for atom1 in file:
            array = np.array([atom1[1], atom1[2], atom1[3]])
            s = atom1[4][12:16]
            s = re.sub("\d+", "", s)
            if s != 'HB':
                try:
                    mass = pt.masses[atom1[0]]
                except KeyError:
                    mass = 0
                atom2 = atom(atom1[0], mass, array, ring.group([],0),atom1[4])
                #print(atom2.ID, atom2.pdb[6:11])
                #atom2.ID = int(atom1[4][6:11])
#                if atom2.symbol == 'N' and atom2.pdb[17:20] == 'PRO':
#                    atom2.max_valence = 4
    #            print('{} atom symbol valence{}\n'.format(atom2.symbol, atom2.max_valence))
                atom_list.append(atom2)
                
        return n, atom_list

def read_in_molecule_data_pdb(file):
    """read in the atoms and molecular data for a molecule"""
    n, atoms = read_pdbfile(file)
    Molecule = molecule(atoms)
    return Molecule, n


#Step 1. In python, analyse topology of system to be modelled, returning:
#
#Atomic positions (trivial), polarizabilities (not so trivial), partial charges 
#   (hardest, I will have to look into how to best do this)
#
#-Bonded atom pairs -> bonds, bond force constants
#
#-Bonded atom triplets -> angles, bond angle parameters (and symmetrized linear combinations)
#
#-Bonded atom quartets -> dihedral angles, dihedral angle parameters (and symmetrized linear combinations)
#
#Note: See attached papers for definition of natural internal coordinates, 
#a procedure for removing redundancies amongst coordinates, using only local connectivity information. 
#Also, a simple algorithm for determining which atoms are in topologically similar or different environments.

 
# Class
#Also, I particularly recommend setting up an Atom class with properties 
#symbol, mass, position, radius, charge. And perhaps a molecule class with properties 
#related to the molecular connectivity – bonds, angles, dihedrals, and of course each atom can be ‘stored’ 
#in the molecule


'==============================================================================='

  
def atomic_brakedown(molecule):
    """takes atoms from molecule and then returns a list directery with the 
    number of each atom type present"""
    atom_types = {}
    for atom in molecule.atoms:
        atom_types[(atom.symbol,atom.ringfinder)] = 1 + atom_types.get((atom.symbol,atom.ringfinder), 0)
    return atom_types
    


def molecule_onclaves(molecule):
    """Function finds onclaves in a molecule by generating groups based upon 
    atomic connectivity and then returns those groups """
    atom_set = []
    all_atoms = set()
    
        

    while not len(all_atoms) == len(molecule.atoms):
        for atom5 in molecule.atoms:
            if not atom5 in all_atoms:
                start = atom5
            else:
                continue
            atoms = set([start])
            new_atoms = [start]
            n = 0
            while len(atoms) > n:
                n = len(atoms)
                for atom1 in new_atoms:
                    for atom2 in ring.get_atoms(atom1.connectivity):
                        atoms.add(atom2)
                new_atoms = [atom4 for atom4 in atoms]

            atom_set.append(atoms)
            all_atoms = atom_count(atom_set)

    return atom_set
    
def atom_count(atom_set):
    """Inputs a list of sets containing atom is and outputting A set with
    all the atoms from the sets into it"""
    all_atoms = set()
    for set1 in atom_set:
        
        all_atoms.update(list(set1))
    
    return all_atoms
    
    

# =========================================================================
    
# BOND ORDER CODE 

def bond_processing(molecule):
    """setup bond valence and options and sparevalence"""
    sparevalence = []
    for atom in molecule.atoms:
        atom.Valence = len(atom.connectivity)
        atom.SpareValence = atom.max_valence - atom.Valence
        if atom.SpareValence > 0:
            sparevalence.append(atom)
#        print('SpareValence {} = max_valence {}  Valence {}'.format(atom.SpareValence, atom.max_valence, atom.Valence))
    
    for atom in sparevalence:
        optionsatoms = []
        for bond in atom.connectivity:
            for atom1 in bond.atoms:
                if atom1 != atom and atom1.SpareValence > 0:
                    optionsatoms.append(atom1.SpareValence) 
#        print('optionsatoms {}'.format(optionsatoms))
        atom.options = len(optionsatoms)
#        print('atom.iD {} = {}'.format(atom1.ID,atom1.options))
    
#    for atom1 in molecule.atoms:
#        print('atom.iD {} = {}'.format(atom1.ID,atom1.options))
    
    # after finish graphs turn back on
    
    unsaturated_bond(molecule)
    
    bond_rings(molecule)
#    bond_rings(molecule)
    processing_charges(molecule)
    unsaturated_chains(molecule)
#    processing_charges(molecule)
    
    carbonylate_ring(molecule)
    
    return None

## Increment bond order for isolated unsaturated bonds
#for bond in molecule.Bonds:
#	if bond.Atoms[0].Options == bond.Atoms[1].Options == 1:
#  if bond.Atoms[0].SpareValence == bonds.Atoms[1].SpareValence:
#    bond.Order += bond.Atoms[0].SpareValence
#    for atom in bond.Atoms:
#	atom.SpareValence = 0
#      atom.FoundValence += bond.Atoms[0].SpareValence
#else: 
#	SV = [atom.SpareValence for atom in bond.Atoms]
#	bond.Order += min(SV)
#	for atom in bond.Atoms:
#		atom.SpareValence = -min(SV)
#			atom.FoundValence += min(SV)	
def unsaturated_bond(molecule):
    """finds isolated unsaturated bonds"""
    
    
    
    for bond in molecule.bonds:
#        print('bond.atoms[0].options {} == bond.atoms[1].options {}'.format(bond.atoms[0].options,bond.atoms[1].options))
        
        if bond.atoms[0].options == bond.atoms[1].options == 1:
#            print('if bond.atoms[0].options == bond.atoms[1].options == 1:')
            if bond.atoms[0].SpareValence == bond.atoms[1].SpareValence:
#              print('if bond.atoms[0].SpareValence == bond.atoms[1].SpareValence:')
              bond.bondorder += bond.atoms[0].SpareValence
#              if bond.atoms[1].SpareValence > 1:
#                  print('bondorder {} = 1 + SpareValence {}'.format(bond.bondorder,bond.atoms[0].SpareValence))
              newsparevalence = bond.atoms[0].SpareValence  
              for atom in bond.atoms:
                  atom.SpareValence = 0
                  atom.Valence += newsparevalence
            else: 
              SV = [atom.SpareValence for atom in bond.atoms]
              bond.bondorder += min(SV)
              for atom in bond.atoms:
        		     atom.SpareValence = -min(SV)
        		     atom.Valence += min(SV)
    for atom in molecule.atoms:
        for bond in atom.connectivity:
            for bond2 in  molecule.bonds:
                bond.atoms.sort(key=lambda atom: atom.ID)
                bond2.atoms.sort(key=lambda atom: atom.ID)
#                print('{} == {}'.format(bond.atoms, bond2.atoms))
                if bond.atoms == bond2.atoms:
                    atom.connectivity.remove(bond)
                    atom.connectivity.append(bond2)
    
    return None



## Adjust bond orders for conjugated systems, starting with rings
#for atom in molecule.Atoms:
#  atom.SpareValence = atom.SpareValence/atom.Options if atom.Options != 0
#
#for ring in molecule.Rings:
#  for atom in ring:
#    for bond in atom.Bonds:
#      bond.Order += atom.SpareValence
#      for bound_atom in bond.Atoms:
#        bound_atom.SpareValence = 0
#        bound_atom.FoundValence += atom.SpareValence
#
#

def bond_rings(molecule):
    """sorts  out the bond order of ring bonds"""
    for atom in molecule.atoms:
        if atom.options != 0:
            atom.SpareValence = atom.SpareValence/atom.options 
        
    for ring1 in molecule.rings:
        sparevalencenum = []
        for atom in ring1:
            if atom.symbol != 'Q':
                sparevalencenum.append(atom.SpareValence)
#        print('sparevalencenum {}'.format(sparevalencenum))
        sparevalencenum = list(set(sparevalencenum))
        
#        print('len(sparevalencenum) == {} and sparevalencenum {} == 0.5'.format(len(sparevalencenum),sparevalencenum[0]))
        if len(sparevalencenum) == 1 and sparevalencenum[0] == 0.5:
            for atom in ring1:
                
                for bond in atom.connectivity:
                    if all(i in ring1 for i in bond.atoms):
                        bond.bondorder += 0.5
                        for bound_atom in bond.atoms:
                            bound_atom.SpareValence = 0
                            bound_atom.Valence += atom.SpareValence
        elif len(sparevalencenum) == 2 and all(i in sparevalencenum for i in [1/3,0.5]):
            for atom in ring1:
                if atom.SpareValence == 1/3:
                    atom.SpareValence = 0.5
#                    print('atom id {}'.format(atom.ID))
                    
            for atom in ring1:
                
                for bond in atom.connectivity:
                    if all(i in ring1 for i in bond.atoms):
                        bond.bondorder += 0.5
                        for bound_atom in bond.atoms:
                            bound_atom.SpareValence = 0
                            bound_atom.Valence += atom.SpareValence

## Then moving on to alternating unsaturated chains
#WholeSpareValences = Count(atoms in molecule with SpareValence = 1)
#while WholeSpareValences > 0:
#  for bond in molecule.Bonds: 
#    for i,atom in enumerate(bond.Atoms):
#	if atom.SpareValence = 1:
#        bond.Order += 1
#        other_atom = bond.Atoms[i-1]
#        other_atom.FoundValence = len(other_atom.Connectivity) + 1
#        other_atom.SpareValence = other_atom.ExpectedValence - 
#                                  other_atom.FoundValence
#	  for bond in other_atom.Bonds:
#     if bond.Atom[1].symbol != H: # not other_atom
#bond.Atom[1].SpareValence = 1 # sets next atoms in the chain          to 1 to the next bond in the conjugated group can be processed in the next iteration.
#  WholeSpareValences = Count(atoms in molecule with SpareValence = 1)
#



#def make_group_unsaturated_chains(molecule):
#    """ """
#    group_set = set()
#    for atom in molecule.atoms:
#        
#        if atom.SpareValence > 0:
#            
#            atoms = set([atom])
#            current_atoms = [atom]
#            n = 0
#            while len(atoms) > n:
#                n = len(atoms)
#                
#                for atom2 in current_atoms:
#                    for bond in atom2.connectivity:
#                        for atom3 in bond.atoms:
#                            if atom2 != atom3 and atom3.SpareValence > 0:
#                                atoms.add(atom3)
#                            
#                current_atoms = [atom4 for atom4 in atoms]
#            
#                
#            atoms = [atom5 for atom5 in atoms]
#            check = [] 
#            for atom4 in atoms:
#                check.append(atom4.SpareValence)
#            if all(i in [1,0.5] for i in check):
#                group_set.add(group(atoms, 0.5))
#            
#    for group1 in group_set:
#        (group1.atoms).sort(key = lambda x: x.__hash__())
#    
#            
#    group_list = [group1 for group1 in group_set]
#    group_list = list(set(group_list))
#    bonds_lists = []
#    for group1 in group_list:
#        bonds = []
#        for atom2 in group1.atoms:
#            for bond1 in atom2.connectivity:
#                if all(i in group1.atoms for i in bond1.atoms):
#                    bonds.append(bond1)
#        bonds = list(set(bonds))
#        for bond2 in bonds:
#            for atom3 in bond2.atoms:
#                for bond3 in atom3.connectivity:
#                    if bond3 != bond2 and bond3 in bonds:
#                        bond2.connections.append(bond3)
#            bond2.connections = list(set(bond2.connections)) 
#        bonds_lists.append(bonds)
#    
#    
#    return bonds_lists #list(set(group_list))



def unsaturated_chains(molecule):
    """the finds the unsaturated chains bond orders"""

       
    WholeSpareValences = [0]
    n = 0
    num = float('inf')
    while len(WholeSpareValences) > 0:
        for atom in molecule.atoms:
            optionsatoms = []
            for bond in atom.connectivity:
                for atom1 in bond.atoms:
                    if atom1 != atom and atom1.SpareValence > 0:
                        optionsatoms.append(atom1.SpareValence) 
            atom.options = len(optionsatoms)
            if atom.SpareValence == 0.5 and atom.options == 1:
                atom.SpareValence = 1 
        
        for bond in molecule.bonds:
    #        sp = []
            
            k = 0
            bond.SP_list = []
            for atom1 in bond.atoms:
                
                bond.SP_list.append(atom1.SpareValence)
                k += 1
            
            bond.SP_list.sort()
#            if [0.5,1] == bond.SP_list or [1,1] == bond.SP_list: #all(i in [1,0.5] for i in bond.SP_list):
#                print('[0.5, 1] == {}, {}'.format(bond.SP_list,bond))
            
            if [0.5,1] == bond.SP_list or [1,1] == bond.SP_list:
                bond.bondorder += 1
                for atom1 in bond.atoms:
                        atom1.SpareValence = 0
                        atom1.Valence += 1 
                        atom.options = 0
        WholeSpareValences = [atom.SpareValence for atom in molecule.atoms if atom.SpareValence > 0]
#        print('WholeSpareValences {}'.format(len(WholeSpareValences)))
        if len(WholeSpareValences) < num:
            num = len(WholeSpareValences)
        else: 
            n += 1
            if n > 20:
                print("failed")
                break
    
    return None       

    
def printinf(molecule,num):
    for bond in molecule.bonds:
            str1 = ''
            ids = []
            for atom2 in bond.atoms:
                str1 += '{} '.format(atom2.ID)
                ids.append(atom2.ID)
            ids.sort()
#            if all(i in [[7,8],[9,10],[11,12],[13,14],[15,16],[17,18],[19,20]] for i in [ids]):
#                print('{} bond.atoms.ID [{}], bondorder {}'.format(num,str1,bond.bondorder))      
    return None


# try:
#    
##        print('optionsatoms {}'.format(optionsatoms))
#        
##        if atom.ID == 21 or atom.ID == 5:
##            print('atom{}.options = {}'.format(atom.ID,atom.options))
##    for atom in molecule.atoms:
##        
##        atom.option = len([atom for atom in atom.connectivity if atom.SpareValence > 0])
#    
#    for atom in molecule.atoms:
#        if atom.SpareValence == 0.5 and atom.options == 1:
#            atom.SpareValence = 1
#    
#    WholeSpareValencesatoms = [atom for atom in molecule.atoms if atom.SpareValence == 1]
#    WholeSpareValences = [atom.SpareValence for atom in molecule.atoms if atom.SpareValence == 1]
#    while  all(i in [[0,0],[1,1]] for i in [WholeSpareValences]): #len(WholeSpareValences) > 0:
#        
#        
#        for atom in WholeSpareValencesatoms:
##            print('1atom{}.SpareValence = {}'.format(atom.ID,atom.SpareValence))
#            for bond in atom.connectivity:
#                str1 = ''
#                for atom2 in bond.atoms:
#                    str1 += '{} '.format(atom2.SpareValence)
##                print('all(atom.SpareValence > 0 , bond.atoms = [{}]'.format(str1))
#                if all(atom.SpareValence > 0 for atom in bond.atoms):
##                    print('2atom{}.SpareValence = {}'.format(atom.ID,atom.SpareValence))
#                    str1 = ''
#                    ids = []
#                    for atom2 in bond.atoms:
#                        str1 += '{} '.format(atom2.ID)
#                        ids.append(atom2.ID)
#                    ids.sort()
#                    if all(i in [[7,8],[9,10],[11,12],[13,14],[15,16],[17,18],[19,20]] for i in [ids]):
#                        print('1bond.atoms.ID [{}], bondorder {}'.format(str1,bond.bondorder))
#                    
#                    
#                    bond.bondorder = 2
#                    str1 = ''
#                    ids = []
#                    for atom2 in bond.atoms:
#                        str1 += '{} '.format(atom2.ID)
#                        ids.append(atom2.ID)
#                    ids.sort()
#                    if all(i in [[7,8],[9,10],[11,12],[13,14],[15,16],[17,18],[19,20]] for i in [ids]):
#                        print('2bond.atoms.ID [{}], bondorder {}'.format(str1,bond.bondorder))
#                    for atom1 in bond.atoms:
##                        print('2atom{}.SpareValence = {}'.format(atom1.ID,atom1.SpareValence))
#                        atom1.SpareValence = 0
#                        atom1.Valence += 1
#                        
##                    for bond1 in atom.connectivity:
##                        for atom2 in bond1.atoms:
##                            if atom2 =
#                        
#      
#                
#        for atom in molecule.atoms:
#            optionsatoms = []
##            for bond in atom.connectivity:
##                str1 = ''
##                for atom2 in bond.atoms:
##                    str1 += '{} '.format(atom2.ID)
##                str1 = ''
##                ids = []
##                for atom2 in bond.atoms:
##                    str1 += '{} '.format(atom2.ID)
##                    ids.append(atom2.ID)
##                ids.sort()    
##                
##                if all(i in [[7,8],[9,10],[11,12],[13,14],[15,16],[17,18],[19,20]] for i in [ids]):
##                    print('2bond.atoms.ID [{}], bondorder {}'.format(str1,bond.bondorder))
#                
#            for atom1 in bond.atoms:
#                if atom1 != atom and atom1.SpareValence > 0:
#                    optionsatoms.append(atom1.SpareValence) 
#                atom.options = len(optionsatoms)
#                
#           for atom in molecule.atoms:
#                if atom.SpareValence == 0.5 and atom.options == 1:
#                    atom.SpareValence = 1     
#            
##                str1 = ''
##                for atom2 in bond.atoms:
##                    str1 += '{} '.format(atom2.ID)
#                
##                str1 = ''
##                ids = []
##                for atom2 in bond.atoms:
##                    str1 += '{} '.format(atom2.ID)
##                    ids.append(atom2.ID)
##                ids.sort() 
##                if all(i in [[7,8],[9,10],[11,12],[13,14],[15,16],[17,18],[19,20]] for i in [ids]):
##                    print('3bond.atoms.ID [{}], bondorder {}'.format(str1,bond.bondorder))
#    #        print('optionsatoms {}'.format(optionsatoms))
#            
#        
#        
#        
#        
#        
#        WholeSpareValencesatoms = [atom for atom in molecule.atoms if atom.SpareValence == 1]
#        
#        WholeSpareValences = [atom.SpareValence for atom in molecule.atoms if atom.SpareValence == 1]
##        for bond in molecule.bonds: 
##            for i,atom in enumerate(bond.atoms):
##                if atom.SpareValence == 1:
##                    bond.bondorder += 1
##                    other_atom = bond.atoms[i-1]
##                    other_atom.Valence = len(other_atom.connectivity) + 1
##                    other_atom.SpareValence = other_atom.max_valence - other_atom.Valence
##                    for bond in other_atom.connectivity:
##                        if bond.atoms[1].symbol != 'H': # not other_atom
##                             bond.atoms[1].SpareValence = 1 # sets next atoms in the chain          to 1 to the next bond in the conjugated group can be processed in the next iteration.
#        
##        WholeSpareValences = len([atom for atom in molecule.atoms if atom.SpareValence == 1])
##        print('WholeSpareValences {}'.format(WholeSpareValences))
##        for atom in molecule.atoms:
##            if atom.SpareValence == 0.5 and atom.option == 1:
##                atom.SpareValence = 1
#    
#    for bond in molecule.bonds:
#        str1 = ''
#        ids = []
#        for atom2 in bond.atoms:
#            str1 += '{} '.format(atom2.ID)
#            ids.append(atom2.ID)
#        ids.sort()
#        if all(i in [[7,8],[9,10],[11,12],[13,14],[15,16],[17,18],[19,20]] for i in [ids]):
#            print('3bond.atoms.ID [{}], bondorder {}'.format(str1,bond.bondorder))
#    
#    WholeSpareValences = len([atom for atom in molecule.atoms if atom.SpareValence == 1])
##    except AttributeError:
##        print('molecule {}'.format())













## Derive formal charges for each atom
#for atom in molecule.Atoms:
#  atom.SpareValence = atom.ExpectedValence - atom.Bond.Order
#  atom.Charge = -atom.SpareValence
#if atom.symbol == 'N': atom.Charge += 1
## Bond orders for carboxylate rings
#
#
    
def processing_charges(molecule):   
    """charge find for each atoms"""
    for atom in molecule.atoms:
        if atom.options == 1:
            atom.charge = -atom.SpareValence
            if atom.symbol == 'N' and atom.max_valence == 4: 
                atom.charge += 1
    
    
    return None
# Bond orders for carboxylate rings
#for ring in molecule.Rings:
#	if len(ring) == 4:
#		atom_symbol = [atom.symbol for atom in ring]
#	if atom_symbol.sort() == [‘C’,’O’,’O’,’Q’]
#		for atom in ring:
#			if atom.symbol == ‘Q’:
#				for bond in atom.Bonds:
#					if bond.Atoms[0] == bond.Atoms[0] == ‘O’:
#						bond.Bondorder = 0
#					else:
#					      bond.Bondorder = 1.5


def carbonylate_ring(molecule):
    """processing the carbonylate rings to curretly find the bond type and 
    charges on the Oxygen atoms"""
    rings4 = []
    for ring1 in molecule.rings:
    	if len(ring1) == 4:
            atom_symbol = [atom.symbol for atom in ring1]
            atom_symbol.sort()
            coo = ['C','O','O','Q']
            coo.sort()
            cno = ['C','N','O','Q']
            cno.sort()
            cnn = ['C','N','N','Q']
            cnn.sort()
#            rings4.append(atom_symbol)
            
        		  #collections.Counter(atom_symbol.sort()) == collections.Counter(['C','O','O','Q'].sort())
#            print('atom_symbol {} == {}, Bool {}'.format(atom_symbol.sort(),['C','O','O','Q'].sort(),atom_symbol.sort() == ['C','O','O','Q'].sort()))
            if atom_symbol == coo:
#                print('1 atom_symbol {}'.format(atom_symbol))
                for atom in ring1:
                    
                    if atom.symbol != 'Q':
                        for bond in atom.connectivity:
                            if bond.atoms[0].symbol == bond.atoms[1].symbol == 'O':
                                bond.bondorder = 0
                                for atom1 in  bond.atoms:
                                    atom1.charge = -0.5
                                molecule.bonds.remove(bond)
                                bond.atoms[0].connectivity.remove(bond)
                                bond.atoms[1].connectivity.remove(bond)
                            elif all(i in ring1 for i in bond.atoms) and all(i in [bond.atoms[0].symbol,bond.atoms[1].symbol] for i in ['C','O']):
                                bond.bondorder = 1.5
#                                print('1 bond {}, bondorder {}'.format(bond,bond.bondorder))
                                
                                
#            for atom in ring1:
#                for bond in atom.connectivity:
#                    if all(i in ring1 for i in bond.atoms) and all(i in [bond.atoms[0].symbol,bond.atoms[1].symbol] for i in ['C','O']):
#                        print('2 bond {}, bondorder {}'.format(bond,bond.bondorder))
            elif atom_symbol == cno:
#                print('2 atom_symbol {}'.format(atom_symbol))
                for atom in ring1:
                    
                    if atom.symbol != 'Q':
                        for bond in atom.connectivity:
                            if all(i in [bond.atoms[0].symbol,bond.atoms[1].symbol] for i in ['N','O']):
                                bond.bondorder = 0
                                
                                molecule.bonds.remove(bond)
                                bond.atoms[0].connectivity.remove(bond)
                                bond.atoms[1].connectivity.remove(bond)
                            elif all(i in [bond.atoms[0].symbol,bond.atoms[1].symbol] for i in ['C','O']):
                                bond.bondorder = 2
                                
                            else:
                                bond.bondorder = 1
#            for atom in ring1:
#                for bond in atom.connectivity:
#                    if all(i in ring1 for i in bond.atoms) and all(i in [bond.atoms[0].symbol,bond.atoms[1].symbol] for i in ['C','O']):
#                        print('3 bond {}, bondorder {}'.format(bond,bond.bondorder))
#                                
            elif atom_symbol == cnn:
#                print('3 atom_symbol {}'.format(atom_symbol))
                h_atoms_symbol = []
                
                for atom in ring1:
                    n_atom_con_sybol =[]
                    if atom.symbol == 'N':
                        
                        for bond1 in atom.connectivity:
                            for atom2 in bond1.atoms:
                                if atom2.symbol == 'H':
                                    h_atoms_symbol.append(atom2.symbol)
                                
                                    
                if len(h_atoms_symbol) == 4:
                    for atom in ring1:
                        n_atom_con_sybol =[]
                        if atom.symbol == 'N':
                            
                            for bond1 in atom.connectivity:
                                for atom2 in bond1.atoms:
                                     if atom2.symbol == 'H':
                                         n_atom_con_sybol.append(atom2.symbol)
#                            print('n_atom_con_sybol {}'.format(n_atom_con_sybol))
                            if len(n_atom_con_sybol) == 2:
                                atom.charge = 0.5
#                                print('atom {}, {}'.format(atom.ID,atom.charge))
                            
                    for atom in ring1:
#                        if atom.symbol != 'C':
#                            atom.charge = 0
                        if atom.symbol != 'Q':
                            for bond in atom.connectivity:
                                if  bond.atoms[0].symbol == bond.atoms[1].symbol == 'N':
                                    bond.bondorder = 0
                                    
                                    molecule.bonds.remove(bond)
                                    bond.atoms[0].connectivity.remove(bond)
                                    bond.atoms[1].connectivity.remove(bond)
                                elif all(i in ring1 for i in bond.atoms) and all(i in [bond.atoms[0].symbol,bond.atoms[1].symbol] for i in ['C','N']):
                                    bond.bondorder = 1.5
#    for ring2 in molecule.rings:
#        for atom in ring2:
#            for bond in atom.connectivity:
#                if all(i in ring2 for i in bond.atoms) and all(i in [bond.atoms[0].symbol,bond.atoms[1].symbol] for i in ['C','O']):
#                    print('3 bond {}, bondorder {}'.format(bond,bond.bondorder))
#            if atom.ID == 78 or atom.ID == 79:
#                print('3atom {} , {}'.format(atom.ID,atom.charge))
                        




