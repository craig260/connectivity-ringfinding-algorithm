# -*- coding: utf-8 -*-
"""
Created on Mon Jan 22 14:12:26 2018

ring finding algothem v 2

@author: Craig Melton
"""

from mpl_toolkits.mplot3d import Axes3D
import matplotlib.pyplot as plt

import force_field_code_final as FF

import numpy as np




atomtype = {'C': 'b', 'O': 'r', 'N': 'g', 'H': 'w', 'S': 'y'}
grouptype = {2: 'b', 1: 'r', 3: 'g', 0: 'w', 4: 'y'}

def group_graph(molecule, groups_list, colortype):
    """construting connectivit graph"""
    for group_list in groups_list:
        groups = []
        
        for group in group_list:
            try:
                if group.valence > 0:
                    groups.append(group)
            except AttributeError:
                continue
#                for ring1 in groups_list:
#                    print(ring1)
#                    rg = ring.group(ring1,20)
                    #print(rg)
#                    groups.append(rg)
        plt.style.use('dark_background')
        fig = plt.figure()
        ax = fig.add_subplot(111, projection='3d')
        
        
        for group in groups:
            
            for atom in group.atoms:
                if colortype == 'group':
                    try:
                        colour = grouptype[atom.ringfinder]
                    except KeyError:
                        colour = 'm'
                if colortype == 'atom':
                    try:
                        colour = atomtype[atom.symbol]
                    except KeyError:
                        colour = 'm'
                x = atom.position[0]
                y = atom.position[1]
                z = atom.position[2]
            
                ax.scatter(x, y, z, c=colour ,marker='o')
                ax.text(x, y, z, atom.symbol)
            
        
        
        for group in groups:
            
            for atom2 in group.atoms:
                for bond in atom2.connectivity:
                    try:
                        clr = 'w-'
                        x1 = atom2.position[0]
                        y1 = atom2.position[1]
                        z1 = atom2.position[2]
                        x2 = bond[1].position[0]
                        y2 = bond[1].position[1]
                        z2 = bond[1].position[2]
                        if bond[0] >= 2.0:
                            clr = 'r-'
                        ax.plot([x1,x2], [y1,y2], [z1,z2], clr)
    
                    except AttributeError:
                        continue
        ax.set_xlabel('x')
        ax.set_ylabel('y')
        ax.set_zlabel('z')
        
        
        plt.show()
    return None

def count_non_zeroatoms(molecule):
    """Counts the number of atoms that have a nonzero ring find a number and 
    returns that count"""
    non_zero = []
    for atom1 in molecule.atoms:
        if atom1.ringfinder > 0:
            non_zero.append(atom1)
    return len(non_zero)


def atom_problem(molecule):
    """finds and prints atoms with atomic postions that are an order magnitude
    high"""
    i = 0
    for atom in molecule.atoms:
        for num in atom.position:
            if num > 100:
                print(atom)
                i = 1
    r = False
    if i == 1:
        r = True
    return r

def get_atoms(bonds):
    """get the atoms from the connectivity list of atoms which have the form
    [(bond_length,atom), ...]"""
    atoms1 = []
    for atom in bonds:
#        print(atom)
#        print(type(atom))
        try:
            atoms1.append(atom[1])
        except TypeError:
            atoms1.append(atom.atoms[1])
    atoms1 = set(atoms1)
    return list(atoms1)



class group:
    """contains list of atoms that are connect and have the same number of 
    nonhydrogen bonds along with other useful properieties of these groups"""
    def __init__(self, atoms, valence):
        """properties of these groups"""
        atoms.sort(key = lambda x: x.__hash__())
        self.atoms = atoms
        self.size = len(atoms)
        if len(atoms) > 0:
            self.valence =  valence     #atoms[0].ringfinder
            self.expectedconnections = 2*(self.valence-1)+(self.size-2)*(self.valence-2)
            self.position = np.array([0.0,0.0,0.0])
            for atom in (self.atoms):
                self.position +=  atom.position
            self.position = (self.position/self.size)
            
        else:
            self.valence = 0
            self.expectedconnections = 0
        self.exteronalconnections = 0
        
        self.atomicconnections = []   # construct with sets so atoms are unqie
        self.groupconnection = set()
        if self.expectedconnections == self.exteronalconnections:
            self.confirmedchain = True
        else:
            self.confirmedchain = False
        self.groupnumber = 0
        
        
            
    def __repr__(self):
        """evaluate atom class returns string"""

        return '({}, {}, {})'.format(self.valence, self.size, self.__hash__())
        
    def __str__(self):
        """what the atom class returns when printing"""

        return '({}, {}, {})'.format(self.valence, self.size, self.__hash__())

    def __key(self):
        """ """
        return (self.valence, self.size, tuple(self.atoms))
    def __eq__(self, other):
        """ """
        if  isinstance(other, self.__class__) and other.__key() == self.__key():
            return True
        else:
            return False
    def __hash__(self):
        """hash number for atom"""
        return hash(self.__key())           
        
    def internal_bonds(self):
        """Takes the group and returns a truple of the atom with a list of
        other atoms bonded to this atom in the same group"""
        atoms = self.atoms
        internal_bonds = []
        for atom in atoms:
            bonds = get_atoms(atom.connectivity)
            bondatoms = set()
            for atom2 in bonds:
                if atom2 in atoms:
                    bondatoms.add(atom2)
            internal_bonds.append((atom,bondatoms))
        return internal_bonds

        
def connections_num(molecules):
        """take connections data from the atoms and find the number of non H 
        atoms connected to each atoms then changes non_H_connections porpertie
        of the atom"""
        atoms =  molecules.atoms
        for atom in atoms:
            if atom.symbol == ('H' or  'Q'):
                continue
            if atom.symbol == 'S':
                for atom2 in get_atoms(atom.connectivity):
                    if not atom2.symbol in ['H' , 'Q' , 'S']:
                        atom.ringfinder += 1
                

            else:
                for atom2 in get_atoms(atom.connectivity):
                    if not atom2.symbol in ['H' , 'Q']:
                        atom.ringfinder += 1
                
def n_zeroatoms(atoms, n=0):
    """takes a list of atom and returns a count of the atoms with non zero 
    ringfinder number"""
    n_atoms = []
    for atom in atoms:
        if atom.ringfinder == n:
            n_atoms.append(atom)
    return n_atoms



def chain_end(molecule):
    """fuction that find the chain  the atom that have 1 as the value for 
    ringfinder"""
    while non_zeroatoms(molecule.atoms, n=1) > 0:
        chainends = n_zeroatoms(molecule.atoms, n=1)
        
        atom_connection = []
        for atom in chainends:
            if atom.ringfinder == 1:
                atom.ringfinder -= 1
                
            
            atom_connection.append(get_atoms(atom.connectivity))
            current_atoms = []
            for connections in atom_connection:
                for atom2 in connections:
                    if atom2.symbol != ('H' or 'Q'):
                        current_atoms.append(atom2)
        for atom in current_atoms:
            if atom.ringfinder > 0:
                atom.ringfinder -= 1
                
            
                       


# intergrated identifetion and remove of chain ends

#def chain_remove(molecule):
#    """find the chain end in a molecule then indivily constructs sets and then 
#    remove these atoms from the molecule by setting ringfinder to zero"""
#    atoms = molecule.atoms
#    chainends = []
#    for atom in atoms:
#        if atom.ringfinder == 1:
#            chainends.append(atom)
#    
#    for atom in chainends:
#        chain = set()
#        chain.add(atom)
#        newatom = atom
#        i = 0           # i is on off switch for if new atom as been add to chain
#        while i < 1:    # loops the addtion of atoms until no new atoms are add to chain
#            for atom2 in get_atoms(newatom.connectivity):
#                if atom2.ringfinder == 2:
#                    i = 1
#                    if not (atom2 in chain):
#                        newatom = atom2
#                        i = 0
#                    chain.add(atom2) 
#        
#        for atom in chain:
#            atom.ringfinder = 0
            
# ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# find the atoms in these groups that are directly connected to other groups

def group_linkers(group_list):
    """takes the lists of atomic groups and finds the atoms attaching these groups
    to other groups"""
    for group in group_list:
        #atoms = group.atoms
        connections = set()
        for atom in group.atoms:
            
            for atom2  in get_atoms(atom.connectivity):
                if atom2 in group.atoms:
                    continue
                if atom2.ringfinder > 0:
                    connections.add(atom2)
        (group.atomicconnections) = [atom for atom in connections]
        group.exteronalconnections = len(group.atomicconnections)
        if group.expectedconnections == group.exteronalconnections:
            group.confirmedchain = True
    for group in group_list:
        atoms2 = group.atomicconnections
        
        for atom in atoms2:
            for group2 in group_list:
                if atom in group2.atoms:
                    (group.groupconnection).add(group2)
    return None
                    
                    
def homogesis_atom_connections(molecule):
    """takes the atom list and finds the groups of atom that have the same 
    number of nonhydrogen bonded atoms that are connected and puts them into a
    list of groups"""
    
   
    group_set = set()
    for atom in molecule.atoms:
        
        if atom.ringfinder > 0:
            
            atoms = set([atom])
            current_atoms = [atom]
            n = 0
            while len(atoms) > n:
                n = len(atoms)
                
                for atom2 in current_atoms:
                    for atom3 in get_atoms(atom2.connectivity):
                        if atom.ringfinder == atom3.ringfinder:
                            atoms.add(atom3)
                            
                current_atoms = [atom4 for atom4 in atoms]
            
                
            atoms = [atom5 for atom5 in atoms]
            
            
            
            group_set.add(group(atoms, atom.ringfinder))
            
    for group1 in group_set:
        (group1.atoms).sort(key = lambda x: x.__hash__())
    
            
    group_list = [group1 for group1 in group_set]
    m = 0
    for groups in group_list:
        groups.groupnumber = m
        for atom6 in groups.atoms:
            atom6.groupnumber = m
        m += 1
                
    return list(set(group_list)), len(list(set(group_list)))

def atom_groups(group_list):
    """add groups to a list in atoms structure"""
    for group1 in group_list:
        for atom1 in group1.atoms:
            (atom1.groups) = group1
            



def simply_ring(molecule):
    """finds simple rings madeup of two groups and with 2 and 3 as atomic 
    ringfinder and the 3 as size 1 numbers to create new chain ends""" 
    group_of_v2 = [group for group in molecule.group_list if group.valence == 2]
    
    new_ring = []
    for group1 in group_of_v2:
        i = 0
        if len(group1.atomicconnections) == 1:
            
            for atom1 in group1.atomicconnections:
                for atom2 in atom1.connectivity: 
                    if atom2[1] in group1.atoms:
                        i += 1
#            if group1.size == 64:
#                    group_graph(molecule,[[group1]+list(group1.groupconnection)], 'group')
#                
#            print('i = {}'.format(i))
#            print('i = {}'.format(i))
            if i  > 1:
                if i == 1:
                    print('i = {} error'.format(i))
#                print((group1,group1.groupconnection))
#                if group1.size == 64:
#                    print('i = {}'.format(i))
        #                group_graph(molecule,[list(group1.groupconnection)], 'atom')
        #                molecule.intersting_groups.append()
                
    #                    print(1)
                for atom in group1.atoms:
                    atom.ringfinder = 0
                group1.atomicconnections[0].ringfinder -= 2
                        
                        
                (group1.atoms).append(group1.atomicconnections[0])
                rings = group1.atoms
                    
                new_ring.append(rings)
                
            else:
                continue
        
    ringprocesser(new_ring, molecule)
    
    

def ringprocesser(new_ring, molecule):
    """process the a ring that as been identitifyied and defines a q atom with
    approperate valence and position in space then add them to the ring and 
    atom list"""
    for ring in new_ring:
        
        size = len(ring)
        if size > 0:
            q_position = [0,0,0]
            if size < 3 or size == 8:
                 molecule.ring_errors.append(ring)
                 
            elif size == 3:
                molecule.carbonyl_group.append(ring)
                 
            for atom in ring:
                
                q_position[0] += (atom.position)[0]
                q_position[1] += (atom.position)[1]
                q_position[2] += (atom.position)[2]
            
            pdbline = '{:<6}{:>5}{:<4}{:<1}{:<3} {:<1}{:<4}{:<1}   {:>8}{:>8}{:>8}{:<6}{:6}           {:>2}{:>2}'.format('ATOM',0,'Q','', '' ,'','',q_position[0] ,q_position[1], q_position[2],'','','','Q','')        
            q_position = np.array([x/size for x in q_position])
            
            q = FF.atom('Q', 0 , q_position, group([],-1),pdbline)
            q.max_valence = size
            
            for atom2 in ring:
                    length = np.linalg.norm(q.position - atom2.position)
                    (q.connectivity).append((length,atom2)) #FF.Bond
                    
            
            molecule.atoms.append(q)
            
            ring.append(q)
            
            if len(ring) > 3:
                molecule.rings.append(ring)

def carbonyl_groups(molecule):
    """finds the carbonyl groups"""
    for C_atoms in molecule.atoms:
        carbonyl_atoms = [C_atoms]
        O_num = 0
        for O_atoms in C_atoms.connectivity:
            if O_atoms[1].symbol == 'O':
                H_num = 0
                for H_atoms in O_atoms[1].connectivity:
                    if H_atoms[1].symbol == 'H':
                         H_num += 1
                if H_num == 0:
                    O_num += 1
                    carbonyl_atoms.append(O_atoms[1])
        if O_num == 2:
            molecule.carbonyl_group.append(carbonyl_atoms)

def non_zeroatoms(atoms, n=0):
    """takes a list of atom and returns a count of the atoms with non zero 
    ringfinder number"""
    non_zero_atoms = []
    for atom in atoms:
        if atom.ringfinder == n:
            non_zero_atoms.append(atom)
    return len(non_zero_atoms)

def group_one( molecule, n):
    """finding rings make of two group"""
    group_of_v2 = [group for group in molecule.group_list if group.valence == 2]
    ring = []
    
    for groups in group_of_v2:
        i = 0
        
        for atom1 in groups.atomicconnections:
                for atom2 in atom1.connectivity: 
                    if atom2[1] in groups.atoms:
                        i += 1
        
        if len(groups.groupconnection) == n:
#            if groups.size == 64:
#                print('i = {}'.format(i))
            #print((groups,groups.groupconnection))
            #if groups.size == 64:
#                group_graph(molecule,[[groups]], 'atom')
#                print(2)
            if i >= 2:
                shortest_paths = atom_chains(groups.atomicconnections)
                if shortest_paths == 'Failed_to_converage':
                    return shortest_paths
                
                for atom1 in groups.atomicconnections:
                    atom1.ringfinder -= 1
                    
                for atom2 in groups.atoms:
                    atom2.ringfinder = 0
                    
                for path in shortest_paths:
                    for atom3 in path:
                        (groups.atoms).append(atom3)
                        
                ring.append(groups.atoms)
            else:
                continue
            
        
    
    
    ringprocesser(ring, molecule)
    return None

def isolated_rings(molecule):
    """finds isolated rings and processes them using ringprocesser function"""
    for group1 in molecule.group_list:
        if len(group1.atomicconnections) == 0:
#            print(group1)
#            if group1.size == 64:
#                group_graph(molecule,[[group1]], 'atom')
#                print(0)
            #print(group1.atoms)
            for atom1 in group1.atoms:
                atom1.ringfinder = 0
            
            if group1.valence != -1 and group1.valence == 2:
                ringprocesser([group1.atoms], molecule)
                    
                
    
def group_path(group_paths, n):
    """find paths though groups that form rings"""
    for path in group_paths:
        if len(path) == n:
            return False
    return True



def group_seach(molecule):
    """finds increasingly complex rings when """
    results = None
    n = 0   
    ring_number = len(molecule.rings)
#    print(n)
    if n == 0:
        
        isolated_rings(molecule)
        if len(molecule.rings) == ring_number:
            n += 1
        
    if n == 1:
        
        simply_ring(molecule)
        if len(molecule.rings) == ring_number:
            n += 1
        
    if n == 2:
        
        results = group_one(molecule, n-1)
        
#        if len(molecule.rings) == ring_number:
#            n += 1
        
#    else:
#        group_rings, n = group_ring(group_list, molecule, n, group_rings=[])
#        groups_lists.append(group_rings)
    return results


#
#def group_ring(group_list, molecule, n, group_rings=[]):
#    """finds rings in the form of groups"""
#    max_n = len(group_list) - 2
#    group_of_v2 = [group for group in group_list if group.valence == 2]
#
#    for  group2 in group_of_v2:
#        paths = [(set([group2]), [group2])]
#        r = 0
#        while length_of_paths(paths, n, max_n):
#            new_paths = []
#            for path in paths:  
#                new_group = path[1][0].groupconnection
#                new_group2 = []
#                for group5 in new_group:
#                    if not group5 in path[0]:
#                        new_group2.append(group5)
#                new_group = new_group2
#                
#                
#                for group3 in new_group:
#                    path3 = set()
#                    for group4 in path[0]:
#                        path3.add(group4)
#                    path3.add(group3)
#                    new_paths.append((path3, [group3]))
#                
#                
#                        
#            paths = []    
#            for path3 in new_paths:
#                paths.append(path3)
#            
#            r += 1
#            
#            if r == 5:
#                
#                break
#        i = 1
#        for path5 in paths:
#            l = 0
#            for path6 in path5[0]:
#                for path7 in path6.groupconnection:
#                    if path7 in path5[0]:
#                        l += 1
#                        
#            if l == 2*len(path5[0]):
#                group_rings.append(ring_group(path5[0])) 
#        
#        k = 1
#        for path1 in paths:
#            for path2 in paths[k:]:
#                if path1[1] == path2[1]:
#                    rings = set()
#                    for path4 in [path1, path2]:
#                        for group4 in path4[0]:
#                            rings.add(group4)
#                    group_rings.append(ring_group(rings))
#            k += 1 
#
#            
#            
#    for group6 in group_rings:
#        (group6.atoms).sort(key = lambda x: x.__hash__())
#    group_rings = list(set(group_rings)) 
#    group_linkers(group_rings)
#    for group7 in group_rings:
#        for atom1 in group7.atomicconnections:
#            for group8 in group_list:
#                if atom1 in group8.atoms:
#                    group7.groupconnection.add(group8)
#                   
#    confirmed_chains = []
#    i = 1
#    j = 0
#    for group9 in group_rings:
#        for group10 in group9.groupconnection:
#            for group11 in group_rings[i:]:
#                for group12 in group11.groupconnection:
#                    if group10 == group12:
#                        confirmed_chains.append(group10)
#                        j += 1
#        i += 1
#    
#    
#    
#    if j > 0:
#        for group14 in confirmed_chains:
#            for atom2 in group14.atoms:
#                atom2.ringfinder = 0
#                atom2.ringfinders.append(atom2.ringfinder)
#            for atom3 in group14.atomicconnections:
#                atom3.ringfinder -= 1
#                atom3.ringfinders.append(atom3.ringfinder)
#            group_rings = []
#            n = 0
#    else:
#        n += 1
#     
#    return group_rings, n

def number_of_grouplinks(group_rings):
    """input group list and returns the number connection between groups in 
    that list"""
    k = 1
    count = 0
    for group1 in group_rings:
        for group2 in group_rings[k:]:
            if group2 in group1.groupconnection:
                count += 1
        k += 1
    return count
def interrconnected_rings(group_rings, group_list):
    """checks for ring_groups that are connected and returns new ring groups 
    containing both rings"""
    groups = group_rings
    group_linkers(groups)
    
    while number_of_grouplinks(groups) > 0:
        new_groups = []
        k = 1
        for group1 in groups:
            for group2 in groups[k:]:
                i = 0
                for atom in group1:
                    if atom in group2.atoms:
                        i = 1

                                    
                                
                if i == 1:
                        
                    atom_list = group1.atoms + group2.atoms
                    new_groups.append(group(atom_list, group1.valence))
                    
            k += 1
        groups = new_groups
        group_linkers(groups)
    
    for ring in groups:
        size = len(ring.atoms)
        q_position = [0,0,0]
        
        
        for atom in ring.atoms:
            
            q_position[0] += (atom.position)[0]
            q_position[1] += (atom.position)[1]
            q_position[2] += (atom.position)[2]
        
                
        q_position = np.array([x/size for x in q_position])
        
        q = FF.atom('Q', 0 , q_position, group([],0))
        q.max_valence = size
        
        for atom2 in ring:
                length = np.linalg.norm(q.position - atom2.position)
                (q.connectivity).append(length, atom2)
        (ring.atoms).append(q)
    
    return groups




   
    
def length_of_paths(paths, n, max_n):
    """finds the length of all paths"""
    length = []
    for path in paths:
        k = len(path[0])
        length.append(k)
    
   
    if len(set(length)) == 1 and length == n:
        return False
    elif length[0] == max_n:
        #print("chains not found")
        return False
    else:
        return True
    for path in paths:
        if len(path[0]) > n:
#            print('error path to long \n')
#            print(path[0])
            break
        
def ring_group(rings):
    """creates a group out of ring paths"""
    atoms = []
    for group1 in rings:
        for atom1 in group1.atoms:
            atoms.append(atom1)
    
         
    valence = 0
    for group2 in rings:
        if group2.valence > valence:
            valence = group2.valence
    
    
    return group(atoms, valence)


            

    
            
            
           
def atom_chains(atoms):
    """shortest atom paths connecting groups of valance 2.
    finds the shortest path between to atoms in a group"""
    
    start = atoms[0]
    groups = (start.groups)
    finish = atoms[-1]
    paths = [(set([start]),[start])]
    for bond in start.connectivity:
        if bond[1] == finish:
            return [[start,finish]]
    i = 0
    while convegaes_of_paths(paths, finish):
        
        for path in paths:
            
            
            new_atom = get_atoms(path[1][0].connectivity)
            
            new = set()
            for atom2 in new_atom:
                if (atom2.groups) == groups:
                    
                    if not atom2 in path[0]:
                        new.add(atom2)
            
            
            new_path =set()
            n = 0
            for atom3 in new:
                if n == 0:
                    
                    path[0].add(atom3)
                    
                    path[1][0] = atom3
                    
                    n = 1
                else:
                    for atom4 in path[0]:
                        new_path.add(atom4)
                            
                    new_path.add(atom3)
                    paths.append((new_path, [atom3]))
        i += 1
        if i == 100:
            return 'Failed_to_converage'
    shortest_paths = []
    for path in paths:
        
        if finish in path[0]:
            
            shortest_paths.append(list(path[0]))
    shortest = len(shortest_paths[0])
    for path2 in shortest_paths:
        if len(path2) < shortest:
            shortest = len(path2)
    shortest_paths = [path3 for path3 in shortest_paths if len(path3) == shortest]
    
    return shortest_paths
    


def atoms_in_rings(molecule):
    """finds atoms that are in mutpial rings"""
    j = 1
    same_atoms =[]
    for ring1 in molecule.rings:
        for atom5 in ring1:
            for ring2 in molecule.rings[j:]:
                for atom6 in ring2:
                    if atom5 == atom6:
                        same_atoms.append([atom5, [ring1, ring2]])   
        j += 1



def aromatic_group(molecule):
    """connects all bonded aromatic atoms into groups"""
    aromate_group = []
    for atom1 in molecule.aromatic:
        atoms = {atom1}
        n = 0
        current_atoms = [atom1]
        while len(atoms) > n:
            n = len(atoms)
            for atom2 in current_atoms:
                for atom3 in get_atoms(atom2.connectivity):
                    if atom3 in molecule.aromatic:
                        atoms.add(atom3)
            current_atoms = [atom4 for atom4 in atoms]
        list(atoms)
        aromate_group.group(atoms, 0)
        
        for group1 in aromate_group:
            (group1.atoms).sort(key = lambda x: x.__hash__())
        molecule.aromatic = list(set(aromate_group))
    




def interconnected_rings(molecule):
    """Connects rings that have Q atoms that I with each other's sphere defined
    by longest bond from the Q atoms generating a new ring with all atoms and 
    wondering and 1Q atoms. """
    longest_bond = 0
    interconnected_rings = []
    for ring1 in molecule.rings:
        for atom in ring1:
            if atom.symbol == 'Q':
                (atom.connectivity).sort(key = lambda x: x[0])
                longest_bond = atom.connectivity[-1][0]
                for atom1 in molecule.atoms:
                    lenght = np.linalg.norm(atom.position - atom1.position)
                    if (lenght < longest_bond) and (atom1.symbol == 'Q'):
                        for ring2 in molecule.rings:
                            if atom1 in ring2:
                                interconnected_rings.append([ring1+ring2])
    
    for rings in interconnected_rings:
        for atom3 in rings:
            if atom3.symbol == 'Q':
                rings.remove(atom3)
    ringprocesser(interconnected_rings, molecule)                                





def convegaes_of_paths(paths, finish):  
    """finds the shortest path between to atoms in a group. finds if the finial
    atom is in once of the groups"""      
    
    for path in paths:
        
        if finish in path[0]:
            
            return False
    
    return True
    
