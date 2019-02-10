# -*- coding: utf-8 -*-
"""
Created on Sat Jan 27 20:40:39 2018

controls the force_field_code and ringalgothem to automaticely process molecules
for MD simulations.

@author: Craig Melton
"""

from mpl_toolkits.mplot3d import Axes3D
import matplotlib.pyplot as plt
import time as t

import expected_ring_number as EV
import force_field_code_final as MD
import ringalgothem_final_2 as ring
#import numpy as np

#pdb reader

#def ring_number_lists(molecule,Residue_number=[],Residue=[]):
#    """creates lists for ringfinder number for each Residue"""
#    
#    if len(Residue_number) < 1:
#        Residue = []
#        for atom in molecule.atoms:
#                Residue.append(atom.pdb[22:26])
#        Residue = list(set(Residue))
#        
#        for res_mun in Residue:
#            Residue_number.append((res_mun,[]))
##        print(Residue, ' ', Residue_number)
#    else:
#        for atom in molecule.atoms:
#            for res_mun in Residue_number:
#                if atom.pdb[22:26] == res_mun[0]:
#                    if len(res_mun[1]) < 1:
#                        res_mun[1].append((atom.ID,[atom.ringfinder]))
#                    else:
#                        if res_mun[0] == atom.ID:
#                            res_mun[1][1].append(atom.ringfinder)
#    
#    
#    return Residue_number, Residue
                        

            
# main function

def Main(file):#amino_acids):
    """This function brings together the force_field_code and ringalgorithm scripts 
    to fully automatedly processor molecule for molecular dynamics simulations """
    
    cycle_times =[]
    chainremovials =[]
    groupcreations =[]
    group_setups =[]
    ring_finding_times =[]
    times_lists = [cycle_times,chainremovials,groupcreations,group_setups,ring_finding_times]
    group_nums =[]
    
    start = t.time()
    t1 = start
    last_time = t1
    
    
    
    
    molecule, n = MD.read_in_molecule_data_pdb(file)
    
    print('atoms # {}'.format(len(molecule.atoms)))
    print('timeout time {}'.format((60*2 + 60*(len(molecule.atoms)/1000)**2)))
    
    bonds, errors = molecule.connectivity(n)
#    molecule.Hydrogen_duplite_remove()
    
    molecule.non_H_connections_num()
    molecule.Extendedconnections()
    
    
    
    t2 = t.time()
    connectivity_time = t2 - last_time
    last_time = t2
    
    
    print(connectivity_time)
    
    ring.connections_num(molecule)
    
#    Residue_number, Residue = ring_number_lists(molecule)
    
    
    groups = []
    
    error_list = []

    groups_lists = []
    
    non_zero = []
    for atom1 in molecule.atoms:
        if atom1.ringfinder > 0:
            non_zero.append(atom1)
    ring.carbonyl_groups(molecule)
    j = 0
    
    
    for atom in molecule.atoms:
            atom.ringfinder_list.append(atom.ringfinder)
    
    
    
    
    t3 = t.time() 
    ringsetup = t3 - last_time
    last_time = t3
    
    
    while len(non_zero) > 0 :
        print(len(non_zero))
        cycle_start = t.time()
        last_time = cycle_start
        
        ring.chain_end(molecule)
        
        pesent_time = t.time()
        chainremovial = pesent_time - last_time
        times_lists[1].append(chainremovial)
        last_time = pesent_time
        
        for atom in molecule.atoms:
            atom.ringfinder_list.append(atom.ringfinder)
        
        
        
        molecule.group_list, groupnum = ring.homogesis_atom_connections(molecule)
        group_nums.append(groupnum)
        
        pesent_time = t.time()
        groupcreation = pesent_time - last_time
        times_lists[2].append(groupcreation)
        last_time = pesent_time
        
        atoms_2 = []
        for atom in molecule.atoms:
            if atom.ringfinder == 0:
                atoms_2.append(atom)
        molecule.group_list.append(ring.group(atoms_2, -1))
        
        groups_lists.append(molecule.group_list)
        
       
        
        ring.atom_groups(molecule.group_list)
        
        
        ring.group_linkers(molecule.group_list)
        
        pesent_time = t.time()
        group_setup = pesent_time - last_time 
        times_lists[3].append(group_setup)
        last_time = pesent_time

        error = ring.group_seach(molecule)
        
        pesent_time = t.time()
        ring_finding_time = pesent_time - last_time 
        times_lists[4].append(ring_finding_time)
        last_time = pesent_time
        
        for atom in molecule.atoms:
            atom.ringfinder_list.append(atom.ringfinder)
        
        molecule.ring_cycles += 1
        
        cycle_finish = t.time()
        cycle_time = cycle_finish - cycle_start
        times_lists[0].append(cycle_time)
        last_time = cycle_finish
        
        if error == 'Failed_to_converage':
            t4 = t.time()
            ring_time = t4-t3
            molecule.metal_connections()
            
            MD.bond_hack(molecule)
            MD.bond_processing(molecule)
            
            
            
            finish = t.time()
            total_time = finish - start
            return (molecule, ring.count_non_zeroatoms(molecule), error_list, ring_time, connectivity_time, total_time, 'Failed_to_converage',ringsetup,times_lists,group_nums)

        non_zero = []
        for atom1 in molecule.atoms:
            if atom1.ringfinder > 0:
                non_zero.append(atom1)
        

#                for atom1 in atoms:
#                    for  bond2 in atom1.connectivity:
#                        for atom3 in bond2.atoms:
#                            atoms.append(atom3)
#                atoms = list(set(atoms))
                    
            
        
        
        t4 = t.time()
        ring_time = t4-t3
        if molecule.group_list != groups:
            groups = molecule.group_list
        else:
            j += 1
            if j == 2:
                print('shorter limit used')
                for group in molecule.group_list:
                    if group.valence > 0:
                        for atom in group.atoms:
#                            if atom.symbol in ['C', 'N', 'O']:
                            for bond in atom.connectivity:
#                                    if bond[1].symbol in ['C', 'N', 'O']:
                                if bond[0] >= 2.00: 
                                    if bond[1].symbol != 'H':                             
                                        atom.ringfinder -= 1
                                    
#                for atom in molecule.atoms:
#                    for bond in atom.connectivity:
#                        if bond[1].ringfinder > 0:
#                            atom.ringfinder2 += 1
#                
#                for atom in molecule.atoms:
#                    if atom.symbol != 'H':
#                        atom.ringfinder = atom.ringfinder2
                    
                for atom in molecule.atoms:
                    atom.ringfinder_list.append(atom.ringfinder)
                
            if j == 6:
#                print('Failed_to_converage\n')
                molecule.metal_connections()
                
                MD.bond_hack(molecule)
                MD.bond_processing(molecule)
                
                
                
                finish = t.time()
                total_time = finish - start
                return (molecule, ring.count_non_zeroatoms(molecule), error_list, ring_time, connectivity_time, total_time, 'Failed_to_converage',ringsetup,times_lists,group_nums)
        print('current_ring_time {}'.format(t4-t3)) 
#        if t4-t3 > (60*2 + 60*(len(molecule.atoms)/1000)**2):
#            print(t4-t3)
##            print("timed out error")
#            error_list.append((t4-t3,"timed out error"))
#            molecule.metal_connections()
#            MD.bond_hack(molecule)
#            
##            timeout = ring.count_non_zeroatoms(molecule)
#            finish = t.time()
#            total_time = finish - start 
#            return (molecule, ring.count_non_zeroatoms(molecule), error_list, ring_time, connectivity_time, total_time, 'TimeoutError')
       
    molecule.metal_connections()
    
    MD.bond_hack(molecule)
    MD.bond_processing(molecule)
    
    
    for atom in molecule.atoms:
            atom.ringfinder_list.append(atom.ringfinder)
    
    try:
        ring_time = t4-t3
    except UnboundLocalError:
        t4 = t.time()
        ring_time = t4-t3
    finish = t.time()
    total_time = finish - start 
   
    
    return (molecule, ring.count_non_zeroatoms(molecule), error_list, ring_time, connectivity_time, total_time, 'succefull_convergation',ringsetup,times_lists,group_nums) #(molecule, groups_lists)



grouptype = {2: 'b', 1: 'r', 3: 'g', 0: 'w', 4: 'y', 5:'m'}

def main(molecule, groups_list):
    """construting connectivit graph"""
    for group_list in groups_list:
        plt.style.use('dark_background')
        fig = plt.figure()
        ax = fig.add_subplot(111, projection='3d')
        
        
        for group in group_list:
            colour = grouptype[group.valence]
            for atom in group.atoms:
                x = atom.position[0]
                y = atom.position[1]
                z = atom.position[2]
            
                ax.scatter(x, y, z, c=colour ,marker='o')
            
        
        
        for bond in list(molecule.bonds):
            try:
                x1 = bond.atoms[0].position[0]
                y1 = bond.atoms[0].position[1]
                z1 = bond.atoms[0].position[2]
                x2 = bond.atoms[1].position[0]
                y2 = bond.atoms[1].position[1]
                z2 = bond.atoms[1].position[2]
                ax.plot([x1,x2], [y1,y2], [z1,z2], 'w-')
    #            if bond.atom[0].symbol == bond.atom[1].symbol:
    #                color = atomtype[bond.atom[0].symbol] + '-'
    #                ax.plot([x1,x2], [y1,y2], [z1,z2], color)
    #            else:
    #                color1 = atomtype[bond.atom[0].symbol] + '-'
    #                color2 = atomtype[bond.atom[1].symbol] + '-'
    #                ax.plot([x1,(x2-x1)/2], [y1,(y2-y1)/2], [z1,(z2-z1)/2], color1)
    #                ax.plot([(x2-x1)/2,x2], [(y2-y1)/2,y2], [(z2-z1)/2,z2], color2)
            except AttributeError:
                continue
        ax.set_xlabel('x')
        ax.set_ylabel('y')
        ax.set_zlabel('z')
        
        
        plt.show()
    return None

#Asn-Gln-Leu-Pro-Tyr2.xyz
# new_tests1.xyz
# test2.xyz
#test1  = Main('Asn-Gln-Leu-Pro-Tyr2.pdb')
#
#main(test1, group)



def group_brakedown(group_list):
    """takes atoms from molecule and then returns a list directery with the 
    number of each atom type present"""
    group_types = {}
    for group1 in group_list:
        group_types[(group1.size,group1.valence)] = 1 + group_types.get((group1.size,group1.valence), 0)
    return group_types



def bond_brakedown(atoms):
    """ """
    i = 1
    bond_types = {}
    for atom1 in atoms:
        for atom3 in atom1.connectivity:
                for atom2 in atoms[i:]:
                    if atom3 == atom2:
                        for atom4 in atom2.connectivity:
                            if atom4 == atom1:
                                symbols = [atom1.symbol, atom2.symbol]
                                symbols.sort()
                                bond = '{}-{}'.format(symbols)
                                bond_types[bond] = 1 + bond_types.get(bond, 0)






