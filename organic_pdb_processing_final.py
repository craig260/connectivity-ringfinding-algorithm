# -*- coding: utf-8 -*-
"""
Created on Tue Apr  3 11:10:56 2018

organic molecules reading pdb code for the 

@author: Craig Melton
"""

import MD_main_final_2 as MD
#import force_field_code_final as FF
import expected_ring_number as EV
import time as t
import pdb_formating_data as pdb
#import connectivity_graph3 as cg
import ringalgothem_final_2 as ring

from mpl_toolkits.mplot3d import Axes3D
import matplotlib.pyplot as plt

#class Expected_rings:
#    """finds and creates groups of atom that expected ring where not found in"""
#    def __init__(self, amino_acid):
#        """amino_acid = ([atom_ID's],piptide(str)) """
#        self.amino_acid = amino_acid
#        self.expected = 0
#        self.confirmed = False
#        self.found = 0
#        self.group = ring.group([], -10)
#        
##    def check_for_rings(self):
##        """find and calculates the number of rings in the amino acid"""
##        for atom in self.group.atoms:
##            if atom.symbol == 'Q':
##                self.found += 1
##            
##            if self.found == self.expected:
##                self.confirmed = True
##        return None
expected_ringnumber = {'HIS': 1, 'PRO': 1, 'PHE': 1, 'TYR': 1, 'TRP':2}        
class Expected_rings:
    """finds and creates groups of atom that expected ring where not found in"""
    def __init__(self, amino_acid):
        """amino_acid = ([atom_ID's],piptide(str)) """
        self.amino_acid = amino_acid
        self.atoms = []
        self.expected = []
        self.confirmed = False
        self.found = 0
        self.group = ring.group([], -10)
        self.q_atoms = []
        self.bonds = []
        self.rings = []
        
    def check_for_rings(self,molecule):
        """find and calculates the number of rings in the amino acid"""
#        for atom in self.group.atoms:
#            if atom.symbol == 'Q':
#                self.found += 1
        
        
        try:
            self.expected = EV.expected_ringnumber[self.amino_acid[1]]
        except KeyError:
            self.expected = 0 
        if self.found >= self.expected[0]:
            self.confirmed = True
        
        return None
atomtype = {'C': 'b', 'O': 'r', 'N': 'g', 'H': 'w', 'S': 'y'}       
  
def unfound_rings(molecule,amino_acids):
    """create groups of amino_acids with unfound ring return list if Expected_rings
    class with these group appended"""
    groups1 = []
    try:    
        for amino_acid1 in amino_acids:
            atoms = []
        
            for ID in amino_acid1.amino_acid[0]:
                for atom1 in molecule.atoms:
                    if ID == atom1.ID:
                        atoms.append(atom1)
                    
                for atom1 in molecule.atoms:
                    for bond2 in atom1.connectivity:
    #                        for atom3 in bond2.atoms:
    ##                            atoms.append(atom3)
                        if atom1.symbol == "Q":
                            n = 0
                            for bond1 in atom1.connectivity:
                                for atom2 in bond1.atoms:
                                    if atom2 in atoms:
                                        n = 1
                            if n == 1:
                                atoms.append(atom1) 
                amino_acid1.group = ring.group(atoms, 20)
            if len(atoms) > 0:
                groups1.append(ring.group(atoms, 20))
    except AttributeError:
        print(amino_acids)
    return amino_acids

def pdb_file_reader(name):
    """read the data from the pdb file and then return atoms data in nested 
    list of lists"""
    data = open(name)
        
            
    xyz_data = []
    line = data.readline()
    
    check = line[:6].strip()
   
    
    expected = 0
    i = int(line[22:26])
    
    amino_acids = []
    amino_acid = []
    
    while check != 'END':
        
        if (check ==  'ATOM') or (check == ('HETATM')):
            datalist = [line[76:78].strip(),float(line[30:38]),float(line[38:46]),float(line[46:54]),line]
            piptide = line[17:20].strip().upper() 
#            if i == -1:
#                #piptide = line[17:20].strip().upper() 
#                id1 = int(line[6:11])
#                amino_acid.append(id1)
            
            if int(line[22:26]) == i:
                id1 = int(line[6:11])
                amino_acid.append(id1)
                
            
            
            xyz_data.append(datalist)

            
            
            
        line = data.readline()
        
        check = line[:6].strip()
        if (check ==  'ATOM') or (check == ('HETATM')):
            if int(line[22:26]) != i:
                amino_acids.append((amino_acid,piptide))
                amino_acid = []
                id1 = int(line[6:11])
                amino_acid.append(id1)
            i = int(line[22:26])
    amino_acids.append((amino_acid,piptide))
    amino_acids_2 = []
    for amino_acid in amino_acids:
        amino_acids_2.append(Expected_rings(amino_acid))
#    print(amino_acids)
    ring_amino_acids = []
    for tuple1 in amino_acids_2:
        try:
            if EV.expected_ringnumber[tuple1.amino_acid[1]][0] > 0 :
                tuple1.expected = EV.expected_ringnumber[tuple1.amino_acid[1]]
                
                ring_amino_acids.append(tuple1)
        except KeyError:
            continue
    for peptide1 in ring_amino_acids:
        expected += EV.expected_ringnumber[peptide1.amino_acid[1]][1]
#    print(ring_amino_acids)
    data.close()
    return xyz_data, ring_amino_acids

def ring_sizes(molecule):
    """takes atoms from molecule and then returns a list directery with the 
    number of each atom type present"""
    ring_sizes = {}
    for ring1 in molecule.rings:
        ring_sizes[len(ring1)] = 1 + ring_sizes.get(len(ring1), 0)
    return ring_sizes

def ring_sizes_errors(molecule):
    """takes atoms from molecule and then returns a list directery with the 
    number of each atom type present"""
    ring_sizes = {}
    for ring1 in molecule.ring_errors:
        ring_sizes[len(ring1)] = 1 + ring_sizes.get(len(ring1), 0)
    return ring_sizes

def amino_acide_atoms(molecule, ID_list, expected):
    """creates list of atom in amino acids that have rings"""
#    amino_acid = []
#    for IDs in ID_list:
        
#        atoms = []
#        for atom in molecule.atoms:
#            if atom.ID in IDs[0]:
#                atoms.append(atom)
#        
#        amino_acid.append(atoms)
    connected_rings = []
    k = 0
    for ring1 in molecule.rings:
        k += 1
        for ring2 in molecule.rings[k:]:
            if ring1 != ring2:
                if any(i in ring1 for i in ring2):
                    ring3 = ring1 + ring2
                    connected_rings.append(list(set(ring3)))
    for ring2 in connected_rings:
        print(len(ring2))
    for ring1 in  molecule.rings:
        connected_rings.append(ring1)
    connected_ID = []
    for ring1 in  connected_rings:
        ring_id = []
        for atom1 in ring1:
            if atom1.symbol != 'Q':
                ring_id.append(atom1.ID)
        connected_ID.append(ring_id)
    
    confirmed_rings = 0
    for piptide in ID_list:  
        for ring2 in connected_ID:
            if len(ring2) == piptide.expected[1]:
                if all(i in piptide.amino_acid[0] for i in ring2):
                    
                    if not ring2  in piptide.rings:
                        piptide.rings.append(ring2)
                        piptide.found += 1
                        
    for piptide in ID_list:
        confirmed_rings += piptide.found
                    
    print('connected_rings {} , rings {}\n'.format(len(connected_rings),len(molecule.rings)))
#    for ring in connected_rings:
#        print(len(ring))
    
    
        
#    for piptide in ID_list:
#        for ring1 in molecule.rings:
#            size = 0
#            for atom1 in ring1:
#                if atom1.symbol != 'Q':
#                    size += 1
#            num = 0
#            if size == piptide.expected[1]:
#                for atom1 in ring1:
#        #                print(atom1.ID in piptide[0])
#        #                print('{} in {}'.format(atom1.ID, piptide[0]))
#                    if atom1.ID in piptide.amino_acid[0]:
#                        num += 1
#                    if num == size:
#                        piptide.found += 1
##                        if piptide.found >= piptide.expected[0]:
##                            piptide.confirmed = True
#        if num == size:
#            confirmed_rings += 1
    
    
    
    confirmed2 = 0     
    bool_list = []
    for piptide in ID_list:
        confirmed2 += piptide.found
        piptide.check_for_rings(molecule)
        bool_list.append(piptide.confirmed)
        
    print('confirmed {} == confirmed2 {}'.format(confirmed_rings,confirmed2))
    
    if all(bool1 == True for bool1 in bool_list):  # confirmed_rings >= expected:
        return (True, confirmed_rings, ID_list)
    else:
        return (False, confirmed_rings, ID_list)
    return None
  
def ring_errorhander(molecule, ID_list, expected, error_numbers, name, ring_errors, ring_dic):
    """takes in the simple ring data from the molecule and writes out and counts  
    the errors into a file"""
    small_ring = ring_sizes_errors(molecule[0])
    Number_invaild_rings = 0
    if len(molecule[0].ring_errors) > 0:
        
        for i in list(small_ring.keys()):
            
            error_numbers[4] += small_ring[i]
            Number_invaild_rings += small_ring[i]
        error_numbers[5] += 1
    num = 0
    for i in [0,1,2,3,8]:
        try:
            if ring_dic[i] > 0:
                num += ring_dic[i]
        except KeyError:
            num += 0 
            
    text2 = 'ring_dic {}\n'.format(ring_dic)
    text2 += 'amino_acids_dic {}\n'.format(molecule[0].amino_acids_dic)
    text2 += 'ring_cycles {}\nnumber_of_atoms {}\n'.format(molecule[0].ring_cycles, len(molecule[0].atoms))
    if num > 0:
        text2 += 'impossible_rings_found {}\n'.format(num)
        error_numbers[4] += 1
    rings_expected = amino_acide_atoms(molecule[0], ID_list, expected)
    text2 += 'rings number found {} and expected {}\nconfirmed_rings found {}\n'.format(len(molecule[0].rings),expected,rings_expected[1])
    text2 += "error_in_ring_size_dic {} total number {}\n".format(small_ring , Number_invaild_rings)
    text2 += 'carbonyl_group {}\n'.format(len(molecule[0].carbonyl_group))
    if rings_expected[0]:
        if len(molecule[0].rings) == expected:
            header = 'No_ring_error_found\nconfirmed_rings found {}\n'
            
        elif len(molecule[0].rings) > expected:
            header = 'Minor_ring_error (extra_rings_found)\n'
            header2 = 'file name\n{}\n'.format(name)
            #print(str(header) + str(header2) + str(text2))
            ring_errors.write(header + header2 + text2)
            error_numbers[0] += 1
            error_numbers[1] += 1
        
    else:
        if molecule[6] == 'TimeoutError':
            header = 'Expected_ring_is_not_found_with_timeout_error\n'
            header2 = 'file name\n{}\n'.format(name)
            ring_errors.write(header + header2 + text2)
            error_numbers[0] += 1
            error_numbers[6] += 1
            
        elif molecule[6] == 'Failed_to_converage':
            header = 'Expected_ring_is_not_found_with_Failed_to_converage\n'
            header2 = 'file name\n{}\n'.format(name)
            ring_errors.write(header + header2 + text2)
            error_numbers[0] += 1
            error_numbers[7] += 1
        
        elif len(molecule[0].rings) > expected:
            header = 'Expected_ring_is_not_found (extra_rings_found)\n'
            header2 = 'file name\n{}\n'.format(name)
            ring_errors.write(header + header2 + text2)
            error_numbers[0] += 1
            error_numbers[2] += 1
            
        elif len(molecule[0].rings) <= expected:
            header = 'Expected_ring_is_not_found\n'
            header2 = 'file name\n{}\n'.format(name)
            ring_errors.write(header + header2 + text2)
            error_numbers[0] += 1
            error_numbers[3] += 1
            
    
    return header + text2, rings_expected[2]
                
grouptype = {2: 'b', 1: 'r', 3: 'g', 0: 'w', 4: 'y'}
bondtypes = {1:'w-',2:'c-',3:'K-',1.5:'w--',0:'k--'}

def main(filename, errors):
    """"takes filename of file contain list of names and then readin this names 
    into the file processes the files one by one."""
    timestart = t.time()
    pdbs = open(filename)
    
    pdb_names = pdbs.readlines()
    pdbs.close()
#    pdb_names2 = []
    
    ring_errors = open(r'C:\Users\Craig_Melton\Desktop\Final 4\organic trails\trial_9\ring_errors1.error','w')
    error = open(r'C:\Users\Craig_Melton\Desktop\Final 4\organic trails\trial_9' + '\\' + errors,'w')
    results = open(r'C:\Users\Craig_Melton\Desktop\Final 4\organic trails\trial_9\organic_connectivity1.results','w')
    
#    ring_errors = open(r'C:\Users\Craig\Downloads\Final 2\organic trails\test trails\testing new confirmed ring method\ring_errors6.error','w')
#    error = open(r'C:\Users\Craig\Downloads\Final 2\organic trails\test trails\testing new confirmed ring method' + '\\' + errors,'w')
#    results = open(r'C:\Users\Craig\Downloads\Final 2\organic trails\test trails\testing new confirmed ring method\organic_connectivity6.results','w')
    
    j = 0
    i = 0
    
    total_complete = 0
    total_exstinened_convergance = 0
    total_errors = 0
    timeout_errors  = 0
    failed_to_converage = 0
    
    totalring_errors = 0
    more_ring_with_expected = 0
    more_ring_without_expected = 0
    to_few_rings = 0
    invald_rings_found = 0
    invald_molecules_rings = 0
    timeout_error_without_expected = 0
    Failed_to_converage_ringerror = 0
    error_numbers = [totalring_errors,more_ring_with_expected,more_ring_without_expected,to_few_rings, invald_rings_found,invald_molecules_rings,timeout_error_without_expected,Failed_to_converage_ringerror]
    
   
    
    for pdb_name in pdb_names:
        # processing code 
#        pdb_name = r'C:\Users\Craig\Downloads\Final 2\bondtype_testing2.pdb' #r'C:\Users\Craig\Downloads\Final 2\NMR_structure_library_pdb\6azf.pdb' #5twi
        name = r'C:\Users\Craig_Melton\Desktop\Final 4\NMR_structure_library_pdb' + '\\' + pdb_name.rstrip('\n')#r'C:\Users\Craig\Desktop\final\NMR_structure_library_pdb' + '\\'+ pdb_name.rstrip('\n') # r'C:\Users\Craig\Documents\work\New folder\connectivity_code\NMR_structure_library_pdb' +'\\' +
#        name = pdb_name.rstrip('\n')
        
        print(name)
        
        xyz_data, amino_acids1 = pdb_file_reader(name)
#        for amino_acid1 in  amino_acids1:
#            print(amino_acid1.amino_acid)
        
#        molecule, n = FF.read_in_molecule_data_pdb(xyz_data)
#        molecule.connectivity(n)
        
        result = MD.Main(xyz_data)
        ring_dic = ring_sizes(result[0])
        result[0].amino_acids(amino_acids1)
        
        expected = 0
        for key in list(result[0].amino_acids_dic.keys()):
            try:
                expected += (EV.expected_ringnumber[key][0])*result[0].amino_acids_dic[key]
            except KeyError:
                expected += 0
        
#        results.append((result, pdb_name2))
#        print('expected {} and i {}'.format(expected, i))
        
        j += 1
#        print('{} index'.format(j))
        
        
#        if len(result[0].rings) != expected:
#            ring_errors.write('file name\n{}\n'.format(name))
#            ring_errors.write('rings number found {} and expected {}\n'.format(len(result[0].rings),expected))
#            incorrect_ring_count += 1
            
        
            
#     write errors to files
        if result[6] == 'Failed_to_converage':
            text = 'Failed_to_converage\n'
            failed_to_converage += 1
            total_errors += 1
        else:
            if result[6] == 'TimeoutError':
                text = 'TimeoutError\n'
                timeout_errors  += 1
                total_errors += 1
            
            elif result[5] > 60:
                text = 'Extendedconvergence\n'
                total_exstinened_convergance += 1
            else:
                text = 'succefull_convergation\n'
                total_complete += 1
        
        text += 'file {}'.format(j) + '\npdb_name \n' +str(name) + '\nerror \n' 
        text += str(result[2]) +'\n' + 'ring_time  ' + str(result[3]) + ' connectivity_time  ' 
        text += str(result[4]) +' total_time ' + str(result[5]) +'\n'
        ring_errors1 = ring_errorhander(result, amino_acids1, expected, error_numbers, name, ring_errors, ring_dic)
        amino_acids = unfound_rings(result[0],ring_errors1[0])
        text += ring_errors1[0] + '\n'
        if result[6] == 'Failed_to_converage':
            end = '.FTC'
            file = open(r'C:\Users\Craig_Melton\Desktop\Final 4\organic trails\trial_9' + '\\' + pdb_name.split('\\')[-1].rstrip('\n') + end,'w')#r'C:\Users\Craig\Desktop\final\succefull_convergation_under' + '\\' + pdb_name.rstrip('\n') + end,'w')
            file.write(text)
            text2 = pdb.pdb_formated_data(result[0])
            file.write('\n')
            file.write(text2)
            file.close()
            error.write(text)
        else:
            if len(result[2]) > 0:
                end = '.Timeout'
                file = open(r'C:\Users\Craig_Melton\Desktop\Final 4\organic trails\trial_9' + '\\' + pdb_name.split('\\')[-1].rstrip('\n')  + end,'w')#r'C:\Users\Craig\Desktop\final\succefull_convergation_under' + '\\' + pdb_name.rstrip('\n') + end,'w')
                file.write(text)
                text2 = pdb.pdb_formated_data(result[0])
                file.write('\n')
                file.write(text2)
                file.close()
                error.write(text)
            elif result[5] > 60:
                end = '.EC'
                file = open(r'C:\Users\Craig_Melton\Desktop\Final 4\organic trails\trial_9' + '\\' + pdb_name.split('\\')[-1].rstrip('\n')  + end,'w')#r'C:\Users\Craig\Desktop\final\succefull_convergation_under' + '\\' + pdb_name.rstrip('\n') + end,'w')
                file.write(text)
                text2 = pdb.pdb_formated_data(result[0])
                file.write('\n')
                file.write(text2)
                file.close()
                error.write(text)
            
            else:
                end = '.SC'
                file = open(r'C:\Users\Craig_Melton\Desktop\Final 4\organic trails\trial_9' + '\\' + pdb_name.split('\\')[-1].rstrip('\n')  + end,'w')#r'C:\Users\Craig\Desktop\final\succefull_convergation_under' + '\\' + pdb_name.rstrip('\n') + end,'w')
                file.write(text)
                text2 = pdb.pdb_formated_data(result[0])
                file.write('\n')
                file.write(text2)
                file.close()
        
        
            
        print('\n'+ text + '\n')
        
        results.write(text)
#        
        
    finishtime = t.time()
    geraition = finishtime - timestart
    print(' time {} s\n'.format(geraition))
    # write state of the result of the trial
    stats = 'stats\n' + 'TimeoutError_total ' + str(timeout_errors) + '  '  
    stats += 'Failed_to_converage ' + str(failed_to_converage) + '  '
    stats += 'Total_Errors '+ str(total_errors) + '  '
    stats += '\nExtendedconvergence ' + str(total_exstinened_convergance) + '  '
    stats += 'succeful_convergation ' + str(total_complete) + '\n'
    stats += 'invald_rings_found_ {}\n'.format(invald_rings_found)
    stats += 'numder_of_molecules_with_invaild_rings {}\n'.format(invald_molecules_rings)
    stats += 'File Number ' + str(total_complete + total_exstinened_convergance + total_errors)
    stats += '\nEnd'
    results.write(stats)
    results.write(' time {} s\n'.format(geraition))
    results.close()
    
    stat = 'stats\n' + 'TimeoutError_total ' + str(timeout_errors) + '  '  
    stat += 'Failed_to_converage ' + str(failed_to_converage) + '  '
    stat += 'Total_Errors '+ str(total_errors) + '\n'
    error.write(stat)
    error.write(' time {} s\n'.format(geraition))
    error.close()
    
    ring_errors.write('\n Ring_Count_Error_Number_Stat \n')
    stat2 = 'total_errors {} \nMinor_ring_error (extra_rings_found) {}\n'
    stat2 += 'Expected_ring_is_not_found (extra_rings_found) {}\nExpected_ring_is_not_found {}\n'
    stat2 += 'less_an_expected_with_timeout_error {}\nFailed_to_converage_ringerror {}\n'
    stat3 = stat2.format(error_numbers[0],error_numbers[1],error_numbers[2],error_numbers[3],error_numbers[6],error_numbers[7])
    ring_errors.write(stat3)
    ring_errors.close()
    
    
    
    return geraition



time = main(r'C:\Users\Craig_Melton\Desktop\Final 4\pdb_files.list', 'organic_errors1.errors')
#time = main(r'C:\Users\Craig\Desktop\Final 2\organic trails\trail_5_results\Failed_to_converage.list',1502,1503, 'organic_errors6.errors')

#time, molecules = main('pdbs.list', 0, -1, 'errors_find.errors')


#xyz_data, amino_acids, expected = pdb_file_reader(r'C:\Users\Craig\Desktop\final\NMR_structure_library_pdb\5yfg.pdb')
#        
#result = MD.Main(xyz_data)
#ring_dic = ring_sizes(result[0])
#cg.group_graph(result[0],result[0].group_list, 'atom')        




#        pdb_names2.append(name)
#    
#    
#    for pdb_name2 in pdb_names2[n:m]:
#        print(name)





#num, time = main(r'C:\Users\Craig\Documents\work\New folder\connectivity_code\pdb_files.list',0,100, 'errorsbtw_0to100.errors')
#num1, time1 = main(r'C:\Users\Craig\Documents\work\New folder\connectivity_code\pdb_files.list',100,200, 'errorsbtw_100to200.errors')
#num2, time2 = main(r'C:\Users\Craig\Documents\work\New folder\connectivity_code\pdb_files.list',200,300, 'errorsbtw_200to300.errors')
#num3, time3 = main(r'C:\Users\Craig\Documents\work\New folder\connectivity_code\pdb_files.list',300,400, 'errorsbtw_300to400.errors')
#num4, time4 = main(r'C:\Users\Craig\Documents\work\New folder\connectivity_code\pdb_files.list',400,500, 'errorsbtw_400to500.errors')
#num5, time5 = main(r'C:\Users\Craig\Documents\work\New folder\connectivity_code\pdb_files.list',500,600, 'errorsbtw_500to600.errors')
#num6, time6 = main(r'C:\Users\Craig\Documents\work\New folder\connectivity_code\pdb_files.list',600,700, 'errorsbtw_600to700.errors')
#num7, time7 = main(r'C:\Users\Craig\Documents\work\New folder\connectivity_code\pdb_files.list',700,800, 'errorsbtw_700to800.errors')
#num8, time8 = main(r'C:\Users\Craig\Documents\work\New folder\connectivity_code\pdb_files.list',800,900, 'errorsbtw_800to900.errors')




#num10, time10 = main(r'C:\Users\Craig\Documents\work\New folder\connectivity_code\pdb_files.list',1000,-1, 'errorsbtw_1000to-1.errors')
#num9, time9 =
#print(num)
#cg.main(num[0][0][0])
##
#cg.group_graph(num[0][0][0], [num[0][0][0].group_list], 'group')


#C:\Users\Craig Melton\Documents\New folder\connectivity_code\NMR_structure_library_pdb\1bf0.pdb


#    if line[0] == 'TER':
#        molecule, n = FF.read_in_molecule_data_pdb(xyz_data)
#            result = MD.Main(molecule, n)
#        if result == 1:
#            errors.append(pdb_name[0])
#        else:
#            molecules.append(result)
#       
#    result = MD.Main(pdb_name[0])
#    results = []   
#    for pdb_name in pdb_names2:
#        result = MD.Main(pdb_name)
#        results.append(result)
#
#   
#    M,n= FF.read_in_molecule_data_pdb(pdb_name[0])