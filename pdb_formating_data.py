# -*- coding: utf-8 -*-
"""
Created on Fri May 11 14:17:49 2018

gerenates pdb file info for data to be 

@author: Craig
"""

def pdb_formated_data(molecule):
    """takes molecule class and generates string in pdb format to be written 
    to file"""
    # ATOM line
    atom_lines = ''
    conect_line = ''
    
    for atom in molecule.atoms:
        
        atom_lines += atom.pdb
#        atom_lines += '\n'
        
        atoms2 = []
        for bond in atom.connectivity:
            try:
                atoms2.append(bond.atoms[1])
            except AttributeError:
                atoms2.append(bond[1])
        line = 'CONECT{:>5}'.format(atom.ID)
        for atom2 in atoms2:
            line += '{:>5}'.format(atom2.ID)
        conect_line += line
        conect_line += '\n'
        
    text = atom_lines + conect_line
        
    return text
        
# if # something
#            line = 'ATOM  {:>5}{:^4}{:>3}{:1}{:3}{:1}{:4}{:1}   {:>8}{:>8}{:>8}'
#            line += '{:>6}{:>6}' +' '*10 + '{:>2}{:>2}'
#            line = line.format(atom.ID, atom.name,Character,Residue,chainID,resSeq,iCode,x,y,z,occupancy,tempFactor,atom.symbol,charge)
#            atom_lines.append(line) 


#    
#    for line in atom_lines:
#        text += line
#        
#    for line in conect_line:
#        text += line      