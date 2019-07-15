"""
aimall_utils.py v0.0
L. J. Duarte, XXXXXXX, P. L. A. Popelier 

Library with function to submit job to Gaussian and get properties values from output_files
GAUSSIAN version: G09 / G16
Check for updates at github.com/ljduarte
For details about the method, please see XXXXXXX

Please, report bugs and issues to leo.j.duarte@hotmail.com.br
coded by L. J. Duarte
"""

"""
###########################################################################################################
FUNCTION: get_atom_list_wfn_g09
          get atomic labels from g09 wfn file

INPUT: wfn_file
    wfn_file = Any wfn file of the desired PES
           
OUTPUT: atom_list
    list of each atom label for all atoms in molecule 
    
ERROR:
    "Atomic labels not found" : Atom list does not exist in wfn_file  
###########################################################################################################
"""

def get_atom_list_wfn_g09(wfn_file):
    #INTERNAL VARIABLES:
    atom_list = []
    
    #OPEN FILE:
    file = open(wfn_file, "r")
    lines = file.readlines() #Convert file into a array of lines
    file.close() #Close file
    
    #ERRORS:
    if "(CENTRE " not in lines[2]:
        raise ValueError("Atomic labels not found")  #Checks if atomic list exist inside file
        
    #GET ATOM LIST:
    for i in range(len(lines)):
        if "(CENTRE " in lines[i]:
            split_line = lines[i].split()
            atom_list.append(split_line[0].lower() + str(split_line[1])) # uppercase to lowercase
    
    return atom_list

"""
###########################################################################################################
FUNCTION: get_atom_list_wfx_g09
          get atomic labels from g09 wfn file

INPUT: wfn_file
    wfx_file = Any double wavefunction file of the desired PES
           
OUTPUT: atom_list
    list of each atom label for all atoms in molecule 
    
ERROR:
    "Atomic labels not found" : Atom list does not exist in wfn_file  
###########################################################################################################
"""

def get_atom_list_wfx_g09(wfx_file):
    #INTERNAL VARIABLES:
    atom_list = []
    
    #OPEN FILE:
    file = open(wfx_file, "r")
    lines = file.readlines() #Convert file into a array of lines
    file.close() #Close file
    
    #ERRORS:
    if " <Nuclear Names>" not in lines[33]:
        raise ValueError("Atomic labels not found")  #Checks if atomic list exist inside file
        
    #GET ATON LINES NUMBER
    for i in range(len(lines)):
        if "<Nuclear Names>" in lines[i]:
            start_line = i+1
        if "</Nuclear Names>" in lines[i]:
            end_line = i

    #GET ATOM LIST:
    for i in range(start_line, end_line):
        split_line = lines[i].split()
        atom_list.append(split_line[0].lower())
    
    return atom_list

"""
###########################################################################################################
FUNCTION: get_atom_list_wfn_g16
          get atomic labels from g16 wfn file

INPUT: wfn_file
    wfn_file = Any wfn file of the desired PES
           
OUTPUT: atom_list
    list of each atom label for all atoms in molecule 
    
ERROR:
    "Atomic labels not found" : Atom list does not exist in wfn_file  
###########################################################################################################
"""

def get_atom_list_wfn_g16(wfn_file):
    #INTERNAL VARIABLES:
    atom_list = []
    
    #OPEN FILE:
    file = open(wfn_file, "r")
    lines = file.readlines() #Convert file into a array of lines
    file.close() #Close file
    
    #ERRORS:
    if "(CENTRE " not in lines[2]:
        raise ValueError("Atomic labels not found")  #Checks if atomic list exist inside file
        
    #GET ATOM LIST:
    for i in range(len(lines)):
        if "(CENTRE " in lines[i]:
            split_line = lines[i].split()
            atom_list.append(split_line[0].lower()) # uppercase to lowercase
    
    return atom_list

