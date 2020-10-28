"""
gaussian_utils.py v0.1
L. J. Duarte, F.Falcioni, P. L. A. Popelier

Library with function to submit job to Gaussian and get properties values from output_files
GAUSSIAN version: G09 / G16
Check for updates at https://github.com/FabioFalcioni/REG.py
For details about the method, please see XXXXXXX

Please, report bugs and issues to fabioslefou@gmail.com
coded by L. J. Duarte and F. Falcioni
"""

import re

def get_atom_list_wfn_g09(wfn_file):
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


def get_atom_list_wfx_g09(wfx_file):
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


def get_atom_list_wfn_g16(wfn_file):
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

def get_control_coordinates_IRC_g16(output_file):
    '''
    ###########################################################################################################
    FUNCTION: get_control-coordinates
        get control coordinates from g16 output file
    INPUT: output_file
        output_file = Gaussian16 output file for Intrinsic Reaction Coordinate (IRC) scan.

    OUTPUT: coordinates
        list of the control coordinate of the IRC scan

    ERROR:
        "Control coordinates not found. Please, check that you have the g16 IRC output in the running folder"
    ###########################################################################################################

    '''
    #INTERNAL VARIABLES
    coordinates = []
    start = 0
    end = 0
    found = False

    #WORKING IN THE FILE
    with open(output_file, 'r') as f:
        lines = f.readlines()
        # GET THE FIRST/LAST COORDINATE POSITION
        for i in range(len(lines)):
            if "Summary of reaction path following" in lines[i]:
                start = i + 3
                found = True

        # ERRORS
        if found is False:
            raise ValueError("Control coordinates not found. Please, check that you have the g16 IRC output in the running folder")

        for j in range(start, len(lines)):
            if "---" in lines[j]:
                end = j
                break

        # GET COORDINATES LIST
        for line in lines[start: end]:
            coordinates.append(float(line.split()[2]))

    return coordinates

def get_mp2_corr(file_list):
    energies = []
    for f in file_list:
        temp = open(f, 'r').read()
        temp = re.sub(r"[\n\t\s]*", "", temp)
        hf = float(temp[temp.find('\HF=-')+4:temp.find('\HF=-')+15])
        mp2 = float(temp[temp.find('\MP2=-')+5:temp.find('\MP2=-')+16])
        energies.append(mp2-hf)
        
    return energies

