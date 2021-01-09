"""
gaussian_utils.py v0.1
L. J. Duarte, F.Falcioni, P. L. A. Popelier

Library with function to submit job to Gaussian and get properties values from output_files
GAUSSIAN version: G09 / G16
Check for updates at https://github.com/FabioFalcioni/REG.py
For details about the method, please see XXXXXXX

Please, report bugs and issues to fabio.falcioni@manchester.ac.uk
coded by L. J. Duarte and F. Falcioni
"""

import re
import numpy as np

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


def get_control_coordinates_PES_Scan(g16_output_filelist, atom_list):
    '''
    ###########################################################################################################
    FUNCTION: get_control_coordinates_PES_Scan
        get REG control coordinates from g16 single point energy output files
    INPUT:
        g16_output_filelist = Gaussian16 output file for Single Point Energy
        atom_list = list of atoms (as strings) used for the ModRedundant option in Gaussian16

    OUTPUT:
        cc = control coordinates list

    ERROR:
        "Single Point Energy File not found"
    ###########################################################################################################
    '''
    cc = []
    xyz_files = []
    for i in range(0, len(g16_output_filelist)):
        xyz_files.append(get_xyz_file(g16_output_filelist[i]))

    # Control Coordinates search for BOND movement
    for i in range(0, len(xyz_files)):
        f1 = open(xyz_files[i], 'r')
        coordinates1_list = f1.readlines()[2:]  # temporary remove the first 2 lines of xyz file
        atom1 = [float(c) for c in re.findall(r"[-+]?\d*\.\d+|\d+", coordinates1_list[atom_list[0]])]
        atom2 = [float(c) for c in re.findall(r"[-+]?\d*\.\d+|\d+", coordinates1_list[atom_list[1]])]
        x = 0
        y = 1
        z = 2
        cc.append(np.sqrt((atom2[x] - atom1[x]) ** 2 + (atom2[y] - atom1[y]) ** 2 + (atom2[z] - atom1[z]) ** 2))
    return cc

def get_xyz_file(g16_single_point_output):
    '''
    ###########################################################################################################
    FUNCTION: get_xyz_file
        get xyz coordinates from g16 single point energy output file
    INPUT: g16_single_point_output
        g16_single_point_output = Gaussian16 output file for Single Point Energy

    OUTPUT: xyz_output
        xyz_output = xyz file format

    ERROR:
        "Single Point Energy File not found"
    ###########################################################################################################
    '''
    start = 0
    end = 0
    xyz_output = str(g16_single_point_output[:-3]) + "xyz"

    openold= open(g16_single_point_output, 'r')
    opennew= open(xyz_output, 'w')
    rlines = openold.readlines()
    for i in range (len(rlines)):
            if "Standard orientation:" in rlines[i] or "Input orientation:" in rlines[i]:
                start = i
    for m in range (start + 5, len(rlines)):
        if "---" in rlines[m]:
            end = m
            break
    opennew.write("\n{}\n\n".format(str(end-start-5)))
    ## Convert to Cartesian coordinates format
    ## convert atomic number to atomic symbol
    for line in rlines[start+5 : end] :
        words = line.split()
        word1 = int(words[1])
        word3 = str(words[3])
        ## Periodic table supported
        ## Your molecule should comprise only atoms between 1-90

        if   word1 ==   1 : word1 = "H"
        elif word1 ==   2 : word1 = "He"
        elif word1 ==   3 : word1 = "Li"
        elif word1 ==   4 : word1 = "Be"
        elif word1 ==   5 : word1 = "B"
        elif word1 ==   6 : word1 = "C"
        elif word1 ==   7 : word1 = "N"
        elif word1 ==   8 : word1 = "O"
        elif word1 ==   9 : word1 = "F"
        elif word1 ==  10 : word1 = "Ne"
        elif word1 ==  11 : word1 = "Na"
        elif word1 ==  12 : word1 = "Mg"
        elif word1 ==  13 : word1 = "Al"
        elif word1 ==  14 : word1 = "Si"
        elif word1 ==  15 : word1 = "P"
        elif word1 ==  16 : word1 = "S"
        elif word1 ==  17 : word1 = "Cl"
        elif word1 ==  18 : word1 = "Ar"
        elif word1 ==  19 : word1 = "K"
        elif word1 ==  20 : word1 = "Ca"
        elif word1 ==  21 : word1 = "Sc"
        elif word1 ==  22 : word1 = "Ti"
        elif word1 ==  23 : word1 = "V"
        elif word1 ==  24 : word1 = "Cr"
        elif word1 ==  25 : word1 = "Mn"
        elif word1 ==  26 : word1 = "Fe"
        elif word1 ==  27 : word1 = "Co"
        elif word1 ==  28 : word1 = "Ni"
        elif word1 ==  29 : word1 = "Cu"
        elif word1 ==  30 : word1 = "Zn"
        elif word1 ==  31 : word1 = "Ga"
        elif word1 ==  32 : word1 = "Ge"
        elif word1 ==  33 : word1 = "As"
        elif word1 ==  34 : word1 = "Se"
        elif word1 ==  35 : word1 = "Br"
        elif word1 ==  36 : word1 = "Kr"
        elif word1 ==  37 : word1 = "Rb"
        elif word1 ==  38 : word1 = "Sr"
        elif word1 ==  39 : word1 = "Y"
        elif word1 ==  40 : word1 = "Zr"
        elif word1 ==  41 : word1 = "Nb"
        elif word1 ==  42 : word1 = "Mo"
        elif word1 ==  43 : word1 = "Tc"
        elif word1 ==  44 : word1 = "Ru"
        elif word1 ==  45 : word1 = "Rh"
        elif word1 ==  46 : word1 = "Pd"
        elif word1 ==  47 : word1 = "Ag"
        elif word1 ==  48 : word1 = "Cd"
        elif word1 ==  49 : word1 = "In"
        elif word1 ==  50 : word1 = "Sn"
        elif word1 ==  51 : word1 = "Sb"
        elif word1 ==  52 : word1 = "Te"
        elif word1 ==  53 : word1 = "I"
        elif word1 ==  54 : word1 = "Xe"
        elif word1 ==  55 : word1 = "Cs"
        elif word1 ==  56 : word1 = "Ba"
        elif word1 ==  57 : word1 = "La"
        elif word1 ==  58 : word1 = "Ce"
        elif word1 ==  59 : word1 = "Pr"
        elif word1 ==  60 : word1 = "Nd"
        elif word1 ==  61 : word1 = "Pm"
        elif word1 ==  62 : word1 = "Sm"
        elif word1 ==  63 : word1 = "Eu"
        elif word1 ==  64 : word1 = "Gd"
        elif word1 ==  65 : word1 = "Tb"
        elif word1 ==  66 : word1 = "Dy"
        elif word1 ==  67 : word1 = "Ho"
        elif word1 ==  68 : word1 = "Er"
        elif word1 ==  69 : word1 = "Tm"
        elif word1 ==  70 : word1 = "Yb"
        elif word1 ==  71 : word1 = "Lu"
        elif word1 ==  72 : word1 = "Hf"
        elif word1 ==  73 : word1 = "Ta"
        elif word1 ==  74 : word1 = "W"
        elif word1 ==  75 : word1 = "Re"
        elif word1 ==  76 : word1 = "Os"
        elif word1 ==  77 : word1 = "Ir"
        elif word1 ==  78 : word1 = "Pt"
        elif word1 ==  79 : word1 = "Au"
        elif word1 ==  80 : word1 = "Hg"
        elif word1 ==  81 : word1 = "Tl"
        elif word1 ==  82 : word1 = "Pb"
        elif word1 ==  83 : word1 = "Bi"
        elif word1 ==  84 : word1 = "Po"
        elif word1 ==  85 : word1 = "At"
        elif word1 ==  86 : word1 = "Rn"
        elif word1 ==  87 : word1 = "Fe"
        elif word1 ==  88 : word1 = "Ra"
        elif word1 ==  89 : word1 = "Ac"
        elif word1 ==  90 : word1 = "Th"

        ## copy from atom list.

        opennew.write("{}{}\n".format(word1,line[30:-1]))
    openold.close()
    opennew.close()
    return xyz_output




