#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
aimall_utils.py v0.0
L. J. Duarte, A. Cador, P. L. A. Popelier 

Library with functions to submit jobs to AIMAll and get properties values from output
AIMAll version: 17.11.14
Check for updates at github.com/ljduarte
For details about the method, please see XXXXXXX

Please, report bugs and issues to leo.j.duarte@hotmail.com.br or ael.cador@univ-rouen.fr
coded by L. J. Duarte and A. Cador
"""

import numpy as np

"""
###########################################################################################################

PROCEDURE LIST

###########################################################################################################

model_check
          check if the model used in AIMAll is correct
get_atom_list
          get atomic labels from wfn file
get_cc_list_wfx
          get control coordinate values from wfx file
get_list_wfx
          get control coordinate values from wfx file and paths to folders
      get_list_wfx_dyn
          get control coordinate values from wfx file (tailored for explicit solvent study)
get_atom_list_wfx
          get atomic labels from wfx file
get_atom_coords_wfx
          get atomic coordinates from wfx file
get_aimall_wfn_energies
          get all wfn energies from the wavefunction files (wfn)
get_aimall_wfx_energies
          get all wfx energies from the wavefunction files (wfx)
intra_property_from_int_file
          get IQA intra-atomic properties from int files
      intra_property_from_int_file_dyn
          get IQA intra-atomic properties from int files (tailored for explicit solvent study)
inter_property_from_int_file
          get IQA interatomic properties from int files
      inter_property_from_int_file_dyn
          get IQA interatomic properties from int files (tailored for explicit solvent study)
disp_property
          get dispersion energies from the D3 calculations in Gaussian outputs

###########################################################################################################
"""

"""
###########################################################################################################
FUNCTION: model_check
          check if the model used in AIMAll is correct

INPUT: wfx_files, model, folders, atom_list
    wfx_files = List of the wfx files of the desired PES
    model = model used in AIMAll calculations (e.g. RB3LYP)
    folders = path to _atomicfiles folders
    atom_list = list of atomic labels e.g.: [n1, c2, h3, ...]
    verbose = Boolean to ask the procedure to say everything or be quiet      
           
OUTPUT: 
    "Correct model used"
    
ERROR:
    "No model found in .wfx" : Model tag does not exist in wfn_file  
    "Incorrect model found in .wfx" : Model tag exists in wfn_file but indicates the wrong model
    
Added by Aël Cador (13-09-21)
###########################################################################################################
"""

def model_check(wfx_files, model, folders, atom_list, SYS, verbose=False):
    print("Checking .wfx files of ",SYS," for Model tag...")
    for wfx_file in wfx_files:
        #OPEN FILE:
        file = open(wfx_file, "r")
        lines = file.readlines() #Convert file into a array of lines
        file.close() #Close file
        flag_model=False
        
        for i in range(len(lines)):
            if "<Model>" in lines[i]:
                flag_model = True
                if model not in lines[i+1]:
                    raise ValueError("No model found in",wfx_file)  #Checks which model is used in the file
                else:
                    if verbose == True:
                        print("Correct model found in",wfx_file)
                    break
                
        if flag_model == False:
            raise ValueError("No model found in",wfx_file)  #Checks if <Model> tag exists inside file
    
    print("Correct model used in all .wfx files!")
    print("Checking .int files of ",SYS," for model...")
    
    for folder in folders:
        for atom in atom_list:
            file = open(folder + "/" + atom + ".int", "r")
            intfilename = folder + "/" + atom + ".int"
            lines = file.readlines()
            file.close()
            flag_model=False
            
            for i in range(len(lines)):
                if "Model:  "+model[0]+"estricted "+model[1:] in lines[i]: #should work even for Unrestricted models
                    flag_model = True
                    break
                
            if flag_model == False:
                raise ValueError("No model found in ",intfilename)  #Checks if <Model> tag exists inside file
            else:
                if verbose == True:
                    print("Correct model found in ",intfilename)
    
    for path in folders:
        for i in range(len(atom_list)):
            atom1 = atom_list[i]
            for j in range(i+1, len(atom_list)):
                atom2 = atom_list[j]
                file = open(path +"/" + atom1 + "_" + atom2 + ".int", "r")
                intfilename = path +"/" + atom1 + "_" + atom2 + ".int"
                lines = file.readlines()
                file.close()
                flag_model=False
            
            for i in range(len(lines)):
                if "Model:  "+model[0]+"estricted "+model[1:] in lines[i]: #should work even for Unrestricted models
                    flag_model = True
                    break
                
            if flag_model == False:
                raise ValueError("No model found in ",intfilename)  #Checks if <Model> tag exists inside file
            else:
                if verbose == True:
                    print("Correct model found in ",intfilename)
                
    print("Correct model used in all .int files!")
"""
###########################################################################################################
FUNCTION: get_atom_list
          get atomic labels from wfn file

INPUT: wfn_file, verbose
    wfn_file = Any wfn file of the desired PES
    verbose = Boolean to ask the procedure to say everything or be quiet
       
OUTPUT: atom_list
    list of each atom label for all atoms in molecule 
    
ERROR:
    "Atomic labels not found" : Atom list does not exist in wfn_file  
###########################################################################################################
"""

def get_atom_list(wfn_file, verbose=False):
    print("Reading atom list in ",wfn_file,"...")
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
    
    if verbose==True:
        print("Atom list:",atom_list)
    print("Atom list successfully found!")
    return atom_list

"""
###########################################################################################################
FUNCTION: get_cc_list_wfx
          get control coordinate values from wfx file
          (this requires that the title of the Gaussian input file contains the control coordinate,
           so that it is written to the .wfx file as well)
  
INPUT: wfx_files, verbose
    wfx_files = List of the .wfx files of the desired PES
    verbose = Boolean to ask the procedure to say everything or be quiet
    
OUTPUT: cc_list
    list of the control coordinates of the PES
    
ERROR:
    "Control coordinate not found in file_name" : control coordinate not present in wfx_files
Added by A. Cador (April 2021)
###########################################################################################################
"""

def get_cc_list_wfx(wfx_files, verbose=False):
    print("Creating control coordinate list...")
    #INTERNAL VARIABLES:
    cc_list = []
    
    for wfx_file in wfx_files:
        if verbose==True:
            print("Reading file:",wfx_file)
        #OPEN FILE:
        file = open(wfx_file, "r")
        lines = file.readlines() #Convert file into a array of lines
        file.close() #Close file
        
        #ERRORS:
        if "coord" not in lines[1]:
            raise ValueError("Control coordinate not found in",wfx_file)  #Checks if coordinate value exists inside file
        
        for i in range(len(lines)):
            if "coord" in lines[i]:
                if verbose==True:
                    print("Control coordinate found in file:",wfx_file)
                split_line = lines[i].split()
                cc_list.append(float(split_line[4]))
    print("Control coordinate list successfully created!")
    if verbose==True:
        print(cc_list)        
    return(cc_list)

"""
###########################################################################################################
FUNCTION: get_list_wfx
          get control coordinate values from wfx file and paths to folders
  
INPUT: path, sys, way, F, R, inv, verbose=False
    path    = path to input file folder
    sys     = system considered
    F       = point number list for Forward path
    R       = point number list for Reverse path
    inv     = Boolean to indicate if the F/R reaction paths need to be inverted or not to recover the actual Reactant to Product pathway
                (since Gaussian can sometimes invert the Forward and Reverse paths).
                True = invert
                False = don't invert'
    verbose = Boolean to ask the procedure to say everything or be quiet
    
OUTPUT: cc_a, folders, folders_disp, wfx_files
    cc_a = list of the control coordinates of the PES 
    folders = list of the paths to the folders of atomic terms
    folders_disp = path to the dispersion term folders
    wfx_files = list of the paths to the .wfx files    
    
ERROR:
    "Control coordinate not found in file_name" : control coordinate not present in wfx_files
Added by A. Cador (April 2021)
###########################################################################################################
"""

def get_list_wfx(path, sys, way, F, R, inv, verbose=False):
    print("Creating file paths and reading control coordinates...")
    
    if inv==False:
        Rr=[]
        Rr=R.copy()
        Rr.reverse()
        #print(R,Rr)
        wfx_filesR = [path + '/C4-'+sys+'-P'+way+'_R_Pt' + str(i) + '.wfx' for i in Rr]
        wfx_filesF = [path + '/C4-'+sys+'-P'+way+'_F_Pt' + str(i) + '.wfx' for i in F]
        wfx_files = wfx_filesR + wfx_filesF
        out_filesR = [path + '/C4-'+sys+'-P'+way+'_R_SP' + str(i) + '.out' for i in Rr] #added by AC 20-04-21 to get .gjf files
        out_filesF = [path + '/C4-'+sys+'-P'+way+'_F_SP' + str(i) + '.out' for i in F] #added by AC 20-04-21 to get .gjf files
        out_files = out_filesR + out_filesF
        cc_a=get_cc_list_wfx(wfx_files)
        foldersR = [path + '/C4-'+sys+'-P'+way+'_R_Pt' + str(i) + '_atomicfiles' for i in Rr]
        foldersF = [path + '/C4-'+sys+'-P'+way+'_F_Pt' + str(i) + '_atomicfiles' for i in F]
        folders = foldersR + foldersF
    
    if inv==True:
        Fr=[]
        Fr=F.copy()
        Fr.reverse()
        #print(F,Fr)
        wfx_filesR = [path + '/C4-'+sys+'-P'+way+'_R_Pt' + str(i) + '.wfx' for i in R]
        wfx_filesF = [path + '/C4-'+sys+'-P'+way+'_F_Pt' + str(i) + '.wfx' for i in Fr]
        wfx_files = wfx_filesF + wfx_filesR
        out_filesR = [path + '/C4-'+sys+'-P'+way+'_R_SP' + str(i) + '.out' for i in R] #added by AC 20-04-21 to get .gjf files
        out_filesF = [path + '/C4-'+sys+'-P'+way+'_F_SP' + str(i) + '.out' for i in Fr] #added by AC 20-04-21 to get .gjf files
        out_files = out_filesF + out_filesR
        cc_a=get_cc_list_wfx(wfx_files)
        for i in range(len(cc_a)):
            cc_a[i]=-cc_a[i]
        foldersR = [path + '/C4-'+sys+'-P'+way+'_R_Pt' + str(i) + '_atomicfiles' for i in R]
        foldersF = [path + '/C4-'+sys+'-P'+way+'_F_Pt' + str(i) + '_atomicfiles' for i in Fr]
        folders = foldersF + foldersR
    
    folders_disp = path + '/ResDisp_'+sys+'-SP'+way #this file can just be a premade file where you "grep"-ed the D3 energies from the Gaussian output.
    #Don't forget to order the lines correctly!
    
    print("File paths successfully created!")
    if verbose==True:
        print("IQA file folders:",folders)
        print("Dispersion file path:",folders_disp)
        print(".wfx file list:",wfx_files)
        print("Gaussian single-point calculation file list:",out_files)
    
    return [cc_a, folders, folders_disp, wfx_files, out_files]

def get_list_wfx_dyn(path, sys, F, R, inv, nbo, thres="Low", verbose=False):
    print("Creating file paths and reading control coordinates...")
    
    if sys == '6+18DMSO':
        file_prefix = '/6_6-31+Gd'
    if sys == '6_TS-B3LYP-D3_6-31+Gd':
        file_prefix = '/6_TS-B3LYP_6-31+Gd'
    if nbo == True and thres == "Low":
        file_prefix_log = '/NBO/LowThres' + file_prefix
    elif nbo == True and thres == "High":
        file_prefix_log = '/NBO/HighThres' + file_prefix
    else:
        file_prefix_log = file_prefix
    file_prefix_wfx = '/WFX' + file_prefix
    if inv==False:
        Rr=[]
        Rr=R.copy()
        Rr.reverse()
        #print(R,Rr)
        wfx_filesR = [path + file_prefix_wfx + '_R_Pt' + str(i) + '.wfx' for i in Rr]
        wfx_filesF = [path + file_prefix_wfx + '_F_Pt' + str(i) + '.wfx' for i in F]
        wfx_files = wfx_filesR + wfx_filesF
        out_filesR = [path + file_prefix_log + '_R_SP' + str(i) + '.out' for i in Rr] #added by AC 20-04-21 to get .gjf files
        out_filesF = [path + file_prefix_log + '_F_SP' + str(i) + '.out' for i in F] #added by AC 20-04-21 to get .gjf files
        out_files = out_filesR + out_filesF
        cc_a=get_cc_list_wfx(wfx_files)
        foldersR = [path + file_prefix + '_R_Pt' + str(i) + '_atomicfiles' for i in Rr]
        foldersF = [path + file_prefix + '_F_Pt' + str(i) + '_atomicfiles' for i in F]
        folders = foldersR + foldersF
        
    if inv==True:
        Fr=[]
        Fr=F.copy()
        Fr.reverse()
        #print(F,Fr)
        wfx_filesR = [path + file_prefix_wfx + '_R_Pt' + str(i) + '.wfx' for i in R]
        wfx_filesF = [path + file_prefix_wfx + '_F_Pt' + str(i) + '.wfx' for i in Fr]
        wfx_files = wfx_filesF + wfx_filesR
        out_filesR = [path + file_prefix_log + '_R_SP' + str(i) + '.out' for i in R] #added by AC 20-04-21 to get .gjf files
        out_filesF = [path + file_prefix_log + '_F_SP' + str(i) + '.out' for i in Fr] #added by AC 20-04-21 to get .gjf files
        out_files = out_filesF + out_filesR
        cc_a=get_cc_list_wfx(wfx_files)
        for i in range(len(cc_a)):
            cc_a[i]=-cc_a[i]
        foldersR = [path + file_prefix + '_R_Pt' + str(i) + '_atomicfiles' for i in R]
        foldersF = [path + file_prefix + '_F_Pt' + str(i) + '_atomicfiles' for i in Fr]
        folders = foldersF + foldersR
    
    folders_disp = path + '/ResDisp_'+sys #this file can just be a premade file where you "grep"-ed the D3 energies from the Gaussian output.
    #Don't forget to order the lines correctly!
    
    print("File paths successfully created!")
    if verbose==True:
        print("IQA file folders:",folders)
        print("Dispersion file path:",folders_disp)
        print(".wfx file list:",wfx_files)
        print("Gaussian single-point calculation file list:",out_files)
        
    return [cc_a, folders, folders_disp, wfx_files, out_files]

"""
###########################################################################################################
FUNCTION: get_atom_list_wfx
          get atomic labels from wfx file

INPUT: wfx_file, verbose=False
    wfx_file = Any .wfx file of the desired PES
           
OUTPUT: atom_list
    list of each atom label for all atoms in molecule 
    verbose = Boolean to ask the procedure to say everything or be quiet
    
ERROR:
    "Atomic labels not found in file_name" : Atom list does not exist in wfx_file
    
Added by A. Cador (February 2020)   
###########################################################################################################
"""

def get_atom_list_wfx(wfx_file, verbose=False):
    print("Reading atom list in ",wfx_file,"...")
    #INTERNAL VARIABLES:
    atom_list = []
    
    #OPEN FILE:
    file = open(wfx_file, "r")
    lines = file.readlines() #Convert file into a array of lines
    file.close() #Close file
    
    #ERRORS:
    if "<Nuclear Names>" not in lines[33]: #36 when the Model tags are added, 33 without; 37 in rare cases when the cc line overflows
        if "<Nuclear Names>" not in lines[36]:
            if "<Nuclear Names>" not in lines[37]:
                raise ValueError("Atomic labels not found in",wfx_file)  #Checks if atomic list exist inside file
        
    #GET ATOM LIST:
    for i in range(len(lines)):
        if "<Nuclear Names>" in lines[i]:
                i+=1
                while "</Nuclear Names>" not in lines[i]:
                    split_line = lines[i].split()
                    atom_list.append(split_line[0].lower()) # uppercase to lowercase
                    i+=1
   
    if verbose==True:
        print("Atom list:",atom_list)
    print("Atom list successfully found!")
    
    return atom_list

"""
###########################################################################################################
FUNCTION: get_atom_coords_wfx
          get atomic coordinates from wfx file

INPUT: wfx_file, verbose=False
    wfx_files = List of the .wfx files of the desired PES
    verbose = Boolean to ask the procedure to say everything or be quiet
           
OUTPUT: atom_coordinates
    list of lists of coordinates each atom label for all atoms in molecule
    atom → [x,y,z]
    structure of point i → [atom1,atom2,...]
    complete list → [structure1,structure2,...]
    
ERROR:
    "Atomic labels not found in file_name" : Atom list does not exist in wfx_file
    
Added by A. Cador  20-04-21  
###########################################################################################################
"""

def get_atom_coord_wfx(wfx_files, verbose=False):
    print("Reading atom coordinates...")
    #INTERNAL VARIABLES:
    atom_coord_list = []
    
    #OPEN FILE:
    for wfx_file in wfx_files:
        if verbose==True:
            print("Reading file:",wfx_file)
        file = open(wfx_file, "r")
        lines = file.readlines() #Convert file into a array of lines
        file.close() #Close file
        
        #ERRORS:
        if "<Nuclear Cartesian Coordinates>" not in lines[75]: #78 when the Model tags are added, 75 without
            if "<Nuclear Cartesian Coordinates>" not in lines[78]:
                raise ValueError("Atomic coordinates not found in",wfx_file)  #Checks if atomic list exist inside file
            
        #GET ATOM LIST:
        for i in range(len(lines)):        
            if "<Nuclear Cartesian Coordinates>" in lines[i]: #go until the coordinate section is found
                if verbose==True:
                    print("Atom coordinates found in file:",wfx_file)
                i+=1
                atom_coord_list_i=[]
                while "</Nuclear Cartesian Coordinates>" not in lines[i]: #repeat until the end of the coordinate section is found
                    coords = lines[i].split()
                    for j in range(3): #(3 coordinates, X, Y, Z)
                        coords[j]=float(coords[j]) #convert the string to float
                    atom_coord_list_i.append(coords)
                    #print(coords)
                    i+=1
                atom_coord_list.append(atom_coord_list_i)
                if verbose==True:
                    print("Atom coordinates in file:",wfx_file)
                    print(atom_coord_list_i)
    
    atom_coord_list=np.array(atom_coord_list)
    print("Atom coordinates succesfully found!")
    
    return atom_coord_list


"""
###########################################################################################################
FUNCTION: get_aimall_wfn_energies
          get all wfn energies from the wavefunction files (wfn)

INPUT: wfn_files, verbose=False
    wfn_files = list of all wfn files of the PES
    verbose = Boolean to ask the procedure to say everything or be quiet
           
OUTPUT: wfn_energy
    wfn_energy = list of energies for each PES point 
    
ERROR:
    "Energy values not found in file : file_name" : No energy values in file_name  
###########################################################################################################
"""    
      
def get_aimall_wfn_energies(wfn_files, verbose=False):
    print("Reading Gaussian energies from .wfn files...")
    #INTERNAL VARIABLES:
    wfn_energy = [] #List of wfn files

    #READ FILES 
    for path in wfn_files:
        if verbose==True:
            print("Reading file",path)
        file = open(path, "r")
        lines = file.readlines()
        #ERRORS:
        if "TOTAL ENERGY " not in lines[-1]: #Checks if there is an energy value at the end of the .wfn file.
            raise ValueError("Energy values not found in file: ", path)      
        wfn_energy.append(float(lines[-1].split()[3]))
        if verbose==True:
            print("Energy from file",path,":",wfn_energy[-1])
        file.close()
    print("Gaussian energy list successfully created!")
    return wfn_energy

"""
###########################################################################################################
FUNCTION: get_aimall_wfx_energies
          get all wfx energies from the wavefunction files (wfx)

INPUT: wfx_files, verbose=False
    wfx_files = list of all wfx files of the PES
    verbose = Boolean to ask the procedure to say everything or be quiet
       
OUTPUT: wfn_energy
    wfn_energy = list of energies for each PES point 
    
ERROR:
    "Energy values not found in file : file_name" : No energy values in file_name

Added by A. Cador (February 2020)   
###########################################################################################################
""" 

def get_aimall_wfx_energies(wfx_files, verbose=False):
    print("Reading Gaussian energies from .wfx files...")
    #INTERNAL VARIABLES:
    wfx_energy = [] #List of wfn files

    #READ FILES 
    for path in wfx_files:
        if verbose==True:
            print("Reading file",path)
        file = open(path, "r")
        lines = file.readlines()
        #ERRORS:
        if "<Energy = T + Vne + Vee + Vnn>" not in lines[-6]: #Checks if there is an energy value at the end of the .wfx file.
            raise ValueError( "Energy values not found in file: ", path)      
        wfx_energy.append(float(lines[-5].split()[0]))
        if verbose==True:
            print("Energy from file",path,":",wfx_energy[-1])
        file.close()
    print("Gaussian energy list successfully created!")
    return wfx_energy

"""
###########################################################################################################
FUNCTION: intra_property_from_int_file
          get IQA intra-atomic properties from int files

INPUT: folders, prop, atom_list, verbose=False
    folders = path to _atomicfiles folders
    prop = list of IQA for each atoms e.g.: "['T(A)', 'Vee(A,A)', 'Vne(A,A)']" 
    atom_list = list of atomic labels e.g.: [n1, c2, h3, ...]       
    verbose = Boolean to ask the procedure to say everything or be quiet
    
OUTPUT: [intra_properties, contributions_list]
    intra_properties = array of array containing the IQA values for each atom for each geometry
    contributions_list = list of contributions  organized in the same order as in intra_properties
    
ERROR:
    "File is empty or does not exist: file_name" no terms found in file or file doesn't exist    
###########################################################################################################
""" 

def intra_property_from_int_file(folders, prop, atom_list, verbose=False):
    print("Reading Intra terms from .int files...")
    #INTERNAL VARIABLES:
    temp1 = [] #Temporary array
    temp2 = [] #Temporary array
    temp3 = [] #Temporary array
    intra_properties = [] #Output
    contributions_list= [] #Output
    
    #READ PROPERTIES FROM .INT FILES
    for folder in folders:
        if verbose==True:
            print("Looking in folder",folder)
        for atom in atom_list:
            if verbose==True:
                print("Looking in file",folder+"/"+atom+".int")
            file = open(folder + "/" + atom + ".int", "r")
            lines = file.readlines()
            
            file.close()
            for i in lines:
                if 'IQA Energy Components (see "2EDM Note")' in i:
                    start = lines.index(i)
                elif '2EDM Note:' in i:
                    end = lines.index(i)
            if end >= len(lines) : #Checks the .int file.
                raise ValueError( "File is empty or does not exist: " + folder +"/" + atom + ".int") #Corrected from atom1 + "-" + atom2
            lines = [lines[i] for i in range(start+1, end)]
            for term in prop:
                for i in lines:
                    if (term + '           ') in i:
                        temp1.append(float(i.split()[-1]))
                        if verbose==True:
                            print("Intra term in file",folder + "/" + atom + ".int :",temp1[-1])
                                  
    #ORGANIZE ARRAY ORDER                
    for j in range(len(prop)):
        for i in range(j, len(temp1), len(prop)):
            temp2.append(temp1[i])
    for j in range(len(atom_list)):
         temp3.append([temp2[i] for i in range(j, len(temp2), len(atom_list))]) 
        
    start = 0 
    for j in range(len(prop)):
        for atom_prop in temp3:
            intra_properties.append([atom_prop[i] for i in range(start, (j+1)*(len(folders)))])
        start = (j+1)*len(folders)
    
    #CREATE CONTRIBUTIONS LIST ARRAY:   
    for a in prop:
        for b in atom_list:
            contributions_list.append(a + '-' + b)
    
    print("Intra terms successfully found!")
    return intra_properties, contributions_list

def intra_property_from_int_file_dyn(folders, prop, atom_list, verbose=False):
    print("Reading Intra terms from .int files...")
    #INTERNAL VARIABLES:
    temp1 = [] #Temporary array
    temp2 = [] #Temporary array
    temp3 = [] #Temporary array
    intra_properties = [] #Output
    contributions_list= [] #Output
    #ind_rem = ["s54","o55","s56","c66","c81","c83","c91","h173","h174","h175","h197","h198","h199"]
    #avoid = False
    count_av = 0
    #READ PROPERTIES FROM .INT FILES
    for folder in folders:
        if verbose==True:
            print("Looking in folder",folder)
        #count_av=0
        for atom in atom_list:
            #for i in ind_rem:
                #if atom == i:
                    #avoid = True
                    #break
            #if avoid == True:
                #count_av += 1
            #else:
                if verbose==True:
                    print("Looking in file",folder+"/"+atom+".int")
                file = open(folder + "/" + atom + ".int", "r")
                lines = file.readlines()
                
                file.close()
                for i in lines:
                    if 'IQA Energy Components (see "2EDM Note")' in i:
                        start = lines.index(i)
                    elif '2EDM Note:' in i:
                        end = lines.index(i)
                if end >= len(lines) : #Checks the .int file.
                    raise ValueError( "File is empty or does not exist: " + folder +"/" + atom + ".int") #Corrected from atom1 + "-" + atom2
                lines = [lines[i] for i in range(start+1, end)]
                for term in prop:
                    for i in lines:
                        if (term + '           ') in i:
                            temp1.append(float(i.split()[-1]))
                            if verbose==True:
                                print("Intra term in file",folder + "/" + atom + ".int :",temp1[-1])
            #avoid = False

    if verbose==True:
        print("IQA_INTRA_TEMP1",temp1[0],len(temp1))                              
    #ORGANIZE ARRAY ORDER                
    for j in range(len(prop)):
        for i in range(j, len(temp1), len(prop)):
            temp2.append(temp1[i])
    if verbose==True:
        print("IQA_INTRA_TEMP2",temp2[0],len(temp2))
    for j in range(len(atom_list)-count_av):
         temp3.append([temp2[i] for i in range(j, len(temp2), (len(atom_list)-count_av))]) 
    if verbose==True:
        print("IQA_INTRA_TEMP3",temp3[0],len(temp3),len(temp3[0]))    
    start = 0 
    for j in range(len(prop)):
        for atom_prop in temp3:
            intra_properties.append([atom_prop[i] for i in range(start, (j+1)*(len(folders)))])
        start = (j+1)*(len(folders))
    
    if verbose==True:
        print("IQA_INTRA_FINAL",intra_properties[0],len(intra_properties),len(intra_properties[0]))
    
    #CREATE CONTRIBUTIONS LIST ARRAY:   
    avoid = False
    for a in prop:
        for b in atom_list:
            #for i in ind_rem:
                #if b == i:
                    #avoid = True
            if avoid == False: 
                contributions_list.append(a + '-' + b)
            avoid = False
    if verbose==True:
        print("IQA_CONT_FINAL",contributions_list[0],len(contributions_list))   
    print("Intra terms successfully found!")
    
    return intra_properties, contributions_list

"""
###########################################################################################################
FUNCTION: inter_property_from_int_file
          get IQA interatomic properties from int files

INPUT: folders, prop, atom_list, verbose=False
    folders = path to _atomicfiles folders
    prop = list of IQA for each atoms e.g.: "['T(A)', 'Vee(A,A)', 'Vne(A,A)']" 
    atom_list = list of atomic lables e.g.: [n1, c2, h3, ...]       
    verbose = Boolean to ask the procedure to say everything or be quiet

OUTPUT: [intra_properties, contributions_list]
    intra_properties = array of array containing the IQA values for eacha atom for each geometry
    contributions_list = list of contributions  organized in the same order as in intra_properties
    
ERROR:
    "File is empty or does not exist: file_name" no terms found in file or file doesn't exist
###########################################################################################################
""" 

def inter_property_from_int_file(folders, prop, atom_list, verbose=False):
    print("Reading Inter terms from .int files...")
    #INTERNAL VARIABLES:
    temp1 = [] #Temporary array
    temp2 = [] #Temporary array
    temp3 = [] #Temporary array
    inter_properties = [] #Output
    contributions_list= [] #Output
    
    for folder in folders:
        if verbose==True:
            print("Looking in folder",folder)
        for i in range(len(atom_list)):
            atom1 = atom_list[i]
            for j in range(i+1, len(atom_list)):
                atom2 = atom_list[j]
                if verbose==True:
                    print("Looking in file",folder + "/" + atom1 + "_" + atom2 + ".int")
                file = open(folder + "/" + atom1 + "_" + atom2 + ".int", "r")
                lines = file.readlines()
                file.close()
                for i in lines:
                    if ' Energy Components (See "2EDM Note"):' in i:
                        start = lines.index(i)
                    elif '2EDM Note:' in i:
                        end = lines.index(i)
                if end >= len(lines) : #Checks the .int file.
                    raise ValueError("File is empty or does not exist: " + folder +"/" + atom1 + "_" + atom2 + ".int") 
                lines = [lines[i] for i in range(start+1, end)]
                for term in prop:
                    for i in lines:
                        if (term + '  ') in i:
                                temp1.append(float(i.split()[-1]))
                                if verbose==True:
                                    print("Inter term in file",folder + "/" + atom1 + "_" + atom2 + ".int",temp1[-1])
    #ORGANIZE ARRAY ORDER                
    for j in range(len(prop)):
        for i in range(j, len(temp1), len(prop)):
            temp2.append(temp1[i])
    for j in range(int(len(atom_list)*(len(atom_list)-1)/2)):
         temp3.append([temp2[i] for i in range(j, len(temp2), int(len(atom_list)*(len(atom_list)-1)/2))]) 
    start = 0 
    for j in range(len(prop)):
        for atom_prop in temp3:
            inter_properties.append([atom_prop[i] for i in range(start, (j+1)*len(folders))])
        start = (j+1)*len(folders)
        #CREATE CONTRIBUTIONS LIST ARRAY:   
    for a in prop:
        for i in range(len(atom_list)):
            for j in range(i+1, len(atom_list)):        
                contributions_list.append(a + '-' + atom_list[i] + '_' + atom_list[j])   
       
    print("Inter terms successfully found!")
    
    return inter_properties, contributions_list

def inter_property_from_int_file_dyn(folders, prop, atom_list, verbose=False):
    print("Reading Inter terms from .int files...")
    #INTERNAL VARIABLES:
    temp1 = [] #Temporary array
    temp2 = [] #Temporary array
    temp3 = [] #Temporary array
    inter_properties = [] #Output
    contributions_list= [] #Output
    #ind_rem = ["s54","o55","s56","c66","c81","c83","c91","h173","h174","h175","h197","h198","h199"]
    #avoid = False
    #count_av = 0
    solu = 19
    solv_at = 10
    solv_n = 18
    
    for folder in folders:
        if verbose==True:
            print("Looking in folder",folder)
        #count_av=0
        for i in range(0,19): #solute with all
            atom1 = atom_list[i]
            for j in range(i+1, len(atom_list)):
                atom2 = atom_list[j]
                #for k in ind_rem:
                #    if atom2 == k:
                #        avoid = True
                #        break
                #if avoid == True:
                #    count_av += 1  
                #else:
                if verbose==True:
                    print("Looking in file",folder + "/" + atom1 + "_" + atom2 + ".int")
                file = open(folder +"/" + atom1 + "_" + atom2 + ".int", "r")
                lines = file.readlines()
                file.close()
                for i in lines:
                    if ' Energy Components (See "2EDM Note"):' in i:
                        start = lines.index(i)
                    elif '2EDM Note:' in i:
                        end = lines.index(i)
                if end >= len(lines) : #Checks the .int file.
                    raise ValueError( "File is empty or does not exist: " + folder +"/" + atom1 + "_" + atom2 + ".int") 
                lines = [lines[i] for i in range(start+1, end)]
                for term in prop:
                    for i in lines:
                        if (term + '  ') in i:
                            temp1.append(float(i.split()[-1]))
                            if verbose==True:
                                print("Inter term in file",folder + "/" + atom1 + "_" + atom2 + ".int",temp1[-1])
                #avoid = False
        solv_list=[[20,21,66,67,122,123,124,125,126,127],[22,23,81,89,167,168,169,191,192,193],[24,25,56,57,92,93,94,95,96,97],
         [26,27,68,69,128,129,130,131,132,133],[28,29,82,90,170,171,172,194,195,196],[30,31,60,61,104,105,106,107,108,109],
         [32,33,84,87,176,177,178,185,186,187],[34,35,75,76,149,150,151,152,153,154],[36,37,73,74,143,144,145,146,147,148],
         [38,39,70,71,134,135,136,137,138,139],[40,41,77,78,155,156,157,158,159,160],[42,43,62,63,110,111,112,113,114,115],
         [44,45,79,80,161,162,163,164,165,166],[46,47,64,65,116,117,118,119,120,121],[54,55,83,91,173,174,175,197,198,199],
         [48,49,72,88,140,141,142,188,189,190],[50,51,58,59,98,99,100,101,102,103],[52,53,85,86,179,180,181,182,183,184]]
 
        for h in range(len(solv_list)): #intra-solvent terms
            for i in range(len(solv_list[h])):
                atom1 = atom_list[solv_list[h][i]-1]
                for j in range(i+1,len(solv_list[h])):
                    atom2 = atom_list[solv_list[h][j]-1]
                    file = open(folder +"/" + atom1 + "_" + atom2 + ".int", "r")
                    lines = file.readlines()
                    file.close()
                    for i in lines:
                        if ' Energy Components (See "2EDM Note"):' in i:
                            start = lines.index(i)
                        elif '2EDM Note:' in i:
                            end = lines.index(i)
                    if end >= len(lines) : #Checks the .int file.
                        raise ValueError( "File is empty or does not exist: " + folder +"/" + atom1 + "_" + atom2 + ".int") 
                    lines = [lines[i] for i in range(start+1, end)]
                    for term in prop:
                        for i in lines:
                            if (term + '  ') in i:
                                temp1.append(float(i.split()[-1]))
                                if verbose==True:
                                    print("Inter term in file",folder + "/" + atom1 + "_" + atom2 + ".int",temp1[-1])
    
    if verbose==True:
        print("INTER_TEMP1", temp1[0], len(temp1))
    
    #ORGANIZE ARRAY ORDER                
    for j in range(len(prop)):
        for i in range(j, len(temp1), len(prop)):
            temp2.append(temp1[i])
    if verbose==True:
        print("INTER_TEMP2", temp2[0], len(temp2))
    int_tot = int((solu*(solu-1)/2) + (solv_at*(solv_at-1)/2)*solv_n + solu*solv_at*solv_n) #4401/285
    for j in range(int_tot):
         temp3.append([temp2[i] for i in range(j, len(temp2), int_tot)]) 
    if verbose==True:
        print("INTER_TEMP3", temp3[0], len(temp3), len(temp3[0]))
    start = 0 
    for j in range(len(prop)):
        for atom_prop in temp3:
            inter_properties.append([atom_prop[i] for i in range(start, (j+1)*(len(folders)))])
        start = (j+1)*len(folders)
    if verbose==True:
        print("INTER_PROP_FINAL", inter_properties[0], len(inter_properties), len(inter_properties[0]))
    
    #CREATE CONTRIBUTIONS LIST ARRAY:   
    #avoid = False
    for a in prop:
        for i in range(0,19):
            for j in range(i+1, len(atom_list)):
                #for k in ind_rem:
                #    if atom_list[j] == k:
                #        avoid = True
                #        break
                #if avoid == False:
                contributions_list.append(a + '-' + atom_list[i] + '_' + atom_list[j])
                #avoid = False
        for h in range(len(solv_list)): #intra-solvent terms
            for i in range(len(solv_list[h])):
                for j in range(i+1,len(solv_list[h])):  
                    contributions_list.append(a + '-' + atom_list[solv_list[h][i]-1] + '_' + atom_list[solv_list[h][j]-1])
    if verbose==True:
        print("INTER_CONT_FINAL", contributions_list[0], len(contributions_list))
    
    print("Inter terms successfully found!")
    
    return inter_properties, contributions_list


"""
###########################################################################################################
FUNCTION: disp_property
          get dispersion energies from the D3 calculations in Gaussian outputs

INPUT: folder, verbose=False
    folders = path to result file (premade file with dispersion energies for each point of the path)      
    verbose = Boolean to ask the procedure to say everything or be quiet

OUTPUT: [disp_properties, disp_header]
    intra_properties = array of arrays containing the IQA values for each atom for each geometry # 1 array is enough
    disp_header = header of the dispersion term for the final array
    
ERROR:

    Added by Aël Cador 19/06/2020 - ...    
###########################################################################################################
"""     

def disp_property(folder, verbose=False):
    print("Reading Dispersion terms from file",folder,"...")
    disp_properties = []
    temp1 = []
    disp_header = ["Dispersion"]
    
    file = open(folder, "r")
    lines = file.readlines()
    file.close()
    for i in lines:
        temp1.append(float(i.split()[-2]))
        if verbose==True:
            print("Dispersion term in file",i.split()[0],temp1[-1])
    disp_properties.append(temp1)
    print("Dispersion terms successfully found!")
    
    return disp_properties, disp_header
    