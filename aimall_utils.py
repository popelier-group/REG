"""
aimall_utils.py v0.0
L. J. Duarte, XXXXXXX, P. L. A. Popelier 

Library with function to submit job to AIMAll and get properties values from output
AIMAll version: 17.11.14
Check for updates at github.com/ljduarte
For details about the method, please see XXXXXXX

Please, report bugs and issues to leo.j.duarte@hotmail.com.br
coded by L. J. Duarte
"""


def get_atom_list(wfn_file):
    """
    ###########################################################################################################
    FUNCTION: get_atom_list
              get atomic labels from wfn file

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

def get_atom_list_wfx(wfx_file):
    #INTERNAL VARIABLES:
    atom_list = []
    
    #OPEN FILE:
    file = open(wfx_file, "r")
    lines = file.readlines() #Convert file into a array of lines
    file.close() #Close file
    
    #ERRORS:
    if "<Nuclear Names>" not in lines[33]:
        raise ValueError("Atomic labels not found")  #Checks if atomic list exist inside file
        
    #GET ATOM LIST:
    for i in range(len(lines)):
        if "<Nuclear Names>" in lines[i]:
                i+=1
                while "</Nuclear Names>" not in lines[i]:
                    split_line = lines[i].split()
                    atom_list.append(split_line[0].lower()) # uppercase to lowercase
                    i+=1
    
    return atom_list


      
def get_aimall_wfn_energies(A):
    """
    ###########################################################################################################
    FUNCTION: get_aimall_wfn_energies
              get all wfn energies from the wavefunction files (wfn)

    INPUT: A
        A = list of all wfn files of the PES

    OUTPUT: wfn_energy
        wfn_energy = list of energies for each PES point

    ERROR:
        "Energy values not found in file : file_name" : No energy values in file_name
    ###########################################################################################################
    """
    #INTERNAL VARIABLES:
    wfn_energy = [] #List of wfn files

    #READ FILES 
    for path in A:
        file = open(path, "r")
        lines = file.readlines()
        #ERRORS:
        if "TOTAL ENERGY " not in lines[-1]: #Checks if there is an energy value at the end of wfn file.
            raise ValueError( "Energy values not found in file: ", path)      
        wfn_energy.append(float(lines[-1].split()[3]))
        file.close()
    
    return wfn_energy

def get_aimall_wfx_energies(A):
    #INTERNAL VARIABLES:
    wfx_energy = [] #List of wfn files

    #READ FILES 
    for path in A:
        file = open(path, "r")
        lines = file.readlines()
        #ERRORS:
        if "<Energy = T + Vne + Vee + Vnn>" not in lines[-6]: #Checks if there is an energy value at the end of wfn file.
            raise ValueError( "Energy values not found in file: ", path)      
        wfx_energy.append(float(lines[-5].split()[0]))
        file.close()
    
    return wfx_energy



def intra_property_from_int_file(folders, prop, atom_list):
    """
    ###########################################################################################################
    FUNCTION: intra_property_from_int_file
              get IQA intra-atomic properties from int files

    INPUT: A
        folders = path to _atomicfiles folders
        prop = list of IQA for each atoms e.g.: "['T(A)', 'Vee(A,A)', 'Vne(A,A)']"
        atom_list = list of atomic lables e.g.: [n1, c2, h3, ...]

    OUTPUT: [intra_properties, contributions_list]
        intra_properties = array of array containing the IQA values for eacha atom for each geometry
        contributions_list = list of contributions  organized in the same order as in intra_properties

    ERROR:

    ###########################################################################################################
    """
    #INTERNAL VARIABLES:
    temp1 = [] #Temporary array
    temp2 = [] #Temporary array
    temp3 = [] #Temporary array
    intra_properties = [] #Output
    contributions_list= [] #Output
    
    #READ PROPERTIES FROM .INT FILES
    for folder in folders:
        for atom in atom_list:
            file = open(folder + "/" + atom + ".int", "r")
            lines = file.readlines()
            
            file.close()
            for i in lines:
                if 'IQA Energy Components (see "2EDM Note")' in i:
                    start = lines.index(i)
                elif '2EDM Note:' in i: 
                    end = lines.index(i)
                
            if end >= len(lines) : #Checks the .int file.
                raise ValueError( "File is empty or not exists: " + folder + "/" + atom + ".int") 
                    
            lines = [lines[i] for i in range(start+1, end)]
            for term in prop:
                for i in lines:
                    if (term + '           ') in i:
                        temp1.append(float(i.split()[-1]))
                                  
    #ORGANIZE ARRAY ORDER                
    for j in range(len(prop)):
        for i in range(j, len(temp1), len(prop)):
            temp2.append(temp1[i])
    for j in range(len(atom_list)):
         temp3.append([temp2[i] for i in range(j, len(temp2), len(atom_list))]) 
        
    start = 0 
    for j in range(len(prop)):
        for atom_prop in temp3:
            intra_properties.append([atom_prop[i] for i in range(start, (j+1)*len(folders))])
        start = (j+1)*len(folders)
                         
    #CREATE CONTRIBUTIONS LIST ARRAY:   
    for a in prop:
        for b in atom_list:
            contributions_list.append(a + '-' + b)
       
    return intra_properties, contributions_list



def inter_property_from_int_file(folders, prop, atom_list):
    """
    ###########################################################################################################
    FUNCTION: inter_property_from_int_file
              get IQA interatomic properties from int files

    INPUT: A
        folders = path to _atomicfiles folders
        prop = list of IQA for each atoms e.g.: "['T(A)', 'Vee(A,A)', 'Vne(A,A)']"
        atom_list = list of atomic lables e.g.: [n1, c2, h3, ...]

    OUTPUT: [intra_properties, contributions_list]
        intra_properties = array of array containing the IQA values for eacha atom for each geometry
        contributions_list = list of contributions  organized in the same order as in intra_properties

    ERROR:

    ###########################################################################################################
    """
    #INTERNAL VARIABLES:
    temp1 = [] #Temporary array
    temp2 = [] #Temporary array
    temp3 = [] #Temporary array
    inter_properties = [] #Output
    contributions_list= [] #Output
    
    for path in folders:
        for i in range(len(atom_list)):
            atom1 = atom_list[i]
            for j in range(i+1, len(atom_list)):
                atom2 = atom_list[j]
                file = open(path +"/" + atom1 + "_" + atom2 + ".int", "r")
                lines = file.readlines()
                file.close()
                for i in lines:
                    if ' Energy Components (See "2EDM Note"):' in i:
                        start = lines.index(i)
                    elif '2EDM Note:' in i:
                        end = lines.index(i)
                        
                if end >= len(lines) : #Checks the .int file.
                    raise ValueError( "File is empty or not exists: " + path +"/" + atom1 + "_" + atom2 + ".int") 
                
                lines = [lines[i] for i in range(start+1, end)]
                for term in prop:
                    for i in lines:
                        if (term + '  ') in i:
                                temp1.append(float(i.split()[-1]))
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
       
    return inter_properties, contributions_list


    
    
    