"""
aimall_utils.py v0.1
F. Falcioni, L. J. Duarte, P. L. A. Popelier 

Library with function to submit job to AIMAll and get properties values from output
AIMAll version: 19.10.12
Check for updates at github.com/FabioFalcioni
For details about the method, please see XXXXXXX

"""

import numpy as np 
from typing import List

def distance_A_B(xyz_file : str, atom_A : int, atom_B : int) -> float:
    """distance_A_B gets the distance between atom A and B in an XYZ type file

    :return: distance value 
    :rtype: float
    """    
    # INTERNAL VARIABLES:
    all_coordinates = []

    # WORKING IN THE FILE:
    with open(xyz_file) as f:
        coordinates_list = f.readlines()[3:]  # remove the first 2 lines of xyz file
        for i in range(0, len(coordinates_list)):
            coordinates_of_atom = [float(c) for c in coordinates_list[i].split()[1:]]
            all_coordinates.append(coordinates_of_atom)

    coord_atom_A = all_coordinates[atom_A - 1]
    coord_atom_B = all_coordinates[atom_B - 1]
    x = 0
    y = 1
    z = 2
    # GET DISTANCE
    r_AB = np.sqrt((coord_atom_B[x] - coord_atom_A[x]) ** 2 + (coord_atom_B[y] - coord_atom_A[y]) ** 2 + (
                coord_atom_B[z] - coord_atom_A[z]) ** 2)
    return r_AB

def get_atom_list_wfn(wfn_file : str) -> List[str]:
    """get_atom_list_wfn is a function that gets atom labels from a wavefunction file (.wfn format)

    :raises ValueError: Atomic labels not found
    :return: list of atom labels
    :rtype: List
    """    
    # INTERNAL VARIABLES:
    atom_list = []

    # OPEN FILE:
    file = open(wfn_file, "r")
    lines = file.readlines()  # Convert file into a array of lines
    file.close()  # Close file

    # ERRORS:
    if "(CENTRE " not in lines[2]:
        raise ValueError("Atomic labels not found")  # Checks if atomic list exist inside file

    # GET ATOM LIST:
    for i in range(len(lines)):
        if "(CENTRE" in lines[i]:
            split_line = lines[i].split()
            atom_list.append(split_line[0].lower() + str(split_line[1]))  # uppercase to lowercase

    return atom_list

def get_atom_list_wfx(wfx_file : str) -> List[str]:
    """get_atom_list_wfx is a function that gets atom labels from a wavefunction file (.wfx format)

    :raises ValueError: Atomic labels not found
    :return: list of atom labels
    :rtype: List
    """    
    # INTERNAL VARIABLES:
    atom_list = []

    # OPEN FILE:
    file = open(wfx_file, "r")
    lines = file.readlines()  # Convert file into a array of lines
    file.close()  # Close file

    # ERRORS:
    if "<Nuclear Names>" not in lines[33]:
        raise ValueError("Atomic labels not found")  # Checks if atomic list exist inside file

    # GET ATOM LIST:
    for i in range(len(lines)):
        if "<Nuclear Names>" in lines[i]:
            i += 1
            while "</Nuclear Names>" not in lines[i]:
                split_line = lines[i].split()
                atom_list.append(split_line[0].lower())  # uppercase to lowercase
                i += 1

    return atom_list

def get_aimall_wfn_energies(wfn_files : List[str]) -> List[float]:
    """get_aimall_wfn_energies gets the total energy for each of the wavefunction (.wfn) files in a list

    :raises ValueError: Energy values not found in wavefunction file
    :return: list of energies
    :rtype: List
    """
    # INTERNAL VARIABLES:
    wfn_energy = []  # List of wfn files

    # READ FILES
    for path in wfn_files:
        file = open(path, "r")
        lines = file.readlines()
        # ERRORS:
        if "TOTAL ENERGY " not in lines[-1]:  # Checks if there is an energy value at the end of wfn file.
            raise ValueError("Energy values not found in file: ", path)
        wfn_energy.append(float(lines[-1].split()[3]))
        file.close()

    return wfn_energy

def get_aimall_wfx_energies(wfx_files: List[str]) -> List[float]:
    """get_aimall_wfx_energies gets the total energy for each of the wavefunction (.wfx type) files in a list

    :raises ValueError: Energy values not found in wavefunction file
    :return: list of energies
    :rtype: List
    """
    # INTERNAL VARIABLES:
    wfx_energy = []  # List of wfn files

    # READ FILES
    for path in wfx_files:
        file = open(path, "r")
        lines = file.readlines()
        # ERRORS:
        if "<Energy = T + Vne + Vee + Vnn>" not in lines[
            -6]:  # Checks if there is an energy value at the end of wfn file.
            raise ValueError("Energy values not found in file: ", path)
        wfx_energy.append(float(lines[-5].split()[0]))
        file.close()

    return wfx_energy


def intra_property_from_int_file(folders : List[str], prop : List[str], atom_list : List[str]) -> List:
    """intra_property_from_int_file gets IQA intra-atomic energy values from .int files output from AIMAll

    :raises ValueError: File is empty or does not exist
    :return: Lists of intra-atomic energies with the corresponding intra-atomic labeling
    :rtype: List of floats and strings
    """    
    # INTERNAL VARIABLES:
    temp1 = []  # Temporary array
    temp2 = []  # Temporary array
    temp3 = []  # Temporary array
    intra_properties = []  # Output
    contributions_list = []  # Output

    # READ PROPERTIES FROM .INT FILES
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

            if end >= len(lines):  # Checks the .int file.
                raise ValueError('File is empty or does not exist: ' + folder + "/" + atom + ".int")

            lines = [lines[i] for i in range(start + 1, end)]
            for term in prop:
                for i in lines:
                    if (term + '           ') in i:
                        temp1.append(float(i.split()[-1]))

    # ORGANIZE ARRAY ORDER
    for j in range(len(prop)):
        for i in range(j, len(temp1), len(prop)):
            temp2.append(temp1[i])
    for j in range(len(atom_list)):
        temp3.append([temp2[i] for i in range(j, len(temp2), len(atom_list))])

    start = 0
    for j in range(len(prop)):
        for atom_prop in temp3:
            intra_properties.append([atom_prop[i] for i in range(start, (j + 1) * len(folders))])
        start = (j + 1) * len(folders)

    # CREATE CONTRIBUTIONS LIST ARRAY:
    for a in prop:
        for b in atom_list:
            contributions_list.append(a + '-' + b)

    return intra_properties, contributions_list


def inter_property_from_int_file(folders : List[str], prop : List[str], atom_list : List[str]) -> List[List[float],List[str]]:
    """inter_property_from_int_file gets IQA inter-atomic energy values from .int files output from AIMAll

    :raises ValueError: File is empty or does not exist
    :return: Lists of inter-atomic energies with the corresponding intra-atomic labeling
    :rtype: List of floats and strings
    """    
    # INTERNAL VARIABLES:
    temp1 = []  # Temporary array
    temp2 = []  # Temporary array
    temp3 = []  # Temporary array
    inter_properties = []  # Output
    contributions_list = []  # Output

    for path in folders:
        for i in range(len(atom_list)):
            atom1 = atom_list[i]
            for j in range(i + 1, len(atom_list)):
                atom2 = atom_list[j]
                file = open(path + "/" + atom1 + "_" + atom2 + ".int", "r")
                lines = file.readlines()
                file.close()
                for i in lines:
                    if ' Energy Components (See "2EDM Note"):' in i:
                        start = lines.index(i)
                    elif '2EDM Note:' in i:
                        end = lines.index(i)

                if end >= len(lines):  # Checks the .int file.
                    raise ValueError("File is empty or does not exist: " + path + "/" + atom1 + "_" + atom2 + ".int")

                lines = [lines[i] for i in range(start + 1, end)]
                for term in prop:
                    for i in lines:
                        if (term + '  ') in i:
                            temp1.append(float(i.split()[-1]))
    # ORGANIZE ARRAY ORDER
    for j in range(len(prop)):
        for i in range(j, len(temp1), len(prop)):
            temp2.append(temp1[i])
    for j in range(int(len(atom_list) * (len(atom_list) - 1) / 2)):
        temp3.append([temp2[i] for i in range(j, len(temp2), int(len(atom_list) * (len(atom_list) - 1) / 2))])
    start = 0
    for j in range(len(prop)):
        for atom_prop in temp3:
            inter_properties.append([atom_prop[i] for i in range(start, (j + 1) * len(folders))])
        start = (j + 1) * len(folders)
        # CREATE CONTRIBUTIONS LIST ARRAY:
    for a in prop:
        for i in range(len(atom_list)):
            for j in range(i + 1, len(atom_list)):
                contributions_list.append(a + '-' + atom_list[i] + '_' + atom_list[j])

    return inter_properties, contributions_list


def charge_transfer_and_polarisation_from_int_file(folders : List[str], atom_list : List[str], inter_properties : List[str], xyz_files: List[str]) -> List :
    """charge_transfer_and_polarisation_from_int_file gets IQA polarisation and charge-transfer energies
    values by computing them as Vct = qAqB/rAB and Vpl = Vcl - Vct

    :raises ValueError: File is empty or does not exist
    :return: Returns lists of energies and labels for both Vct and Vpl
    :rtype: List of floats and strings
    """    
    # INTERNAL VARIABLES:
    n = len(atom_list)
    f = len(folders)
    temp1 = []  # Temporary array
    temp2 = []  # Temporary array
    temp3 = []  # Temporary array
    net_charges = []
    charge_transfer_properties = []  # Output
    polarisation_properties = []  # Output
    contributions_list_CT = []  # Output
    contributions_list_PL = []  # Output

    # CREATE CONTRIBTUIONS LIST ARRAY
    for i in range(len(atom_list)):
        for j in range(i + 1, len(atom_list)):
            contributions_list_PL.append('Vpl_IQA(A,B)-' + atom_list[i] + '_' + atom_list[j])
            contributions_list_CT.append('Vct_IQA(A,B)-' + atom_list[i] + '_' + atom_list[j])

    # READ NET-CHARGE PROPERTIES FROM .INT FILES
    for folder in folders:
        net_charge_group = []
        for i in range(0, len(atom_list)):
            file = open(folder + "/" + atom_list[i] + ".int", "r")
            lines = file.readlines()
            file.close()

            for line in lines:
                if 'Results of the basin integration:' in line:
                    start = lines.index(line)
                elif '|Dipole|' in line:
                    end = lines.index(line)
            if end >= len(lines):  # Checks the .int file.
                raise ValueError("File is empty or does not exist: " + folder + "/" + atom_list[i] + ".int")

            lines = [lines[j] for j in range(start + 1, end)]
            Q = (float(lines[0].split()[-4]))
            net_charge_group.append(Q)
        net_charges.append(net_charge_group)

    # GET CHARGE TRANSFER TERMS
    for k in range(len(net_charges)):
        for i in range(len(atom_list)):
            for j in range(i + 1, len(atom_list)):
                temp1.append((net_charges[k][i] * net_charges[k][j]) / ((distance_A_B(xyz_files[k], i + 1, j + 1))*1.8897259886))

    # ORGANIZE CHARGE TRANSFER ARRAY ORDER
    for j in range(int(n * (n - 1) / 2)):
        temp2.append([temp1[i] for i in range(j, len(temp1), int(n * (n - 1) / 2))])
    start = 0
    for atom_prop in temp2:
        charge_transfer_properties.append([atom_prop[i] for i in range(start, f)])

    # Isolate Vcl terms
    classical_properties = inter_properties[:len(charge_transfer_properties)]

    # OBTAIN POLARISATION TERMS AS Vpl = Vcl - Vct
    for i in range(len(classical_properties)):
        for j in range(len(classical_properties[i])):
            temp3.append(classical_properties[i][j] - charge_transfer_properties[i][j])

    # ORGANIZE POLARISATION ARRAY ORDER
    polarisation_properties = [temp3[i * f:(i + 1) * f] for i in range((len(temp3) + f - 1) // f)]

    return charge_transfer_properties, contributions_list_CT, polarisation_properties, contributions_list_PL

