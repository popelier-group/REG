"""
dft-d3_utils.py v0.0
F. Falcioni, P. L. A. Popelier

Library with function to submit job to DFT-D3 and get properties values from output
Check for updates at github.com/FabioFalcioni

Please, report bugs and issues to fabio.falcioni@manchester.ac.uk
coded by F.Falcioni
"""

import os


def disp_property_from_dftd3_file(folders, atom_list):
    """
    ###########################################################################################################
    Added by Aël Cador - 17 March 2020
    FUNCTION: disp_property_from_dftd3_file
              get dispersion interatomic properties from DFT-D3 result files (Grimme DFT-D3)

    INPUT:
        folders = path to DFT-D3 result files
        atom_list = list of atomic labels* e.g.: [n1, c2, h3, ...]
        disp_lim = list of atomic label numbers which separate the fragments (the first atom of each fragment starting with the 2nd fragment)

    OUTPUT: [intra_properties, contributions_list]
        disp_properties = array of array containing the dispersion energy values for each atom for each geometry
        contributions_list = list of contributions  organized in the same order as in disp_properties

    ERROR:
        'File is empty or does not exist:'

    ###########################################################################################################
    """
    # INTERNAL VARIABLES:
    temp1 = []  # Temporary array
    temp2 = []  # Temporary array
    #   temp3 = [] #Temporary array
    disp_properties = []  # Output
    contributions_list = []  # Output
    n = len(atom_list)

    for i in range(len(atom_list)):
        for j in range(i + 1, len(atom_list)):
            contributions_list.append(
                "E_Disp(A,B)-" + atom_list[i] + "_" + atom_list[j]
            )

    for path in folders:
        #        for i in range(len(atom_list)):
        #            atom1 = atom_list[i]
        #            for j in range(i+1, len(atom_list)):
        #                atom2 = atom_list[j]
        file = open(path, "r")  # +"/" + atom1 + "_" + atom2 + ".int", "r")
        lines = file.readlines()
        file.close()
        for i in lines:
            if "analysis of pair-wise terms (in kcal/mol)" in i:
                start = lines.index(i) + 2
            elif "distance range (Angstroem) analysis" in i:
                end = lines.index(i) - 1  # check if both are valid*

        if end >= len(lines):  # Checks the DFT-D3 file.
            raise ValueError(
                "File is empty or does not exist: " + path
            )  # +"/" + atom1 + "_" + atom2 + ".int")

        lines = [lines[i] for i in range(start, end)]
        #               for term in prop:
        for i in lines:
            #                   if (term + '  ') in i:
            temp1.append(float(i.split()[-1]))
    # ORGANIZE ARRAY ORDER
    #    for j in range(len(prop)):
    #        for i in range(j, len(temp1), len(prop)):
    #            temp2.append(temp1[i])
    for j in range(int(n * (n - 1) / 2)):
        temp2.append(
            [temp1[i] / 627.503 for i in range(j, len(temp1), int(n * (n - 1) / 2))]
        )  # temp2 ?
    start = 0
    #   for j in range(len(prop)):
    for atom_prop in temp2:
        disp_properties.append(
            [atom_prop[i] for i in range(start, len(folders))]
        )  # (j+1)*
    # start = (j+1)*len(folders)
    # CREATE CONTRIBUTIONS LIST ARRAY:

    return disp_properties, contributions_list


def run_DFT_D3(program_path, xyz_file, functional, BJ_Damping=True):
    """
    ###########################################################################################################
    FUNCTION: run_DFT_D3
              run Grimme DFT-D3 program for selected xyz file given the path of the program

    INPUT:
        program_path = path to DFT-D3 program
        xyz_file = xyz format file
        functional = functional used for
        BJ_Dumping = Becke Johnson dumping (default=True)

    ERRORS:
        'Insert a functional that works with BJ dumping and AIMAll'
        'Insert a functional that works with AIMAll'

    Note (06/12/2020): AIMAll works with LSDA, B3LYP, PBE, PBE0, M062X
    ###########################################################################################################
    """
    if BJ_Damping == True:
        functional_list = ["B3-LYP", "PBE0", "PBE"]
        if functional.upper() in functional_list:
            os.system(
                program_path
                + " "
                + xyz_file
                + " -func "
                + functional.lower()
                + " -bj -anal > dft-d3.log"
            )
        else:
            raise ValueError(
                "Insert a functional that works with BJ dumping and AIMAll"
            )

    else:
        functional_list = ["B3-LYP", "M062X", "PBE", "PBE0"]
        if functional.upper() in functional_list:
            os.system(
                program_path
                + " "
                + xyz_file
                + " -func "
                + functional.lower()
                + " -zero  -anal > dft-d3.log"
            )
        else:
            raise ValueError("Insert a functional that works with AIMAll")
