"""
auto_reg.py v0.1
F. Falcioni, P. L. A. Popelier

Library with function to run a REG analysis
Check for updates at github.com/FabioFalcioni
For details about the method, please see XXXXXXX

Please, report bugs and issues to fabio.falcioni@manchester.ac.uk
coded by F.Falcioni

NOTE: The automatic analysis works if this file is run with python3 inside a folder containing all the REG points (
saved in numbered folders) """

# IMPORT LIBRARIES
import sys
sys.path.insert(1, '/mnt/iusers01/pp01/v69787ff/REG/')  # PLEASE INSERT THE PATH OF REG.py folder installation
import reg
import aimall_utils as aim_u
import numpy as np
import pandas as pd
import reg_vis as rv
import gaussian_utils as gauss_u
import dftd3_utils as disp_u
import re
import os
import time

### STARTING TIMER ###
start_time = time.time()
##############################    VARIABLES    ##################################

SYS = 'SN2_ClBr_B3LYP'  # name of the system

### PES Critical points options ###
POINTS = 4  # number of points for "find_critical" function
AUTO = True  # Search for critical points
turning_points = []  # manually put critical points in the PES if necessary
# NOTE: If analysing only a single segment (i.e. the PES has no critical points) please put AUTO=False and tp=[]

# DEFINE THE DESIRED TERMS:
intra_prop = ['E_IQA_Intra(A)']  # chose the AIMAll intra atomic properties to analyse
intra_prop_names = ['Eintra']  # names of the properties shown in the output
inter_prop = ['VC_IQA(A,B)', 'VX_IQA(A,B)', 'E_IQA_Inter(A,B)']  # chose the AIMAll inter atomic properties to analyse
inter_prop_names = ['Vcl', 'Vxc', 'Einter']  # names of the properties shown in the output

REVERSE = True  # Reverse the REG points

INFLEX = False

### CONTROL COORDINATE OPTIONS ###
CONTROL_COORDINATE_TYPE = ''  # 'Scan' or 'IRC'. If empty ('') then default will be used
Scan_Atoms = [1, 6]  # list of the atoms used for PES Scan (i.e. ModRedundant option in Gaussian)
IRC_output = ''  # insert the g16 output file path if using IRC as control coordinate

CHARGE_TRANSFER_POLARISATION = True  # Split the classical electrostatic term into polarisation and monopolar charge-transfer

DISPERSION = True  # Run DFT-D3 program to consider dispersion
# NOTE: This only works if you have DFT-D3 program installed in the same machine https://www.chemie.uni-bonn.de/pctc/mulliken-center/software/dft-d3/
### DISPERSION OPTIONS ###
DFT_D3_PATH = '/mnt/iusers01/pp01/v69787ff/DFT-D3/dftd3'  # insert path of DFT-D3 program
DISP_FUNCTIONAL = 'B3-LYP'  # insert the functional used for D3 correction
BJ_DAMPING = True  # Becke-Johnson Damping

WRITE = True  # write csv files for energy values and REG analysis
SAVE_FIG = True  # save figures
ANNOTATE = True  # annotate figures
DETAILED_ANALYSIS = True
LABELS = True  # label figures
n_terms = 4  # number of terms to rank in figures and tables

##################################################################################

###############################################################################
#                                                                             #
#                           AUTOMATIC FILES SETUP                             #
#                                                                             #
###############################################################################

# DEFINE PATHS AND FILES AUTOMATICALLY:
cwd = str(os.getcwd())

wfn_filelist = []

# Get wfn files
for k in os.walk(r"%s" % cwd):
    for i in range(0, len(k), 1):
        for j in range(0, len(k[i]), 1):
            if ".wfn" in k[i][j]:
                wfn_filelist.append(k[i][j])

# Finds the root, folders and all the files within the CWD and stores them as variables:
# 'root', 'dirs' and 'allfiles'
folderlist = []
i = 0
for root, dirs, allfiles in os.walk(r'%s' % cwd):
    for folder in dirs:
        folderlist.append(folder)
        i = i + 1

# Sort to ensure order
folderlist = sorted(folderlist)
wfn_filelist = sorted(wfn_filelist)

# Find all the g16 single point energy files
# NOTE: this works if the output files end with ".out"
g16_out_file = []
i = 0
for root, dirs, allfiles in os.walk(r'%s' % cwd):
    for file in allfiles:
        if file.endswith('.out'):
            g16_out_file.append(file)

# Find all the REG points folders and sort them by number
reg_folders = []
j = 0
for i in range(0, len(folderlist)):
    try:
        if os.path.isfile((cwd + '/' + folderlist[i] + '/' + wfn_filelist[j])):
            reg_folders.append(folderlist[i])
            j = j + 1
    except:
        pass
reg_folders.sort(key=lambda f: int(re.sub('\D', '', f)))

if REVERSE:
    reg_folders = reg_folders[::-1]

os.chdir(cwd)  # working directory
# Create results directory
access_rights = 0o755
try:
    os.mkdir(SYS + "_results", access_rights)
except OSError:
    print("Creation of the directory {a}/{b}_results failed or has already been created".format(a=cwd,b=SYS))
else:
    print("Successfully created the directory {a}/{b}_results".format(a=cwd,b=SYS))

# GET ATOM LIST FROM ANY .WFN FILE:
atoms = aim_u.get_atom_list(cwd + '/' + reg_folders[0] + '/' + wfn_filelist[0])

# Arrange files and folders in lists
wfn_files = [cwd + '/' + reg_folders[i] + '/' + wfn_filelist[i] for i in range(0, len(reg_folders))]
atomic_files = [cwd + '/' + reg_folders[i] + '/' + wfn_filelist[i][:-4] + '_atomicfiles' for i in
                range(0, len(reg_folders))]
g16_files = [cwd + '/' + reg_folders[i] + '/' + g16_out_file[i] for i in range(0, len(reg_folders))]
xyz_files = [gauss_u.get_xyz_file(file) for file in g16_files]

# Get control coordinate list
if CONTROL_COORDINATE_TYPE == 'Scan':
    cc = gauss_u.get_control_coordinates_PES_Scan(g16_files, Scan_Atoms)
    X_LABEL = r"Control Coordinate [$\AA$]"
elif CONTROL_COORDINATE_TYPE == 'IRC':
    cc = gauss_u.get_control_coordinates_IRC_g16(IRC_output)
    X_LABEL = r"Control Coordinate r[$\AA$]"
else:
    cc = [int(reg_folders[i]) for i in range(0, len(reg_folders))]
    X_LABEL = "Control Coordinate [REG step]"
cc = np.array(cc)
if REVERSE:
    cc = -cc

### INTRA AND INTER ENERGY TERMS ###


# GET TOTAL ENERGY FROM THE .WFN FILES:
total_energy_wfn = aim_u.get_aimall_wfn_energies(wfn_files)
total_energy_wfn = np.array(total_energy_wfn)
# GET INTRA-ATOMIC TERMS:
iqa_intra, iqa_intra_header = aim_u.intra_property_from_int_file(atomic_files, intra_prop, atoms)
iqa_intra_header = np.array(iqa_intra_header)  # used for reference
iqa_intra = np.array(iqa_intra)
# GET INTER-ATOMIC TERMS:
iqa_inter, iqa_inter_header = aim_u.inter_property_from_int_file(atomic_files, inter_prop, atoms)
iqa_inter_header = np.array(iqa_inter_header)  # used for reference
iqa_inter = np.array(iqa_inter)

###############################################################################
#                                                                             #
#                               REG ANALYSIS                                  #
#                                                                             #
###############################################################################

# INTRA ATOMIC CONTRIBUTION
reg_intra = reg.reg(total_energy_wfn, cc, iqa_intra, np=POINTS, critical=AUTO, inflex=INFLEX,
                    critical_index=turning_points)
# INTER ATOMIC CONTRIBUTION
reg_inter = reg.reg(total_energy_wfn, cc, iqa_inter, np=POINTS, critical=AUTO, inflex=INFLEX,
                    critical_index=turning_points)

# GET CT and PL TERMS:
if CHARGE_TRANSFER_POLARISATION:
    iqa_charge_transfer_terms, iqa_charge_transfer_headers, iqa_polarisation_terms, iqa_polarisation_headers = aim_u.charge_transfer_and_polarisation_from_int_file(
        atomic_files, atoms, iqa_inter, xyz_files)
    iqa_charge_transfer_headers = np.array(iqa_charge_transfer_headers)
    iqa_polarisation_headers = np.array(iqa_polarisation_headers)
    iqa_polarisation_terms = np.array(iqa_polarisation_terms)
    iqa_charge_transfer_terms = np.array(iqa_charge_transfer_terms)
    # CHARGE TRANSFER CONTRIBUTION
    reg_ct = reg.reg(total_energy_wfn, cc, iqa_charge_transfer_terms, np=POINTS, critical=AUTO, inflex=INFLEX,
                     critical_index=turning_points)
    # POLARISATION CONTRIBUTION
    reg_pl = reg.reg(total_energy_wfn, cc, iqa_polarisation_terms, np=POINTS, critical=AUTO, inflex=INFLEX,
                     critical_index=turning_points)

### DISPERSION ANALYSIS ###
if DISPERSION:
    for i in range(0, len(reg_folders)):
        if os.path.isfile((cwd + '/' + reg_folders[i] + '/' + wfn_filelist[i])):
            os.chdir((cwd + '/' + reg_folders[i] + '/'))
            xyz_file = gauss_u.get_xyz_file(g16_out_file[i])
            disp_u.run_DFT_D3(DFT_D3_PATH, xyz_file, DISP_FUNCTIONAL)

    folders_disp = [cwd + '/' + reg_folders[i] + '/dft-d3.log' for i in range(0, len(reg_folders))]
    # GET INTER-ATOMIC DISPERSION TERMS:
    iqa_disp, iqa_disp_header = disp_u.disp_property_from_dftd3_file(folders_disp, atoms)
    iqa_disp_header = np.array(iqa_disp_header)  # used for reference
    iqa_disp = np.array(iqa_disp)
    # REG
    reg_disp = reg.reg(total_energy_wfn, cc, iqa_disp, np=POINTS, critical=AUTO, inflex=INFLEX,
                       critical_index=turning_points)
    total_energy_dispersion = sum(iqa_disp)

# CALCULATE TOTAL ENERGIES
total_energy_iqa = sum(iqa_inter[:(len(atoms)*(len(atoms)-1))]) + sum(iqa_intra[:len(atoms)])  # used to calculate the integration error

# CALCULATE RECOVERY ERROR
if DISPERSION:
    rmse_integration = reg.integration_error(total_energy_wfn, total_energy_iqa + total_energy_dispersion)
else:
    rmse_integration = reg.integration_error(total_energy_wfn, total_energy_iqa)
print('Integration error [kJ/mol](RMSE)')
print(rmse_integration[1])

###############################################################################
#                                                                             #
#                             WRITE CSV FILES                                 #
#                                                                             #
###############################################################################
os.chdir(cwd + '/' + SYS + "_results")
dataframe_list = []

if WRITE:
    # initialise excel file
    writer = pd.ExcelWriter(path=cwd + '/' + SYS + "_results/REG.xlsx", engine='xlsxwriter')
    # ENERGY and CONTROL  COORDINATE ONLY .CSV FILES
    df_energy_output = pd.DataFrame()
    df_energy_output['WFN'] = total_energy_wfn
    df_energy_output['IQA'] = total_energy_iqa
    df_energy_output.index = cc
    if DISPERSION:
        df_energy_output['D3'] = total_energy_dispersion
    df_energy_output.to_csv('total_energy.csv', sep=',')
    df_energy_output.to_excel(writer, sheet_name="total_energy")

    # INTER AND INTRA PROPERTIES RE-ARRANGEMENT
    list_property_final = []
    final_properties_comparison = []
    for i in range(len(reg_inter[0])):
        list_property_sorted = []
        properties_comparison = []
        df_inter = rv.create_term_dataframe(reg_inter, iqa_inter_header,i)
        df_intra = rv.create_term_dataframe(reg_intra, iqa_intra_header, i)
        for j in range(len(inter_prop)):
            df_property = rv.filter_term_dataframe(df_inter, inter_prop[j], inter_prop_names[j])
            if j <= 1:
                properties_comparison.append(df_property)
            df_property.to_csv(inter_prop_names[j] + "_seg_" + str(i + 1) + ".csv", sep=',')
            df_property.to_excel(writer, sheet_name=inter_prop_names[j] + "_seg_" + str(i + 1))
            list_property_sorted.append(
                pd.concat([df_property[-n_terms:], df_property[:n_terms]], axis=0).sort_values('REG'))
        for j in range(len(intra_prop)):
            df_property = rv.filter_term_dataframe(df_intra, intra_prop[j], intra_prop_names[j])
            if j == 0:
                properties_comparison.append(df_property)
            df_property.to_csv(intra_prop_names[j] + "_seg_" + str(i + 1) + ".csv", sep=',')
            df_property.to_excel(writer, sheet_name=intra_prop_names[j] + "_seg_" + str(i + 1))
            list_property_sorted.append(pd.concat([df_property], axis=0).sort_values('REG'))
        list_property_final.append(list_property_sorted)
        final_properties_comparison.append(properties_comparison)

    # DISPERSION OUTPUT
    if DISPERSION:
        df_dispersion_sorted = pd.DataFrame()
        disp_name_old = 'E_Disp(A,B)'
        disp_name_new = 'Vdisp'
        for i in range(len(reg_inter[0])):
            df_disp = rv.create_term_dataframe(reg_disp, iqa_disp_header,i)
            df_disp_new = rv.filter_term_dataframe(df_disp, disp_name_old, disp_name_new)
            df_disp_new.to_csv(disp_name_new + "_seg_" + str(i + 1) + ".csv", sep=',')
            df_disp_new.to_excel(writer, sheet_name=disp_name_new + "_seg_" + str(i + 1))
            df_disp_new.dropna(axis=0, how='any', subset=None,
                           inplace=True)  # get rid of "NaN" terms which have a null REG Value
            df_dispersion_sorted = pd.concat([df_dispersion_sorted.reset_index(drop=True),
                                              pd.concat([df_disp_new[-n_terms:], df_disp_new[:n_terms]],
                                                        axis=0).sort_values(
                                                  'REG').reset_index(drop=True)], axis=1)
        df_dispersion_sorted.to_csv('REG_' + disp_name_new + '_analysis.csv', sep=',')
        df_dispersion_sorted.to_excel(writer, sheet_name="REG_" + disp_name_new)
        rv.pandas_REG_dataframe_to_table(df_dispersion_sorted, 'REG_' + disp_name_new + '_table', SAVE_FIG=SAVE_FIG)

    # CHARGE-TRANSFER and POLARISATION
    if CHARGE_TRANSFER_POLARISATION:
        df_ct_pl_sorted = pd.DataFrame()
        for i in range(len(reg_inter[0])):
            df_pl = rv.filter_term_dataframe(rv.create_term_dataframe(reg_pl, iqa_polarisation_headers,i),
                                             'Vpl_IQA(A,B)',
                                             'Vpl')
            df_ct = rv.filter_term_dataframe(rv.create_term_dataframe(reg_ct, iqa_charge_transfer_headers,i),
                                             'Vct_IQA(A,B)',
                                             'Vct')
            df_pl.to_csv("Vpl_seg_" + str(i + 1) + ".csv", sep=',')
            df_pl.to_excel(writer, sheet_name="Vpl_seg_" + str(i + 1))
            df_ct.to_csv("Vct_seg_" + str(i + 1) + ".csv", sep=',')
            df_ct.to_excel(writer, sheet_name="Vct_seg_" + str(i + 1))
            df_temp = pd.concat([df_pl, df_ct]).sort_values('REG').reset_index(drop=True)
            df_ct_pl_sorted = pd.concat([df_ct_pl_sorted.reset_index(drop=True),
                                         pd.concat([df_temp[-n_terms:], df_temp[:n_terms]], axis=0).sort_values(
                                             'REG').reset_index(drop=True)], axis=1)
        df_ct_pl_sorted.to_csv('REG_Vct-Vpl_analysis.csv', sep=',')
        df_ct_pl_sorted.to_excel(writer, sheet_name='REG_Vct-Vpl')
        rv.pandas_REG_dataframe_to_table(df_ct_pl_sorted, 'REG_Vct-Vpl_table', SAVE_FIG=SAVE_FIG)

    # OUTPUT OF ALL INTER AND INTRA TERMS SELECTED BY THE USER
    all_prop_names = inter_prop_names + intra_prop_names
    for i in range(len(inter_prop) + len(intra_prop)):
        df_property_sorted = pd.DataFrame()
        for j in range(len(reg_inter[0])):
            df_property_sorted = pd.concat([df_property_sorted, list_property_final[j][i]], axis=1)
        df_property_sorted.to_csv('REG_' + all_prop_names[i] + '_analysis.csv', sep=',')
        df_property_sorted.to_excel(writer, sheet_name='REG_' + all_prop_names[i])
        rv.pandas_REG_dataframe_to_table(df_property_sorted, 'REG_' + all_prop_names[i] + '_table', SAVE_FIG=SAVE_FIG)

    # FINAL COMPARISON
    df_final_sorted = pd.DataFrame()
    for i in range(len(reg_inter[0])):
        df_final = pd.DataFrame()
        for j in range(3):
            df_final = pd.concat([df_final, final_properties_comparison[i][j]])
        if DISPERSION:
            df_final = pd.concat([df_final, df_disp_new])
        dataframe_list.append(df_final)
        df_final = df_final.sort_values('REG').reset_index(drop=True)
        df_final_sorted = pd.concat([df_final_sorted.reset_index(drop=True),
                                     pd.concat([df_final[-n_terms:], df_final[:n_terms]], axis=0).sort_values(
                                         'REG').reset_index(drop=True)], axis=1)
        df_final.to_csv('REG_full_comparison_seg_' + str(i + 1) + '.csv', sep=',')
        df_final.to_excel(writer, sheet_name='REG_full_comparison_seg_' + str(i+1))
    df_final_sorted.to_csv('REG_final_analysis.csv', sep=',')
    df_final_sorted.to_excel(writer, sheet_name='REG_final')
    rv.pandas_REG_dataframe_to_table(df_final_sorted, 'REG_final_table', SAVE_FIG=SAVE_FIG)

    writer.save()

    rv.plot_violin([dataframe_list[i]['R'] for i in range(len(reg_inter[0]))], save=SAVE_FIG,
                   file_name='violin.png')  # Violing plot of R vs Segments

###############################################################################
#                                                                             #
#                                   GRAPHS                                    #
#                                                                             #
###############################################################################
if AUTO:
    critical_points = reg.find_critical(total_energy_wfn, cc, min_points=POINTS, use_inflex=INFLEX)
else:
    critical_points = turning_points

rv.plot_segment(cc, 2625.50 * (total_energy_wfn - (sum(total_energy_wfn) / len(total_energy_wfn))), critical_points,
                annotate=ANNOTATE,
                label=LABELS,
                y_label=r'Relative Energy [$kJ.mol^{-1}$]', x_label=X_LABEL, title=SYS,
                save=SAVE_FIG, file_name='REG_analysis.png')

if DETAILED_ANALYSIS:
    for i in range(len(reg_inter[0])):
        rv.generate_data_vis(dataframe_list[i], [dataframe_list[i]['R'] for i in range(len(reg_inter[0]))],
                             n_terms, save=SAVE_FIG, file_name='detailed_seg_' + str(i + 1) + '.png',
                             title=SYS + ' seg. ' + str(i + 1))

###ENDING TIMER ###
print("--- Total time for REG Analysis: {s} minutes ---".format(s=((time.time() - start_time) / 60)))
