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

#DEFINE THE DESIRED TERMS:
intra_prop = ['E_IQA_Intra(A)'] #intra atomic properties
inter_prop = ['VC_IQA(A,B)', 'VX_IQA(A,B)'] #inter atomic properties
REVERSE = True #Reverse the REG points

INFLEX = False

### CONTROL COORDINATE OPTIONS ###
CONTROL_COORDINATE_TYPE = 'Scan'  # 'Scan' or 'IRC'. If empty ('') then default will be used
Scan_Atoms = [1, 6]  # list of the atoms used for PES Scan (i.e. ModRedundant option in Gaussian)
IRC_output = ''  # insert the g16 output file path if using IRC as control coordinate

CHARGE_TRANSFER_POLARISATION = True  # Split the classical electrostatic term into polarisation and monopolar charge-transfer
FULL_INTERATOMIC_CONTRIBUTION = True  # Combine classical and exchange-correlation contributions

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
#Create results directory
access_rights = 0o755
try:
    os.mkdir(SYS + "_results", access_rights)
except OSError:
    print("Creation of the directory %s failed" % cwd)
else:
    print("Successfully created the directory %s" % cwd)

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
    X_LABEL = "Control Coordinate " + r"[$\AA$]"
elif CONTROL_COORDINATE_TYPE == 'IRC':
    cc = gauss_u.get_control_coordinates_IRC_g16(IRC_output)
    X_LABEL = X_LABEL = "Control Coordinate " + r"[$\AA$]"
else:
    cc = [int(reg_folders[i]) for i in range(0, len(reg_folders))]
    X_LABEL = "Control Coordinate [REG step]"

### INTRA AND INTER ENERGY TERMS ###

# DEFINE THE DESIRED TERMS:
intra_prop = ['E_IQA_Intra(A)']  # intra atomic properties
inter_prop = ['VC_IQA(A,B)', 'VX_IQA(A,B)']  # inter atomic properties

# GET TOTAL ENERGY FROM THE .WFN FILES:
total_energy_wfn = aim_u.get_aimall_wfn_energies(wfn_files)
total_energy_wfn = np.array(total_energy_wfn)

# GET INTRA-ATOMIC TERMS:
iqa_intra = aim_u.intra_property_from_int_file(atomic_files, intra_prop, atoms)[0]
iqa_intra_header = np.array(
    aim_u.intra_property_from_int_file(atomic_files, intra_prop, atoms)[1])  # used for reference
iqa_intra = np.array(iqa_intra)
# GET INTER-ATOMIC TERMS:
iqa_inter = aim_u.inter_property_from_int_file(atomic_files, inter_prop, atoms)[0]
iqa_inter_header = np.array(
    aim_u.inter_property_from_int_file(atomic_files, inter_prop, atoms)[1])  # used for reference
iqa_inter = np.array(iqa_inter)

###############################################################################
#                                                                             #
#                               REG ANALYSIS                                  #
#                                                                             #
###############################################################################

# INTRA ATOMIC CONTRIBUTION
reg_intra = reg.reg(total_energy_wfn, cc, iqa_intra, np=POINTS, critical=AUTO, inflex=INFLEX, critical_index=turning_points)
# INTER ATOMIC CONTRIBUTION
reg_inter = reg.reg(total_energy_wfn, cc, iqa_inter, np=POINTS, critical=AUTO, inflex=INFLEX, critical_index=turning_points)

# GET CT and PL TERMS:
if CHARGE_TRANSFER_POLARISATION:
    iqa_charge_transfer_terms, iqa_charge_transfer_headers, iqa_polarisation_terms, iqa_polarisation_headers = aim_u.charge_transfer_and_polarisation_from_int_file(
        atomic_files, atoms, iqa_inter, xyz_files)
    iqa_charge_transfer_headers = np.array(iqa_charge_transfer_headers)
    iqa_polarisation_headers = np.array(iqa_polarisation_headers)
    iqa_polarisation_terms = np.array(iqa_polarisation_terms)
    iqa_charge_transfer_terms = np.array(iqa_charge_transfer_terms)
    # CHARGE TRANSFER CONTRIBUTION
    reg_ct = reg.reg(total_energy_wfn, cc, iqa_charge_transfer_terms, np=POINTS, critical=AUTO, inflex=INFLEX, critical_index=turning_points)
    # POLARISATION CONTRIBUTION
    reg_pl = reg.reg(total_energy_wfn, cc, iqa_polarisation_terms, np=POINTS, critical=AUTO, inflex=INFLEX, critical_index=turning_points)

# GET COMBINED INTER-ATOMIC TERMS:
if FULL_INTERATOMIC_CONTRIBUTION:
    E_iqa_inter_terms = aim_u.E_IQA_from_int_file(atomic_files, atoms)[0]
    E_iqa_inter_headers = np.array(aim_u.E_IQA_from_int_file(atomic_files, atoms)[1])
    E_iqa_inter_terms = np.array(E_iqa_inter_terms)
    # INTER ONLY ATOMIC CONTRIBUTION
    reg_inter_full = reg.reg(total_energy_wfn, cc, E_iqa_inter_terms, np=POINTS, critical=AUTO, inflex=INFLEX, critical_index=turning_points)
### DISPERSION ANALYSIS ###
if DISPERSION:
    for i in range(0, len(reg_folders)):
        if os.path.isfile((cwd + '/' + reg_folders[i] + '/' + wfn_filelist[i])):
            os.chdir((cwd + '/' + reg_folders[i] + '/'))
            xyz_file = gauss_u.get_xyz_file(g16_out_file[i])
            disp_u.run_DFT_D3(DFT_D3_PATH, xyz_file, DISP_FUNCTIONAL)

    folders_disp = [cwd + '/' + reg_folders[i] + '/dft-d3.log' for i in range(0, len(reg_folders))]
    # GET INTER-ATOMIC DISPERSION TERMS:
    iqa_disp = disp_u.disp_property_from_dftd3_file(folders_disp, atoms)[0]
    iqa_disp_header = np.array(disp_u.disp_property_from_dftd3_file(folders_disp, atoms)[1])  # used for reference
    iqa_disp = np.array(iqa_disp)
    # REG
    reg_disp = reg.reg(total_energy_wfn, cc, iqa_disp, np=POINTS, critical=AUTO, inflex=INFLEX, critical_index=turning_points)
    total_energy_dispersion = sum(iqa_disp)

# CALCULATE TOTAL ENERGIES
total_energy_iqa = sum(iqa_inter) + sum(iqa_intra)  # used to calculate the integration error

# CALCULATE RECOVERY ERROR
rmse_integration = reg.integration_error(total_energy_wfn, total_energy_iqa)
print('Integration error [kJ/mol](RMSE)')
print(rmse_integration[1])

###############################################################################
#                                                                             #
#                             WRITE CSV FILES                                 #
#                                                                             #
###############################################################################
os.chdir(SYS + "_results")
dataframe_list = []

if WRITE == True:
    # Initialise necessary dataframes
    df_final = pd.DataFrame()
    df_xc_sorted = pd.DataFrame()
    df_cl_sorted = pd.DataFrame()
    df_dispersion_sorted = pd.DataFrame()
    df_intra_sorted = pd.DataFrame()
    df_ct_pl_sorted = pd.DataFrame()
    df_inter_full_sorted = pd.DataFrame()
    for i in range(len(reg_inter[0])):
        # INTER-ATOMIC ENERGY DATAFRAME
        temp = [reg_inter[0][i], reg_inter[1][i]]
        df_inter = pd.DataFrame(temp).transpose()
        df_inter.columns = ["REG", "R"]
        df_inter.index = iqa_inter_header
        df_inter.sort_values('REG').to_csv("Inter_seg_" + str(i + 1) + ".csv", sep=',')
        df_inter = df_inter.rename_axis('TERM').reset_index()

        # INTRA-ATOMIC ENERGY DATAFRAME
        temp = [reg_intra[0][i], reg_intra[1][i]]
        df_intra = pd.DataFrame(temp).transpose()
        df_intra.columns = ["REG", "R"]
        df_intra.index = iqa_intra_header
        df_intra.sort_values('REG').to_csv("E_Intra_seg_" + str(i + 1) + ".csv", sep=',')
        df_intra = df_intra.rename_axis('TERM').reset_index()

        # DISPERTION ENERGY DATAFRAME
        if DISPERSION:
            temp = [reg_disp[0][i], reg_disp[1][i]]
            df_disp = pd.DataFrame(temp).transpose()
            df_disp.columns = ["REG", "R"]
            df_disp.index = iqa_disp_header
            df_disp.sort_values('REG').to_csv("Disp_seg_" + str(i + 1) + ".csv", sep=',')
            df_disp = df_disp.rename_axis('TERM').reset_index()
            df_disp.dropna(axis=0, how='any', subset=None,
                           inplace=True)  # get rid of "NaN" terms which have a null REG Value

            df_disp = df_disp.sort_values('REG').reset_index()
            df_disp['TERM'] = df_disp['TERM'].str.replace("E_Disp\(A,B\)-", 'Vdisp(')
            df_disp['TERM'] = df_disp['TERM'].str.replace("_", ',')
            df_disp['TERM'] = df_disp['TERM'] + ')'
            df_dispersion_sorted = pd.concat(
                [df_dispersion_sorted, pd.concat([df_disp[-n_terms:], df_disp[:n_terms]], axis=0).sort_values('REG')],
                axis=1)

        if CHARGE_TRANSFER_POLARISATION:
            # POLARISATION ENERGY DATAFRAME
            temp = [reg_pl[0][i], reg_pl[1][i]]
            df_pl = pd.DataFrame(temp).transpose()
            df_pl.columns = ["REG", "R"]
            df_pl.index = iqa_polarisation_headers
            df_pl.sort_values('REG').to_csv("Vpl_seg_" + str(i + 1) + ".csv", sep=',')
            df_pl = df_pl.rename_axis('TERM').reset_index()

            # CHARGE-TRANSFER ENERGY DATAFRAME
            temp = [reg_ct[0][i], reg_ct[1][i]]
            df_ct = pd.DataFrame(temp).transpose()
            df_ct.columns = ["REG", "R"]
            df_ct.index = iqa_charge_transfer_headers
            df_ct.sort_values('REG').to_csv("Vct_seg_" + str(i + 1) + ".csv", sep=',')
            df_ct = df_ct.rename_axis('TERM').reset_index()

            df_temp_ct_pl = pd.concat([df_pl, df_ct]).sort_values('REG').reset_index()

            df_ct_pl_sorted = pd.concat([df_ct_pl_sorted, pd.concat([df_temp_ct_pl[-n_terms:], df_temp_ct_pl[:n_terms]],
                                                                    axis=0).sort_values('REG')], axis=1)

        # FULL INTER-ATOMIC ENERGY DATAFRAME
        temp = [reg_inter_full[0][i], reg_inter_full[1][i]]
        df_inter_full = pd.DataFrame(temp).transpose()
        df_inter_full.columns = ["REG", "R"]
        df_inter_full.index = E_iqa_inter_headers
        df_inter_full.sort_values('REG').to_csv("E_full_inter_seg_" + str(i + 1) + ".csv", sep=',')
        df_inter_full = df_inter_full.rename_axis('TERM').reset_index()

        df_inter_full = df_inter_full.sort_values('REG').reset_index()
        df_inter_full_sorted = pd.concat([df_inter_full_sorted,
                                          pd.concat([df_inter_full[-n_terms:], df_inter_full[:n_terms]],
                                                    axis=0).sort_values('REG')], axis=1)

        # TEMPORARY DATAFRAME FOR ALL CONTRIBUTIONS MAIN CONTRIBUTIONS
        if DISPERSION:
            df_temp = pd.concat([df_intra, df_inter, df_disp]).sort_values('REG').reset_index()
        else:
            df_temp = pd.concat([df_intra, df_inter]).sort_values('REG').reset_index()
        df_temp['TERM'] = df_temp['TERM'].str.replace("VC_IQA\(A,B\)-", 'Vcl(')
        df_temp['TERM'] = df_temp['TERM'].str.replace("VX_IQA\(A,B\)-", 'Vxc(')
        df_temp['TERM'] = df_temp['TERM'].str.replace("E_IQA_Intra\(A\)-", 'Eintra(')
        df_temp['TERM'] = df_temp['TERM'].str.replace("_", ',')
        df_temp['TERM'] = df_temp['TERM'] + ')'
        dataframe_list.append(df_temp)

        # FILTERCONTRIBUTIONS: Vxc
        col = ['TERM', 'REG', 'R']
        df_xc = pd.DataFrame()
        for j in range(0, len(df_inter['TERM'])):
            if 'VX_IQA' in df_inter['TERM'][j]:
                df_xc = df_xc.append(df_inter.iloc[j])
        df_xc = df_xc[col].reset_index()
        df_xc = df_xc.sort_values('REG').reset_index()
        df_xc['TERM'] = df_xc['TERM'].str.replace("VX_IQA\(A,B\)-", 'Vxc(')
        df_xc['TERM'] = df_xc['TERM'].str.replace("_", ',')
        df_xc['TERM'] = df_xc['TERM'] + ')'

        df_xc.to_csv("Vxc_seg_" + str(i + 1) + ".csv", sep=',')
        df_xc_sorted = pd.concat(
            [df_xc_sorted, pd.concat([df_xc[-n_terms:], df_xc[:n_terms]], axis=0).sort_values('REG')], axis=1)

        # FILTER CONTRIBUTIONS: Vcl

        df_cl = pd.DataFrame()
        for j in range(0, len(df_inter['TERM'])):
            if 'VC_IQA' in df_inter['TERM'][j]:
                df_cl = df_cl.append(df_inter.iloc[j])
        df_cl = df_cl[col].reset_index()
        df_cl = df_cl.sort_values('REG').reset_index()
        df_cl['TERM'] = df_cl['TERM'].str.replace("VC_IQA\(A,B\)-", 'Vcl(')
        df_cl['TERM'] = df_cl['TERM'].str.replace("_", ',')
        df_cl['TERM'] = df_cl['TERM'] + ')'
        df_cl.to_csv("Vcl_seg_" + str(i + 1) + ".csv", sep=',')
        df_cl_sorted = pd.concat(
            [df_cl_sorted, pd.concat([df_cl[-n_terms:], df_cl[:n_terms]], axis=0).sort_values('REG')], axis=1)

        # FILTER CONTRIBUTIONS: Eintra
        df_intra = df_intra.sort_values('REG').reset_index()
        df_intra['TERM'] = df_intra['TERM'].str.replace("E_IQA_Intra\(A\)-", 'Eintra(')
        df_intra['TERM'] = df_intra['TERM'].str.replace("_", ',')
        df_intra['TERM'] = df_intra['TERM'] + ')'
        if len(atoms) > n_terms:
            df_intra_sorted = pd.concat(
                [df_intra_sorted, pd.concat([df_intra[-n_terms:], df_intra[:n_terms]], axis=0).sort_values('REG')],
                axis=1)
        else:
            df_intra_sorted = pd.concat([df_intra_sorted, pd.concat([df_intra], axis=0).sort_values('REG')], axis=1)

        df_final = pd.concat([df_final, pd.concat([df_temp[-n_terms:], df_temp[:n_terms]], axis=0).sort_values('REG')],
                             axis=1)

    # SAVE .csv files and tables .png
    df_final.to_csv('IQA_segment_analysis.csv', sep=',')
    rv.pandas_REG_dataframe_to_table(df_final, 'REG_final_table', SAVE_FIG=SAVE_FIG)
    df_cl_sorted.to_csv('IQA_segment_Vcl_analysis.csv', sep=',')
    rv.pandas_REG_dataframe_to_table(df_cl_sorted, 'REG_Vcl_table', SAVE_FIG=SAVE_FIG)
    df_xc_sorted.to_csv('IQA_segment_Vxc_analysis.csv', sep=',')
    rv.pandas_REG_dataframe_to_table(df_xc_sorted, 'REG_Vxc_table', SAVE_FIG=SAVE_FIG)
    df_intra_sorted.to_csv('IQA_segment_Eintra_analysis.csv', sep=',')
    rv.pandas_REG_dataframe_to_table(df_intra_sorted, 'REG_Eintra_table', SAVE_FIG=SAVE_FIG)
    df_inter_full_sorted.to_csv('IQA_segment_full-inter_analysis.csv', sep=',')
    rv.pandas_REG_dataframe_to_table(df_inter_full_sorted, 'REG_E_inter_full_table', SAVE_FIG=SAVE_FIG)
    if CHARGE_TRANSFER_POLARISATION:
        df_ct_pl_sorted.to_csv('IQA_segment_Vcl-Vpl_analysis.csv')
        rv.pandas_REG_dataframe_to_table(df_ct_pl_sorted, 'REG_Vct-Vpl_table', SAVE_FIG=SAVE_FIG)
    if DISPERSION:
        df_dispersion_sorted.to_csv('IQA_segment_Vdisp_analysis.csv', sep=',')
        rv.pandas_REG_dataframe_to_table(df_dispersion_sorted, 'REG_Vdisp_table', SAVE_FIG=SAVE_FIG)

    # ENERGY ONLY .CSV FILES

    df_energy_output = pd.DataFrame()
    df_energy_output['WFN'] = total_energy_wfn
    df_energy_output['IQA'] = total_energy_iqa
    df_energy_output.to_csv('total_energy.csv', sep=',')

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

rv.plot_segment(cc, 2625.50 * (total_energy_wfn - total_energy_wfn[0]), critical_points, annotate=ANNOTATE,
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
