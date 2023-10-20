######GUI INITIATION########
import PySimpleGUI as sg

sg.theme("TealMono")

###FILE SEARCH#####

filepath = sg.popup_get_folder("Enter REG.py installation directory")
working_folder = sg.popup_get_folder("Enter working folder for the analysis")
control_coordinate_file = sg.popup_get_file("Enter control coordinate file")

#####ANALYSIS INPUT VARIABLES########
layout = [
    [
        sg.Checkbox("Auto Critical Point Search", key="-AUTO-"),
        sg.InputText("Enter Number of minimum points for search", key="-POINTS-"),
    ],
    [
        sg.Checkbox("Inflex", key="-INFLEX-"),
        sg.Checkbox("Write CSV Files", key="-WRITE-"),
    ],
    [
        sg.Checkbox("Save Figures", key="-SAVE_FIG-"),
        sg.Checkbox("Annotate Figures", key="-ANNOTATE-"),
        sg.Checkbox("Label Figures", key="-LABELS-"),
    ],
    [
        sg.Checkbox("Detailed Analysis", key="-DETAILED_ANALYSIS-"),
        sg.Text("N Terms Analysis"),
        sg.Input(size=(2, 2), key="-N_TERMS-"),
    ],
    [sg.Button("Done")],
]

# Create the window
window = sg.Window("REG Analysis Input", layout)

# Create an event loop
while True:
    event, values = window.read()
    # End program if user closes window or
    # presses the Close button
    if event == "Done" or event == sg.WIN_CLOSED:
        break

window.close()


############################

# IMPORT MODULES
import sys

sys.path.insert(1, filepath)
import reg
import aimall_utils
import numpy as np
import pandas as pd
import reg_vis as rv
import seaborn as sns
import gaussian_utils
import matplotlib.pyplot as plt
import os

##############################    VARIABLES/OPTIONS    ##################################
# SYS = values['-SYSTEM_NAME-']
POINTS = int(values["-POINTS-"])
INFLEX = values["-INFLEX-"]
WRITE = values["-WRITE-"]
AUTO = values["-AUTO-"]
if AUTO == False:
    quickLayout = [
        [sg.Text("Insert critical points positions separated by a comma")],
        [sg.InputText()],
        [sg.Button("Done")],
    ]
    popup_auto = sg.Window("Critical Points", quickLayout)
    while True:
        event, critical_points = popup_auto.read()
        if event == "Done" or event == sg.WIN_CLOSED:
            break
    popup_auto.close()
    tp = [int(i) for i in critical_points[0].split(",")]

SAVE_FIG = values["-SAVE_FIG-"]
ANNOTATE = values["-ANNOTATE-"]
DETAILED_ANALYSIS = values["-DETAILED_ANALYSIS-"]
LABELS = values["-LABELS-"]
n_terms = int(values["-N_TERMS-"])


###############################################################################
#                                                                             #
#                                 SETUP FILES                                 #
#                                                                             #
###############################################################################


# DEFINE PATHS AND FILES:
def skim_list(index_list, A):
    """
    Function to skim the lists to test REG analysis with fewer points along the control coordinate
    """
    skimmed_list = []
    try:
        for i in index_list:
            skimmed_list.append(A[i - 1])
    except:
        pass
    return skimmed_list


os.chdir(working_folder)
cc = gaussian_utils.get_control_coordinates_IRC_g16(control_coordinate_file)
A = [i for i in range(1, len(cc) + 1)]

wfn_files = [working_folder + "/" + str(i) + "/frame.wfn" for i in A]
folders = [working_folder + "/" + str(i) + "/frame_atomicfiles/" for i in A]

# GET ATOM LIST FROM ANY .WFN FILE:
atoms = aimall_utils.get_atom_list(wfn_files[0])

# DEFINE THE DESIRED TERMS:
intra_prop = ["E_IQA_Intra(A)"]  # intra atomic properties
inter_prop = ["VC_IQA(A,B)", "VX_IQA(A,B)"]  # inter atomic properties

# GET TOTAL ENERGY FROM THE .WFN FILES:
total_energy_wfn = aimall_utils.get_aimall_wfn_energies(wfn_files)

# GET INTRA-ATOMIC TERMS:
iqa_intra_temp = aimall_utils.intra_property_from_int_file(folders, intra_prop, atoms)[
    0
]
iqa_intra_header_temp = np.array(
    aimall_utils.intra_property_from_int_file(folders, intra_prop, atoms)[1]
)  # used for reference
# GET INTER-ATOMIC TERMS:
iqa_inter_temp = aimall_utils.inter_property_from_int_file(folders, inter_prop, atoms)[
    0
]
iqa_inter_header_temp = np.array(
    aimall_utils.inter_property_from_int_file(folders, inter_prop, atoms)[1]
)  # used for reference

# CONVERT LISTS TO ARRAYS.
total_energy_wfn = np.array(total_energy_wfn)
iqa_inter_temp = np.array(iqa_inter_temp)
iqa_intra_temp = np.array(iqa_intra_temp)

# CALCULATED THE IQA TOTAL ENERGY
total_energy_iqa = sum(iqa_inter_temp) + sum(
    iqa_intra_temp
)  # used to calculate the integration error.

###############################################################################
#                                                                             #
#                               REG ANALYSIS                                  #
#                                                                             #
###############################################################################

# INTRA ATOMIC CONTRIBUTION
reg_intra = reg.reg(
    total_energy_wfn, cc, iqa_intra_temp, np=POINTS, critical=AUTO, inflex=INFLEX
)
# INTER ATOMIC CONTRIBUTION
reg_inter = reg.reg(
    total_energy_wfn, cc, iqa_inter_temp, np=POINTS, critical=AUTO, inflex=INFLEX
)

# CALCULATE INTEGRATION (RECOVERY) ERROR:
rmse_integration = reg.integration_error(total_energy_wfn, total_energy_iqa)
print("Recovery Error [kJ/mol](RMSE)")
print(rmse_integration[1])

###############################################################################
#                                                                             #
#                             WRITE CSV FILES                                 #
#                                                                             #
###############################################################################
dataframe_list = []

if WRITE == True:
    df_final = pd.DataFrame()
    for i in range(len(reg_inter[0])):
        temp = [reg_inter[0][i], reg_inter[1][i]]
        df_inter = pd.DataFrame(temp).transpose()
        df_inter.columns = ["REG", "R"]
        df_inter.index = iqa_inter_header_temp
        df_inter.to_csv(
            working_folder + "/Inter_seg_" + str(i + 1) + ".csv", sep=","
        )  # inter
        df_inter = df_inter.rename_axis("TERM").reset_index()
        temp = [reg_intra[0][i], reg_intra[1][i]]  # intra
        df_intra = pd.DataFrame(temp).transpose()  # intra
        df_intra.columns = ["REG", "R"]  # intra
        df_intra.index = iqa_intra_header_temp  # intra
        df_intra.to_csv(
            working_folder + "/Intra_seg_" + str(i + 1) + ".csv", sep=","
        )  # intra
        df_intra = df_intra.rename_axis("TERM").reset_index()  # intra

        df_temp = pd.concat([df_intra, df_inter]).sort_values("REG").reset_index()

        df_temp["TERM"] = df_temp["TERM"].str.replace("G1", "B")
        df_temp["TERM"] = df_temp["TERM"].str.replace("G2", "X")
        df_temp["TERM"] = df_temp["TERM"].str.replace("G3", "N")
        df_temp["TERM"] = df_temp["TERM"].str.replace("G4", "H")

        df_temp["TERM"] = df_temp["TERM"].str.replace("VC_IQA\(A,B\)-", "Vcl(")
        df_temp["TERM"] = df_temp["TERM"].str.replace("VX_IQA\(A,B\)-", "Vxc(")
        df_temp["TERM"] = df_temp["TERM"].str.replace("E_IQA_Intra\(A\)-", "Eintra(")
        df_temp["TERM"] = df_temp["TERM"].str.replace("_", ",")
        df_temp["TERM"] = df_temp["TERM"] + ")"
        df_temp["R"] = df_temp["R"] ** 2
        dataframe_list.append(df_temp)
        df_final = pd.concat(
            [
                df_final,
                pd.concat([df_temp[-n_terms:], df_temp[:n_terms]], axis=0).sort_values(
                    "REG", ascending=False
                ),
            ],
            axis=1,
        )
    df_final.to_csv(working_folder + "/IQA_segment_analysis.csv", sep=",")

# rv.plot_violin([dataframe_list[i]['R'] for i in range(len(reg_inter[0]))], save=SAVE_FIG, file_name = working_folder +'/violin.png' )

###############################################################################
#                                                                             #
#                                   GRAPHS                                    #
#                                                                             #
###############################################################################

if AUTO == True:
    critical_points = reg.find_critical(
        total_energy_wfn, cc, min_points=POINTS, use_inflex=INFLEX
    )
elif AUTO == False:
    critical_points = tp

rv.plot_segment(
    cc,
    2625.50 * total_energy_wfn,
    critical_points,
    annotate=ANNOTATE,
    label=LABELS,
    y_label=r"Energy [$kJ.mol^{-1}$]",
    x_label=r"r (a.u)",
    title=working_folder.split("/")[-1],
    save=SAVE_FIG,
    file_name=working_folder + "/REG_analysis.png",
)

if DETAILED_ANALYSIS == True:
    for i in range(len(reg_inter[0])):
        rv.generate_data_vis(
            dataframe_list[i],
            [dataframe_list[i]["R"] for i in range(len(reg_inter[0]))],
            n_terms,
            save=SAVE_FIG,
            file_name=working_folder + "/detailed_seg_" + str(i + 1) + ".png",
            title=working_folder.split("/")[-1] + " seg. " + str(i + 1),
        )

##### save data ####
intra_properties = pd.DataFrame(iqa_intra_temp * 2625.5).T
intra_properties.columns = iqa_intra_header_temp
intra_properties.columns = intra_properties.columns.str.replace(
    "E_IQA_Intra\(A\)-", "Eintra("
)
intra_properties.columns = intra_properties.columns.str.replace("G1", "B")
intra_properties.columns = intra_properties.columns.str.replace("G2", "X")
intra_properties.columns = intra_properties.columns.str.replace("G3", "N")
intra_properties.columns = intra_properties.columns.str.replace("G4", "H")

intra_properties.columns = intra_properties.columns + ")"

intra_properties.to_csv(working_folder + "/intra_properties.csv")

inter_properties = pd.DataFrame(iqa_inter_temp * 2625.5).T
inter_properties.columns = iqa_inter_header_temp
inter_properties.columns = inter_properties.columns.str.replace(
    "VC_IQA\(A,B\)-", "Vcl("
)
inter_properties.columns = inter_properties.columns.str.replace(
    "VX_IQA\(A,B\)-", "Vxc("
)
inter_properties.columns = inter_properties.columns.str.replace("G1", "B")
inter_properties.columns = inter_properties.columns.str.replace("G2", "X")
inter_properties.columns = inter_properties.columns.str.replace("G3", "N")
inter_properties.columns = inter_properties.columns.str.replace("G4", "H")
inter_properties.columns = inter_properties.columns.str.replace("_", ",")
inter_properties.columns = inter_properties.columns + ")"

inter_properties.to_csv(working_folder + "/inter_properties.csv")
