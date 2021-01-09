# IMPORT MODULES
import reg 
import aimall_utils
import numpy as np
import pandas as pd
import reg_vis as rv

##############################    OPTIONS    ##################################
SYS = 'bbr3'
EQ_DIST = 0
POINTS = 3
INFLEX = False
WRITE = True
AUTO = False
tp = []

SAVE_FIG = True
ANNOTATE = True
DETAILED_ANALYSIS = True
LABELS= True
n_terms = 4
###############################################################################

###############################################################################
#                                                                             #
#                                 SETUP FILES                                 #
#                                                                             #
###############################################################################

#DEFINE PATHS AND FILES:
path = "../rehybridization/" + SYS +'/'
A=[i for i in range(1,12)]

#DEFINE THE CONTROL COORDINATES:
A=[i for i in range(1,12)]

wfn_files = [path + str(format(i, '02d')) + '/' + str(format(i, '02d')) + '.wfn' for i in A]
folders = [path + str(format(i, '02d')) +'/'+ str(format(i, '02d')) +'_atomicfiles' for i in A]


###############################################################################
#   Note: In this case we decided to group some atoms in specific groups.     #
###############################################################################

groups = [['b1'],                 #G1 
          ['br2', 'br3', 'br4']]     #G2     

#GET ATOM LIST FROM ANY .WFN FILE:
atoms = aimall_utils.get_atom_list(wfn_files[0])

#DEFINE THE DESIRED TERMS:
intra_prop = ['E_IQA_Intra(A)'] #intra atomic properties
inter_prop = ['VC_IQA(A,B)', 'VX_IQA(A,B)'] #inter atomic propertier

#GET TOTAL ENERGY FROM THE .WFN FILES:
total_energy_wfn = aimall_utils.get_aimall_wfn_energies(wfn_files)

#GET INTRA-ATOMIC TERMS:
iqa_intra_temp = aimall_utils.intra_property_from_int_file(folders,intra_prop, atoms)[0]
iqa_intra_header_temp = np.array(aimall_utils.intra_property_from_int_file(folders,intra_prop, atoms)[1]) #used for reference
#GET INTER-ATOMIC TERMS:
iqa_inter_temp = aimall_utils.inter_property_from_int_file(folders, inter_prop, atoms)[0]
iqa_inter_header_temp = np.array(aimall_utils.inter_property_from_int_file(folders,inter_prop, atoms)[1]) #used for refenrence

#CONVERT LISTS TO ARRAYS.
total_energy_wfn = np.array(total_energy_wfn)
iqa_inter_temp = np.array(iqa_inter_temp)
iqa_intra_temp = np.array(iqa_intra_temp)


#CALCULATED THE IQA TOTAL ENERGY 
total_energy_iqa = sum(iqa_inter_temp) + sum(iqa_intra_temp) # used to calculate the integration error.

#GROUP INTRA ATOMIC TERMS:
iqa_intra = np.array(reg.group_intra(iqa_intra_temp, iqa_intra_header_temp, groups)[0])
iqa_intra_header = np.array(reg.group_intra(iqa_intra_temp, iqa_intra_header_temp, groups)[1])

#GROUP INTER ATOMIC TERMS:
iqa_inter = np.array(reg.group_inter(iqa_inter_temp, iqa_inter_header_temp, groups)[0])
iqa_inter_header = np.array(reg.group_inter(iqa_inter_temp, iqa_inter_header_temp, groups)[1])

###############################################################################
#                                                                             #
#                               REG ANALYSIS                                  #
#                                                                             #
###############################################################################

#INTRA ATOMIC CONTRIBUTION
reg_intra = reg.reg(total_energy_wfn, cc, iqa_intra, np=POINTS, critical=AUTO, inflex=INFLEX, critical_index=tp)
#INTER ATOMIC CONTRIBUTION
reg_inter = reg.reg(total_energy_wfn, cc, iqa_inter, np=POINTS, critical=AUTO, inflex=INFLEX, critical_index=tp)

#CALCULATE INTEGRATION (RECOVERY) ERROR:
rmse_integration = reg.integration_error(total_energy_wfn, total_energy_iqa)
print('Integration error [kJ/mol](RMSE)')
print(rmse_integration[1])

###############################################################################
#                                                                             #
#                             WRITE CSV FILES                                 #
#                                                                             #
###############################################################################
dataframe_list=[]

if WRITE == True:  
     df_final = pd.DataFrame()   
     for i in range(len(reg_inter[0])):                                 
        temp = [reg_inter[0][i], reg_inter[1][i]]                  #inter
        df_inter = pd.DataFrame(temp).transpose()                                 #inter
        df_inter.columns = ["REG", "R"]                                           #inter
        df_inter.index = iqa_inter_header                                 #inter
        df_inter.to_csv(path + "Inter_seg_" +str(i+1) +".csv", sep=',')   #inter 
        df_inter = df_inter.rename_axis('TERM').reset_index()                    #inter
        temp = [reg_intra[0][i], reg_intra[1][i]]                  #intra
        df_intra = pd.DataFrame(temp).transpose()                                 #intra
        df_intra.columns = ["REG", "R"]                                           #intra
        df_intra.index = iqa_intra_header                                 #intra
        df_intra.to_csv(path + "Intra_seg_" +str(i+1) + ".csv", sep=',')  #intra
        df_intra = df_intra.rename_axis('TERM').reset_index()                    #intra
    
        df_temp = pd.concat([df_intra, df_inter]).sort_values('REG').reset_index()

        df_temp['TERM'] = df_temp['TERM'].str.replace("G1", 'B')       
        df_temp['TERM'] = df_temp['TERM'].str.replace("G2", 'X')
        df_temp['TERM'] = df_temp['TERM'].str.replace("G3", 'N')
        df_temp['TERM'] = df_temp['TERM'].str.replace("G4", 'H')
        
        df_temp['TERM'] = df_temp['TERM'].str.replace("VC_IQA\(A,B\)-", 'Vcl(')
        df_temp['TERM'] = df_temp['TERM'].str.replace("VX_IQA\(A,B\)-", 'Vxc(')
        df_temp['TERM'] = df_temp['TERM'].str.replace("E_IQA_Intra\(A\)-", 'Eintra(')   
        df_temp['TERM'] = df_temp['TERM'].str.replace("_", ',')
        df_temp['TERM'] = df_temp['TERM'] + ')'
        df_temp['R'] = df_temp['R']**2
        dataframe_list.append(df_temp)
        df_final=pd.concat([df_final, pd.concat([df_temp[-5:], df_temp[:5]], axis=0).sort_values('REG', ascending=False)], axis=1)
        
     df_final.to_csv(path + 'IQA_segment_analysis.csv', sep=',')

     rv.plot_violin([dataframe_list[i]['R'] for i in range(len(reg_inter[0]))], save=SAVE_FIG, file_name = path +'violin.png' )


###############################################################################
#                                                                             #
#                                   GRAPHS                                    #
#                                                                             #
###############################################################################

if AUTO == True:
    critical_points = reg.find_critical(total_energy_wfn, cc, min_points=POINTS, use_inflex=INFLEX)
elif AUTO == False:
    critical_points = tp

rv.plot_segment(cc, 2625.50*total_energy_wfn, critical_points, annotate=ANNOTATE, label=LABELS,
                y_label = r'Energy [$kJ.mol^{-1}$]', x_label=r'r (a.u)', title = SYS, save = SAVE_FIG, file_name=path +'REG_analysis.png')       

if DETAILED_ANALYSIS ==True:
    for i in range(len(reg_inter[0])):
        rv.generate_data_vis(dataframe_list[i], [dataframe_list[i]['R'] for i in range(len(reg_inter[0]))],
                                            n_terms, save=SAVE_FIG, file_name=path +'detailed_seg_'+ str(i+1) +'.png', title = SYS +' seg. '+str(i+1))

##### save data ####

intra_properties = pd.DataFrame(iqa_intra*2625.5).T
intra_properties.columns = iqa_intra_header
intra_properties.columns = intra_properties.columns.str.replace("E_IQA_Intra\(A\)-", 'Eintra(')
intra_properties.columns = intra_properties.columns.str.replace("G1", 'B')
intra_properties.columns = intra_properties.columns.str.replace("G2", 'X')
intra_properties.columns = intra_properties.columns.str.replace("G3", 'N')
intra_properties.columns = intra_properties.columns.str.replace("G4", 'H')

intra_properties.columns = intra_properties.columns + ')'

intra_properties.to_csv(path+'intra_properties.csv')

inter_properties = pd.DataFrame(iqa_inter*2625.5).T
inter_properties.columns = iqa_inter_header
inter_properties.columns = inter_properties.columns.str.replace("VC_IQA\(A,B\)-", 'Vcl(')
inter_properties.columns = inter_properties.columns.str.replace("VX_IQA\(A,B\)-", 'Vxc(')
inter_properties.columns = inter_properties.columns.str.replace("G1", 'B')
inter_properties.columns = inter_properties.columns.str.replace("G2", 'X')
inter_properties.columns = inter_properties.columns.str.replace("G3", 'N')
inter_properties.columns = inter_properties.columns.str.replace("G4", 'H')
inter_properties.columns = inter_properties.columns.str.replace("_", ',')
inter_properties.columns = inter_properties.columns + ')'

inter_properties.to_csv(path+'inter_properties.csv')

