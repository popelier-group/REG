# IMPORT MODULES
import reg 
import aimall_utils
import numpy as np
import pandas as pd

##############################    OPTIONS    ##################################
AA = 'arg'
ANGLE = 'phi'
POINTS =5
INFLEX = True
WRITE = False
AUTO = False
tp = [6, 9, 16]
ANNOTATE = True
###############################################################################

###############################################################################
#                                                                             #
#                               PREPARE FILES                                 #
#                                                                             #
###############################################################################

#DEFINE PATHS AND FILES:
path = "../" + AA +"/" + ANGLE + "/"
A=[0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23]
wfn_files = [path + str(format(i, '03d')) + '/' + str(format(i, '03d')) + '.wfn' for i in A]
folders = [path + str(format(i, '03d')) +'/'+ str(format(i, '03d')) +'_atomicfiles' for i in A]

#GET ATOM LIST FROM ANY .WFN FILE:
atoms = aimall_utils.get_atom_list(wfn_files[0])

###############################################################################
#   Note: In this case we decided to group some atoms in specific groups.     #
#   All atoms from the main chain are put into one of the first seven groups. # 
#   Other atoms, thoose belonging to the side chain, are put into G8          #
###############################################################################

groups = [['h1', 'h2', 'h3', 'c4'],     #G1
          ['c5', 'o6'],                 #G2
          ['n7', 'h8'],                 #G3
          ['c9','h10'],                 #G4
          ['c11', 'o12'],               #G5
          ['n13', 'h14'],               #G6
          ['c15', 'h16', 'h17', 'h18'], #G7
          []]                           #G8

#Put all other atoms in G8:
[groups[7].append(atoms[i]) for i in range(len(atoms))]
[groups[7].remove(groups[j][i]) for j in range(0,7) for i in range(len(groups[j]))]

###############################################################################

#DEFINE THE DESIRED TERMS:
intra_prop = ['E_IQA_Intra(A)'] #intra atomic properties
inter_prop = ['VC_IQA(A,B)', 'VX_IQA(A,B)'] #inter atomic propertier

#DEFINE THE CONTROL COORDINATES:
cc_original = [i for i in range (-180, 180, 15)]
cc = [i for i in range(0, 360, 15)]
#GET TOTAL ENERGY FROM THE .WFN FILES:
total_energy_wfn = aimall_utils.get_aimall_wfn_energies(wfn_files)

#GET INTRA-ATOMIC TERMS:
iqa_intra = aimall_utils.intra_property_from_int_file(folders,intra_prop, atoms)[0]
iqa_intra_header = np.array(aimall_utils.intra_property_from_int_file(folders,intra_prop, atoms)[1]) #used for reference
#GET INTER-ATOMIC TERMS:
iqa_inter = aimall_utils.inter_property_from_int_file(folders, inter_prop, atoms)[0]
iqa_inter_header = np.array(aimall_utils.inter_property_from_int_file(folders,inter_prop, atoms)[1]) #used for refenrence


###############################################################################
#   Note: In this case, the control coordinate is cyclic, because of that is  #
#   convinient to star from the energy minima. Doing that, the global minima  #
#   will allways be the first point of the first segment.                     #
###############################################################################

while total_energy_wfn[0] != min(total_energy_wfn):
    total_energy_wfn.append(total_energy_wfn.pop(0))
    cc_original.append(cc_original.pop(0))
    for i in range(len(iqa_intra)):
        iqa_intra[i].append(iqa_intra[i].pop(0))      
    for i in range(len(iqa_inter)):
        iqa_inter[i].append(iqa_inter[i].pop(0))

###############################################################################

#CONVERT LISTS TO ARRAYS.
total_energy_wfn = np.array(total_energy_wfn)
iqa_inter = np.array(iqa_inter)
iqa_intra = np.array(iqa_intra)

#CALCULATED THE IQA TOTAL ENERGY 
total_energy_iqa = sum(iqa_inter) + sum(iqa_intra) # used to calculate the integration error.

#GROUP INTRA ATOMIC TERMS:
grouped_iqa_intra = np.array(reg.group_intra(iqa_intra, iqa_intra_header, groups)[0])
grouped_iqa_intra_header = np.array(reg.group_intra(iqa_intra, iqa_intra_header, groups)[1])

#GROUP INTER ATOMIC TERMS:
grouped_iqa_inter = np.array(reg.group_inter(iqa_inter, iqa_inter_header, groups)[0])
grouped_iqa_inter_header = np.array(reg.group_inter(iqa_inter, iqa_inter_header, groups)[1])

###############################################################################
#                                                                             #
#                               REG ANALYSIS                                  #
#                                                                             #
###############################################################################

#INTRA ATOMIC CONTRIBUTION
reg_intra = reg.reg(total_energy_wfn, cc, iqa_intra, np=POINTS, critical=AUTO, inflex=INFLEX, critical_index=tp)
#INTER ATOMIC CONTRIBUTION
reg_inter = reg.reg(total_energy_wfn, cc, iqa_inter, np=POINTS, critical=AUTO, inflex=INFLEX, critical_index=tp)

#GROUPED INTRA ATOMIC CONTRIBUTION
grouped_reg_intra = reg.reg(total_energy_wfn, cc, grouped_iqa_intra, np=POINTS, critical=AUTO, inflex=INFLEX, critical_index=tp)
#GROUPED INTER ATOMIC CONTRIBUTION
grouped_reg_inter = reg.reg(total_energy_wfn, cc, grouped_iqa_inter, np=POINTS, critical=AUTO, inflex=INFLEX, critical_index=tp)

#CALCULATE INTEGRATION ERROR:
rmse_integration = reg.integration_error(total_energy_wfn, total_energy_iqa)
print('Integration error [kJ/mol](RMSE)')
print(rmse_integration[1])

###############################################################################
#                                                                             #
#                             WRITE CSV FILES                                 #
#                                                                             #
###############################################################################
if WRITE == True:  
    #IQA contributions:
    df_final = pd.DataFrame()  
    for i in range(len(reg_inter[0])):
        temp = [reg_inter[0][i], reg_inter[1][i]]                                  #inter
        df_inter = pd.DataFrame(temp).transpose()                                  #inter
        df_inter.columns = ["REG", "R"]                                            #inter
        df_inter.index = iqa_inter_header                                          #inter 
        df_inter.to_csv(path + "Inter seg_" +str(i+1) +".csv", sep=',')            #inter
        df_inter = df_inter.rename_axis('TERM').reset_index()                      #inter
        temp = [reg_intra[0][i], reg_intra[1][i]]                                  #intra
        df_intra = pd.DataFrame(temp).transpose()                                  #intra
        df_intra.columns = ["REG", "R"]                                            #intra
        df_intra.index = iqa_intra_header                                          #intra
        df_intra.to_csv(path + "Intra seg_" +str(i+1) + ".csv", sep=',')           #intra
        df_intra = df_intra.rename_axis('TERM').reset_index()                      #intra
    
        df_temp = pd.concat([df_intra, df_inter]).sort_values('REG').reset_index()
        df_final=pd.concat([df_final, pd.concat([df_temp[-5:], df_temp[:5]], axis=0).sort_values('REG', ascending=False)], axis=1)
    
    df_final.to_csv(path + 'IQA_segment_analysis.csv', sep=',')
    
    #Grouped IQA contributions:
    df_final = pd.DataFrame()   
    for i in range(len(grouped_reg_inter[0])):                                 
        temp = [grouped_reg_inter[0][i], grouped_reg_inter[1][i]]                  #inter
        df_ginter = pd.DataFrame(temp).transpose()                                 #inter
        df_ginter.columns = ["REG", "R"]                                           #inter
        df_ginter.index = grouped_iqa_inter_header                                 #inter
        df_ginter.to_csv(path + "grouped_Inter seg_" +str(i+1) +".csv", sep=',')   #inter 
        df_ginter = df_ginter.rename_axis('TERM').reset_index()                    #inter
        temp = [grouped_reg_intra[0][i], grouped_reg_intra[1][i]]                  #intra
        df_gintra = pd.DataFrame(temp).transpose()                                 #intra
        df_gintra.columns = ["REG", "R"]                                           #intra
        df_gintra.index = grouped_iqa_intra_header                                 #intra
        df_gintra.to_csv(path + "grouped_Intra seg_" +str(i+1) + ".csv", sep=',')  #intra
        df_gintra = df_gintra.rename_axis('TERM').reset_index()                    #intra
    
        df_temp = pd.concat([df_gintra, df_ginter]).sort_values('REG').reset_index()
        df_final=pd.concat([df_final, pd.concat([df_temp[-5:], df_temp[:5]], axis=0).sort_values('REG', ascending=False)], axis=1)
    
    df_final.to_csv(path + 'Grouped_IQA_segment_analysis.csv', sep=',')

  
