#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Mar  2 11:56:01 2021

@author: Ezan
"""

# IMPORT MODULES
import reg 
import aimall_utils
import numpy as np
#import numpy.polynomial.polynomial as poly
#import matplotlib.pyplot as plt
import pandas as pd
import reg_vis as rv

##############################    OPTIONS    ##################################
SYS0 = 'BH2'
SYS1 = 'BH2-out'
SYS2 = 'BH2-in'
SYS_f1 = 'BH2-fit-out'
SYS_f2 = 'BH2-fit-in'
SYS = 'BH2-diff (in-out)'
path='BH2-diff'

EQ_DIST = 0
POINTS = 4
INFLEX = True
WRITE = True
AUTO = True
AUTO_fit = False 
FORCE = False #If you want the force, or not
COLOR_TERMS = False # If you want your term names to be colored in the segment analysis, or not
#disp_lim=[1,8,20]

SIZE_TITLE = 18
SIZE_LABEL = 16 
SAVE_FIG = True
ANNOTATE = True
DETAILED_ANALYSIS = True
LABELS= True
n_terms = 5 #number of terms to be displayed on each segment graph

###############################################################################
#                                                                             #
#                                 SETUP FILES                                 #
#                                                                             #
###############################################################################

#DEFINE PATHS AND FILES FOR PATH 1:
#path1 = "../BH2/SPout" #Path to the folder with .wfn 
path1 = "BH2-out" #Path to the folder with .wfn 
F1=[i for i in range(0,155,5)]
R1=[i for i in range(5,70,5)]

#DEFINE PATHS AND FILES FOR PATH 2:
#path2 = "../BH2/SPin" #Path to the folder with .wfn 
path2 = "BH2-in" #Path to the folder with .wfn 
F2=[i for i in range(0,130,5)]
R2=[i for i in range(5,65,5)]

#COMMON

#GET FILES AND CONTROL COORDINATES FOR PATH 1:
cc_a1=aimall_utils.get_list_wfx(path1,SYS0,"out",F1,R1,False)[0]
folders1=aimall_utils.get_list_wfx(path1,SYS0,"out",F1,R1,False)[1]
folders_disp1=aimall_utils.get_list_wfx(path1,SYS0,"out",F1,R1,False)[2]
wfx_files1=aimall_utils.get_list_wfx(path1,SYS0,"out",F1,R1,False)[3]

cc_a2=aimall_utils.get_list_wfx(path2,SYS0,"in",F2,R2,False)[0]
folders2=aimall_utils.get_list_wfx(path2,SYS0,"in",F2,R2,False)[1]
folders_disp2=aimall_utils.get_list_wfx(path2,SYS0,"in",F2,R2,False)[2]
wfx_files2=aimall_utils.get_list_wfx(path2,SYS0,"in",F2,R2,False)[3]

#GET ATOM LIST FROM ANY .WFX FILE:
atoms1 = aimall_utils.get_atom_list_wfx(wfx_files1[0])
atoms2 = aimall_utils.get_atom_list_wfx(wfx_files2[0])

#GROUP TERMS
groups=[['b10','h11','h12']]
#for i in range(len(disp_lim)):
#    if i == len(disp_lim)-1 :
#        groups_temp=[atoms[j] for j in range(disp_lim[i]-1,len(atoms))]
#    else:
#        groups_temp=[atoms[j] for j in range(disp_lim[i]-1,disp_lim[i+1]-1)]
#    groups.append(groups_temp)

#DEFINE THE DESIRED TERMS:
intra_prop = ['E_IQA_Intra(A)'] #intra atomic properties
inter_prop = ['VC_IQA(A,B)', 'VX_IQA(A,B)'] #inter atomic properties

#GET TOTAL ENERGY FROM THE .WFN FILES:
total_energy_wfx1 = aimall_utils.get_aimall_wfx_energies(wfx_files1)
total_energy_wfx2 = aimall_utils.get_aimall_wfx_energies(wfx_files2)

#GET INTRA-ATOMIC TERMS:
iqa_intra1 = aimall_utils.intra_property_from_int_file(folders1, intra_prop, atoms1)[0]
iqa_intra_header1 = np.array(aimall_utils.intra_property_from_int_file(folders1,intra_prop, atoms1)[1]) #used for reference
iqa_intra2 = aimall_utils.intra_property_from_int_file(folders2, intra_prop, atoms2)[0]
iqa_intra_header2 = np.array(aimall_utils.intra_property_from_int_file(folders2,intra_prop, atoms2)[1]) #used for reference

#change the term names to shorten them
for i in range(len(iqa_intra_header1)):
        iqa_intra_header1[i] = iqa_intra_header1[i].replace("E_IQA_Intra(A)","E_Intra")
        iqa_intra_header2[i] = iqa_intra_header2[i].replace("E_IQA_Intra(A)","E_Intra")

#GET INTER-ATOMIC TERMS:
iqa_inter1 = aimall_utils.inter_property_from_int_file(folders1, inter_prop, atoms1)[0]
iqa_inter_header1 = np.array(aimall_utils.inter_property_from_int_file(folders1, inter_prop, atoms1)[1]) #used for reference
iqa_inter2 = aimall_utils.inter_property_from_int_file(folders2, inter_prop, atoms2)[0]
iqa_inter_header2 = np.array(aimall_utils.inter_property_from_int_file(folders2, inter_prop, atoms2)[1]) #used for reference

#change the term names to shorten them
for i in range(len(iqa_inter_header1)):
        iqa_inter_header1[i] = iqa_inter_header1[i].replace("_IQA(A,B)","")
        iqa_inter_header2[i] = iqa_inter_header2[i].replace("_IQA(A,B)","")

#GET TOTAL DISPERSION TERMS:
iqa_disp1 = aimall_utils.disp_property(folders_disp1)[0]
iqa_disp_header1 = np.array(aimall_utils.disp_property(folders_disp1)[1]) #used for reference
iqa_disp2 = aimall_utils.disp_property(folders_disp2)[0]
iqa_disp_header2 = np.array(aimall_utils.disp_property(folders_disp2)[1]) #used for reference

#CONVERT LISTS TO ARRAYS.
total_energy_wfx1 = np.array(total_energy_wfx1)
iqa_inter1 = np.array(iqa_inter1)
iqa_intra1 = np.array(iqa_intra1)
iqa_disp1 = np.array(iqa_disp1)
#print(iqa_intra1[0])
total_energy_wfx2 = np.array(total_energy_wfx2)
iqa_inter2 = np.array(iqa_inter2)
iqa_intra2 = np.array(iqa_intra2)
iqa_disp2 = np.array(iqa_disp2)
#print(iqa_intra2[0])

#NORMALIZE THE CONTROL COORDINATES
crit_list1 = reg.find_critical(total_energy_wfx1,cc_a1,use_inflex=True) 
crit_list2 = reg.find_critical(total_energy_wfx2,cc_a2,use_inflex=True)


if len(crit_list1) != len(crit_list2):
    crit_diff=abs(len(crit_list1)-len(crit_list2))
    print("The numbers of critical points differ on each pathway:")
    print("Critical point list for "+SYS1+": ",crit_list1," (length: "+str(len(crit_list1))+")")
    print("Critical point list for "+SYS2+": ",crit_list2," (length: "+str(len(crit_list2))+")")
    if crit_diff == 1:
        remove_points=input("Enter "+str(crit_diff)+" point to be removed from critical list:")
    else:    
        remove_points=input("Enter "+str(crit_diff)+" points to be removed from critical list:")
    list_remove_points=remove_points.split()
    for i in range(len(list_remove_points)):
        list_remove_points[i]=int(list_remove_points[i])
    if len(crit_list1) > len(crit_list2):
        for i in range(len(list_remove_points)): 
            crit_list1.remove(list_remove_points[i])
    if len(crit_list2) > len(crit_list1):
        for i in range(len(list_remove_points)): 
            crit_list2.remove(list_remove_points[i])
            
cc_a1_renorm = reg.normalize(cc_a1,crit_list1)
cc_a2_renorm = reg.normalize(cc_a2,crit_list2)

rv.plot_energy(total_energy_wfx1,total_energy_wfx2,cc_a1,cc_a2,crit_list1,crit_list2,SYS0,SIZE_TITLE,SIZE_LABEL,path) #1 = OUT, 2 = IN
rv.plot_renorm(total_energy_wfx1,total_energy_wfx2,cc_a1,cc_a2,cc_a1_renorm,cc_a2_renorm,crit_list1,SYS0,SIZE_TITLE,SIZE_LABEL,path) #1 = OUT, 2 = IN
"""
plot0o = plt.figure(figsize=(10,7)) #ENERGY OUT ONLY
plot1 = plot0o.add_subplot(111)
plot1.plot(cc_a1,2625.5*(total_energy_wfx1-total_energy_wfx1[0]),'bo--',label='Energy (out)')
plot1.axvline(cc_a1[0], linestyle = 'dotted', color = 'blue')
plot1.axvline(cc_a1[-1], linestyle = 'dotted', color = 'blue')
for i in crit_list1: #split the segments
    plot1.axvline(cc_a1[i], linestyle = 'dotted', color = 'blue')
plot1.legend()
plot1.set_title('Reaction coordinate renormalization',fontsize=18)
plot1.set_xlabel('Reaction coordinate',fontsize=16)
plot1.set_xlim(-6,12)
plot1.set_ylabel('Relative Energy ($kJ.mol^{-1}$)',fontsize=16)
plot1.set_ylim(-100,120)
plot1.tick_params(axis='both',labelsize=15)
plot0o.savefig('BH2-diff/'+SYS1+'_different_reaction_coordinates')

plot0i = plt.figure(figsize=(10,7)) #ENERGY IN ONLY
plot1 = plot0i.add_subplot(111)
plot1.plot(cc_a2,2625.5*(total_energy_wfx2-total_energy_wfx2[0]),'ro--',label='Energy (in)')
plot1.axvline(cc_a2[0], linestyle = 'dotted', color = 'red')
plot1.axvline(cc_a2[-1], linestyle = 'dotted', color = 'red')
for i in crit_list2: #split the segments
    plot1.axvline(cc_a2[i], linestyle = 'dotted', color = 'red')
plot1.legend()
plot1.set_title('Reaction coordinate renormalization',fontsize=18)
plot1.set_xlabel('Reaction coordinate',fontsize=16)
plot1.set_xlim(-6,12)
plot1.set_ylabel('Relative Energy ($kJ.mol^{-1}$)',fontsize=16)
plot1.set_ylim(-100,120)
plot1.tick_params(axis='both',labelsize=15)
plot0i.savefig('BH2-diff/'+SYS2+'_different_reaction_coordinates')
"""
#GENERATE THE FULL CRITICAL POINT LIST (+EXTREMA)
crit_list1_0N=[0]
crit_list1_0N += crit_list1
crit_list1_0N.append(len(cc_a1_renorm)-1)
seg_long1=[0]
for i in range(len(crit_list1_0N)-1):
    seg_long1.append(crit_list1_0N[i+1]-crit_list1_0N[i])
seg_long1.sort()
maxlong1=seg_long1[-1]

crit_list2_0N=[0]
crit_list2_0N += crit_list2
crit_list2_0N.append(len(cc_a2_renorm)-1)
seg_long2=[0]
for i in range(len(crit_list2_0N)-1):
    seg_long2.append(crit_list2_0N[i+1]-crit_list2_0N[i])
seg_long2.sort()
maxlong2=seg_long2[-1]

maxlong=max(maxlong1,maxlong2)

crit_list_fit=[(maxlong-1)*i for i in range(1,len(crit_list1)+1)]

#FITTING
total_energy_wfx1_fit = reg.fit(total_energy_wfx1, crit_list1_0N, maxlong, cc_a1_renorm)[0]
cc_a_sampled = reg.fit(total_energy_wfx1, crit_list1_0N, maxlong, cc_a1_renorm)[1]
iqa_intra1_fit = reg.fit(iqa_intra1, crit_list1_0N, maxlong, cc_a1_renorm)[0]
iqa_inter1_fit = reg.fit(iqa_inter1, crit_list1_0N, maxlong, cc_a1_renorm)[0]
iqa_disp1_fit = reg.fit(iqa_disp1, crit_list1_0N, maxlong, cc_a1_renorm)[0]
#print(iqa_intra1_fit)
total_energy_wfx2_fit = reg.fit(total_energy_wfx2, crit_list2_0N, maxlong, cc_a2_renorm)[0]
iqa_intra2_fit = reg.fit(iqa_intra2, crit_list2_0N, maxlong, cc_a2_renorm)[0]
iqa_inter2_fit = reg.fit(iqa_inter2, crit_list2_0N, maxlong, cc_a2_renorm)[0]
iqa_disp2_fit = reg.fit(iqa_disp2, crit_list2_0N, maxlong, cc_a2_renorm)[0]

rv.plot_fit(total_energy_wfx1,total_energy_wfx2,total_energy_wfx1_fit,total_energy_wfx2_fit,cc_a1_renorm,cc_a2_renorm,cc_a_sampled,
             crit_list_fit,SYS0,SIZE_TITLE,SIZE_LABEL,path)

#SUBSTRACT
total_energy_wfx = total_energy_wfx2_fit - total_energy_wfx1_fit
iqa_inter = iqa_inter2_fit - iqa_inter1_fit
iqa_intra = iqa_intra2_fit - iqa_intra1_fit
iqa_disp = iqa_disp2_fit - iqa_disp1_fit
#print(iqa_intra[0])

rv.plot_substract(total_energy_wfx1_fit,total_energy_wfx2_fit,total_energy_wfx,cc_a_sampled,crit_list_fit,
                  SYS0,SIZE_TITLE,SIZE_LABEL,path)
rv.plot_substract_term(iqa_intra1_fit,iqa_intra1,iqa_intra2_fit,iqa_intra2,iqa_intra,iqa_intra_header1,0,
                       cc_a1_renorm,cc_a2_renorm,cc_a_sampled,crit_list_fit,SYS0,SIZE_TITLE,SIZE_LABEL,path)

#CALCULATE THE IQA TOTAL ENERGY 
total_energy_iqa1 = sum(iqa_inter1) + sum(iqa_intra1) + sum(iqa_disp1) # used to calculate the integration error.
total_energy_iqa2 = sum(iqa_inter2) + sum(iqa_intra2) + sum(iqa_disp2) # used to calculate the integration error.
total_energy_iqa1_fit = reg.fit(total_energy_iqa1, crit_list1_0N, maxlong, cc_a1_renorm)[0]
total_energy_iqa2_fit = reg.fit(total_energy_iqa2, crit_list2_0N, maxlong, cc_a2_renorm)[0]
total_energy_iqa1_fit1 = sum(iqa_inter1_fit) + sum(iqa_intra1_fit) + sum(iqa_disp1_fit) # used to calculate the integration error.
total_energy_iqa2_fit1 = sum(iqa_inter2_fit) + sum(iqa_intra2_fit) + sum(iqa_disp2_fit) # used to calculate the integration error.

diff_total_energy_iqa = sum(iqa_inter) + sum(iqa_intra) + sum(iqa_disp) # used to calculate the integration error.
#print(total_energy_iqa1,total_energy_iqa2,diff_total_energy_iqa)

###############################################################################
#                                                                             #
#                               REG ANALYSIS                                  #
#                                                                             #
###############################################################################

#INTRA ATOMIC CONTRIBUTION
reg_intra1 = reg.reg(total_energy_wfx1, cc_a1, iqa_intra1, np=POINTS, critical=AUTO, inflex=INFLEX)
reg_intra2 = reg.reg(total_energy_wfx2, cc_a2, iqa_intra2, np=POINTS, critical=AUTO, inflex=INFLEX)
reg_intra1_fit = reg.reg(total_energy_wfx1_fit, cc_a_sampled, iqa_intra1_fit, np=POINTS, critical=AUTO_fit, inflex=INFLEX, critical_index=crit_list_fit)
reg_intra2_fit = reg.reg(total_energy_wfx2_fit, cc_a_sampled, iqa_intra2_fit, np=POINTS, critical=AUTO_fit, inflex=INFLEX, critical_index=crit_list_fit)
reg_intra = reg.reg(total_energy_wfx, cc_a_sampled, iqa_intra, np=POINTS, critical=AUTO_fit, inflex=INFLEX, critical_index=crit_list_fit)

#INTER ATOMIC CONTRIBUTION
reg_inter1 = reg.reg(total_energy_wfx1, cc_a1, iqa_inter1, np=POINTS, critical=AUTO, inflex=INFLEX)
reg_inter2 = reg.reg(total_energy_wfx2, cc_a2, iqa_inter2, np=POINTS, critical=AUTO, inflex=INFLEX)
reg_inter1_fit = reg.reg(total_energy_wfx1_fit, cc_a_sampled, iqa_inter1_fit, np=POINTS, critical=AUTO_fit, inflex=INFLEX, critical_index=crit_list_fit)
reg_inter2_fit = reg.reg(total_energy_wfx2_fit, cc_a_sampled, iqa_inter2_fit, np=POINTS, critical=AUTO_fit, inflex=INFLEX, critical_index=crit_list_fit)
reg_inter = reg.reg(total_energy_wfx, cc_a_sampled, iqa_inter, np=POINTS, critical=AUTO_fit, inflex=INFLEX, critical_index=crit_list_fit)

#GLOBAL DISPERSION CONTRIBUTION
reg_disp1 = reg.reg(total_energy_wfx1, cc_a1, iqa_disp1, np=POINTS, critical=AUTO, inflex=INFLEX)
reg_disp2 = reg.reg(total_energy_wfx2, cc_a2, iqa_disp2, np=POINTS, critical=AUTO, inflex=INFLEX)
reg_disp1_fit = reg.reg(total_energy_wfx1_fit, cc_a_sampled, iqa_disp1_fit, np=POINTS, critical=AUTO_fit, inflex=INFLEX, critical_index=crit_list_fit)
reg_disp2_fit = reg.reg(total_energy_wfx2_fit, cc_a_sampled, iqa_disp2_fit, np=POINTS, critical=AUTO_fit, inflex=INFLEX, critical_index=crit_list_fit)
reg_disp = reg.reg(total_energy_wfx, cc_a_sampled, iqa_disp, np=POINTS, critical=AUTO_fit, inflex=INFLEX, critical_index=crit_list_fit)

#GROUP INTER ATOMIC TERMS:
#iqa_disp_g = np.array(reg.group_inter(iqa_disp, iqa_disp_header, groups)[0])
#iqa_disp_header_g = np.array(reg.group_inter(iqa_disp, iqa_disp_header, groups)[1])

#CALCULATE INTEGRATION (RECOVERY) ERROR:
rmse_integration = reg.integration_error(total_energy_wfx1, total_energy_iqa1)
print('Integration error [kJ/mol](RMSE) - Pathway 1')
print(rmse_integration[1])
rmse_integration = reg.integration_error(total_energy_wfx2, total_energy_iqa2)
print('Integration error [kJ/mol](RMSE) - Pathway 2')
print(rmse_integration[1])
rmse_integration = reg.integration_error(total_energy_wfx1_fit, total_energy_iqa1_fit)
print('Integration error [kJ/mol](RMSE) - Pathway 1 (fitted)')
print(rmse_integration[1])
rmse_integration = reg.integration_error(total_energy_wfx2_fit, total_energy_iqa2_fit)
print('Integration error [kJ/mol](RMSE) - Pathway 2 (fitted)')
print(rmse_integration[1])
rmse_integration = reg.integration_error(total_energy_wfx1_fit, total_energy_iqa1_fit1)
print('Integration error [kJ/mol](RMSE) - Pathway 1 (fitted)')
print(rmse_integration[1])
rmse_integration = reg.integration_error(total_energy_wfx2_fit, total_energy_iqa2_fit1)
print('Integration error [kJ/mol](RMSE) - Pathway 2 (fitted)')
print(rmse_integration[1])

###############################################################################
#                                                                             #
#                     WRITE CSV FILES (WITHOUT GROUPS)                        #
#                                                                             #
###############################################################################

if WRITE == True: 
    
    df_final1=reg.csv(reg_inter1,iqa_inter_header1,reg_intra1,iqa_intra_header1,reg_disp1,iqa_disp_header1,SYS1,path1,False)[0]
    df_final2=reg.csv(reg_inter2,iqa_inter_header2,reg_intra2,iqa_intra_header2,reg_disp2,iqa_disp_header2,SYS2,path2,False)[0]
    df_final1_fit=reg.csv(reg_inter1_fit,iqa_inter_header1,reg_intra1_fit,iqa_intra_header1,reg_disp1_fit,iqa_disp_header1,SYS_f1,path1,False)[0]
    df_final2_fit=reg.csv(reg_inter2_fit,iqa_inter_header2,reg_intra2_fit,iqa_intra_header2,reg_disp1_fit,iqa_disp_header2,SYS_f2,path2,False)[0]
    df_final=reg.csv(reg_inter,iqa_inter_header1,reg_intra,iqa_intra_header1,reg_disp,iqa_disp_header1,SYS,path,False)[0]

    dataframe_list1=reg.csv(reg_inter1,iqa_inter_header1,reg_intra1,iqa_intra_header1,reg_disp1,iqa_disp_header1,SYS1,path1,False)[1]
    dataframe_list2=reg.csv(reg_inter2,iqa_inter_header2,reg_intra2,iqa_intra_header2,reg_disp2,iqa_disp_header2,SYS2,path2,False)[1]
    dataframe_list1_fit=reg.csv(reg_inter1_fit,iqa_inter_header1,reg_intra1_fit,iqa_intra_header1,reg_disp1_fit,iqa_disp_header1,SYS_f1,path1,False)[1]
    dataframe_list2_fit=reg.csv(reg_inter2_fit,iqa_inter_header1,reg_intra2_fit,iqa_intra_header1,reg_disp1_fit,iqa_disp_header2,SYS_f2,path2,False)[1]
    dataframe_list=reg.csv(reg_inter,iqa_inter_header1,reg_intra,iqa_intra_header1,reg_disp,iqa_disp_header1,SYS,path,False)[1]
    
    reg.excel(iqa_inter1,iqa_inter_header1,iqa_intra1,iqa_intra_header1,iqa_disp1,iqa_disp_header1,SYS1,path1,False)
    reg.excel(iqa_inter2,iqa_inter_header2,iqa_intra2,iqa_intra_header2,iqa_disp2,iqa_disp_header2,SYS2,path2,False)
    reg.excel(iqa_inter1_fit,iqa_inter_header1,iqa_intra1_fit,iqa_intra_header1,iqa_disp1_fit,iqa_disp_header1,SYS_f1,path1,False)
    reg.excel(iqa_inter2_fit,iqa_inter_header2,iqa_intra2_fit,iqa_intra_header2,iqa_disp2_fit,iqa_disp_header2,SYS_f2,path2,False)
    reg.excel(iqa_inter,iqa_inter_header1,iqa_intra,iqa_intra_header1,iqa_disp,iqa_disp_header1,SYS,path,False)

    #rv.plot_violin([dataframe_list[i]['R'] for i in range(len(reg_inter[0]))], save=SAVE_FIG, file_name = '../' +'violin.png' )
    rv.plot_violin([dataframe_list1[i]['R2'] for i in range(len(reg_inter1[0]))], save=SAVE_FIG, file_name = path1 +'/'+SYS1+'_violin.png' )
    rv.plot_violin([dataframe_list2[i]['R2'] for i in range(len(reg_inter2[0]))], save=SAVE_FIG, file_name = path2 +'/'+SYS2+'_violin.png' )
    rv.plot_violin([dataframe_list1_fit[i]['R2'] for i in range(len(reg_inter1_fit[0]))], save=SAVE_FIG, file_name = path1 +'/'+SYS_f1+'_violin.png' )
    rv.plot_violin([dataframe_list2_fit[i]['R2'] for i in range(len(reg_inter2_fit[0]))], save=SAVE_FIG, file_name = path2 +'/'+SYS_f2+'_violin.png' )
    rv.plot_violin([dataframe_list[i]['R2'] for i in range(len(reg_inter[0]))], save=SAVE_FIG, file_name = path1 +'/'+SYS+'_violin.png' )
    
#if AUTO == True:
#    critical_points = reg.find_critical(total_energy_wfx, cc_a, min_points=POINTS, use_inflex=INFLEX)
#elif AUTO == False:
#    critical_points = tp

rv.plot_segment(cc_a1, 2625.50*(total_energy_wfx1-total_energy_wfx1[0]), crit_list1, annotate=ANNOTATE, label=LABELS,
                y_label = r'Relative Energy [$kJ.mol^{-1}$]', x_label=r'Reaction coordinate (a.m.u)', title = SYS1, save = SAVE_FIG, file_name=path1+'/'+SYS1+'REG_analysis.png',force = FORCE)   
rv.plot_segment(cc_a2, 2625.50*(total_energy_wfx2-total_energy_wfx2[0]), crit_list2, annotate=ANNOTATE, label=LABELS,
                y_label = r'Relative Energy [$kJ.mol^{-1}$]', x_label=r'Reaction coordinate (a.m.u)', title = SYS2, save = SAVE_FIG, file_name=path2+'/'+SYS2+'REG_analysis.png',force = FORCE)   
rv.plot_segment(cc_a_sampled, 2625.50*(total_energy_wfx1_fit-total_energy_wfx1_fit[0]), crit_list_fit, annotate=ANNOTATE, label=LABELS,
                y_label = r'Relative Energy [$kJ.mol^{-1}$]', x_label=r'Reaction coordinate (a.m.u)', title = SYS_f1, save = SAVE_FIG, file_name=path1+'/'+SYS_f1+'REG_analysis.png',force = FORCE)   
rv.plot_segment(cc_a_sampled, 2625.50*(total_energy_wfx2_fit-total_energy_wfx2_fit[0]), crit_list_fit, annotate=ANNOTATE, label=LABELS,
                y_label = r'Relative Energy [$kJ.mol^{-1}$]', x_label=r'Reaction coordinate (a.m.u)', title = SYS_f2, save = SAVE_FIG, file_name=path2+'/'+SYS_f2+'REG_analysis.png',force = FORCE)   
rv.plot_segment(cc_a_sampled, 2625.50*(total_energy_wfx-total_energy_wfx[0]), crit_list_fit, annotate=ANNOTATE, label=LABELS,
                y_label = r'Relative Energy [$kJ.mol^{-1}$]', x_label=r'Reaction coordinate (a.m.u)', title = SYS, save = SAVE_FIG, file_name=path+'/'+SYS+'REG_analysis.png',force = FORCE)   
#AC:Changed the energy plot by substracting the energy of the first point


#GROUP INTRA ATOMIC TERMS:
iqa_intra_g = np.array(reg.group_intra(iqa_intra, iqa_intra_header1, groups)[0])
iqa_intra_header_g = np.array(reg.group_intra(iqa_intra, iqa_intra_header1, groups)[1])

#GROUP INTER ATOMIC TERMS:
iqa_inter_g = np.array(reg.group_inter(iqa_inter, iqa_inter_header1, groups)[0])
iqa_inter_header_g = np.array(reg.group_inter(iqa_inter, iqa_inter_header1, groups)[1])

#INTRA ATOMIC CONTRIBUTION
reg_intra_g = reg.reg(total_energy_wfx, cc_a_sampled, iqa_intra_g, np=POINTS, critical=AUTO_fit, inflex=INFLEX, critical_index=crit_list_fit)
#INTER ATOMIC CONTRIBUTION
reg_inter_g = reg.reg(total_energy_wfx, cc_a_sampled, iqa_inter_g, np=POINTS, critical=AUTO_fit, inflex=INFLEX, critical_index=crit_list_fit)
#INTER ATOMIC CONTRIBUTION
reg_disp_g = reg.reg(total_energy_wfx, cc_a_sampled, iqa_disp, np=POINTS, critical=AUTO_fit, inflex=INFLEX, critical_index=crit_list_fit)

##############################################################################
#                                                                             #
#                           WRITE CSV FILES (GROUPED)                         #
#                                                                             #
###############################################################################
#dataframe_list=[]

if WRITE == True:
    """
    df_final = pd.DataFrame()   
    for i in range(len(reg_inter[0])):                                 
        temp = [reg_inter_g[0][i], reg_inter_g[1][i]]                  #inter
        df_inter = pd.DataFrame(temp).transpose()                                 #inter
        df_inter.columns = ["REG", "R"]                                           #inter
        df_inter.index = iqa_inter_header_g                                 #inter
        df_inter.to_csv('../' + "g_Inter_seg_" +str(i+1) +".csv", sep=',')   #inter 
        df_inter = df_inter.rename_axis('TERM').reset_index()                    #inter
        temp = [reg_intra_g[0][i], reg_intra_g[1][i]]                  #intra
        df_intra = pd.DataFrame(temp).transpose()                                 #intra
        df_intra.columns = ["REG", "R"]                                           #intra
        df_intra.index = iqa_intra_header_g                                 #intra
        df_intra.to_csv('../' + "g_Intra_seg_" +str(i+1) + ".csv", sep=',')  #intra
        df_intra = df_intra.rename_axis('TERM').reset_index()                    #intra
    
        df_temp = pd.concat([df_intra, df_inter]).sort_values('REG').reset_index()
    
        df_temp['R'] = df_temp['R']**2
        dataframe_list.append(df_temp)
        df_final=pd.concat([df_final, pd.concat([df_temp[-5:], df_temp[:5]], axis=0).sort_values('REG', ascending=False)], axis=1)
        
    df_final.to_csv('../' + 'g_IQA_segment_analysis.csv', sep=',')
    """

    df_final=reg.csv(reg_inter_g,iqa_inter_header_g,reg_intra_g,iqa_intra_header_g,reg_disp_g,iqa_disp_header1,SYS,path,True)[0]
    dataframe_list_g=reg.csv(reg_inter_g,iqa_inter_header_g,reg_intra_g,iqa_intra_header_g,reg_disp_g,iqa_disp_header1,SYS,path,True)[1]
    reg.excel(iqa_inter_g,iqa_inter_header_g,iqa_intra_g,iqa_intra_header_g,iqa_disp,iqa_disp_header1,SYS,path,True)
    rv.plot_violin([dataframe_list_g[i]['R2'] for i in range(len(reg_inter_g[0]))], save=SAVE_FIG, file_name = path1 +'/'+SYS+'_violin_g.png' )
   
rv.plot_segment(cc_a_sampled, 2625.50*(total_energy_wfx-total_energy_wfx[0]), crit_list_fit, annotate=ANNOTATE, label=LABELS,
                y_label = r'Relative Energy [$kJ.mol^{-1}$]', x_label=r'Reaction coordinate (a.m.u)', title = SYS, 
                save = SAVE_FIG, file_name=path+'/'+SYS+'REG_analysis.png',force = FORCE)   

if DETAILED_ANALYSIS ==True:
    
    for i in range(len(reg_inter1[0])): #AC change: R to R2.
        rv.generate_data_vis(dataframe_list1[i], [dataframe_list1[i]['R2'] for i in range(len(reg_inter1[0]))], 
                                            n_terms, save=SAVE_FIG, file_name=path1+'/'+SYS1+'_detailed_seg_'+ str(i+1) +'.png', 
                                            title = SYS1 +' seg. '+str(i+1), color = COLOR_TERMS)
    for i in range(len(reg_inter2[0])): #AC change: R to R2.
        rv.generate_data_vis(dataframe_list2[i], [dataframe_list2[i]['R2'] for i in range(len(reg_inter2[0]))], 
                                            n_terms, save=SAVE_FIG, file_name=path2+'/'+SYS2+'_detailed_seg_'+ str(i+1) +'.png', 
                                            title = SYS2 +' seg. '+str(i+1), color = COLOR_TERMS)
    for i in range(len(reg_inter1_fit[0])): #AC change: R to R2.
        rv.generate_data_vis(dataframe_list1_fit[i], [dataframe_list1_fit[i]['R2'] for i in range(len(reg_inter1_fit[0]))], 
                                            n_terms, save=SAVE_FIG, file_name=path1+'/'+SYS_f1+'_detailed_seg_'+ str(i+1) +'.png', 
                                            title = SYS_f1 +' seg. '+str(i+1), color = COLOR_TERMS)
    for i in range(len(reg_inter2_fit[0])): #AC change: R to R2.
        rv.generate_data_vis(dataframe_list2_fit[i], [dataframe_list2_fit[i]['R2'] for i in range(len(reg_inter2_fit[0]))], 
                                            n_terms, save=SAVE_FIG, file_name=path2+'/'+SYS_f2+'_detailed_seg_'+ str(i+1) +'.png', 
                                            title = SYS_f2 +' seg. '+str(i+1), color = COLOR_TERMS)
    for i in range(len(reg_inter[0])): #AC change: R to R2.
        rv.generate_data_vis(dataframe_list[i], [dataframe_list[i]['R2'] for i in range(len(reg_inter[0]))], 
                                            n_terms, save=SAVE_FIG, file_name=path+'/'+SYS+'_detailed_seg_'+ str(i+1) +'.png', 
                                            title = SYS +' seg. '+str(i+1), color = COLOR_TERMS)
    for i in range(len(reg_inter[0])): #AC change: R to R2.
        rv.generate_data_vis(dataframe_list_g[i], [dataframe_list_g[i]['R2'] for i in range(len(reg_inter_g[0]))], 
                                            n_terms, save=SAVE_FIG, file_name=path+'/'+SYS+'_detailed_seg_'+ str(i+1) +'_g.png', 
                                            title = SYS +' grouped seg. '+str(i+1), color = COLOR_TERMS)
                                            
    tot_en_wfx1=pd.DataFrame(total_energy_wfx1).transpose()
    tot_en_wfx2=pd.DataFrame(total_energy_wfx2).transpose()
    tot_en_wfx1_fit=pd.DataFrame(total_energy_wfx1_fit).transpose()
    tot_en_wfx2_fit=pd.DataFrame(total_energy_wfx2_fit).transpose() 
    tot_en_wfx=pd.DataFrame(total_energy_wfx).transpose()
    with pd.ExcelWriter(path + '/' + SYS + '_energy.xlsx') as writer: #writer object needed to have the lists written in different sheets
        #also need the data to be a dataframe to enable writing to XLSX file
        tot_en_wfx1.to_excel(writer, sheet_name='E_WFX_OUT')
        tot_en_wfx2.to_excel(writer, sheet_name='E_WFX_IN')
        tot_en_wfx1_fit.to_excel(writer, sheet_name='E_WFX_FIT_OUT')
        tot_en_wfx1_fit.to_excel(writer, sheet_name='E_WFX_FIT_OUT')
        tot_en_wfx.to_excel(writer, sheet_name='E_WFX_DIFF')
    cc_a1d=pd.DataFrame(cc_a1).transpose() 
    cc_a2d=pd.DataFrame(cc_a2).transpose() 
    cc_a1rd=pd.DataFrame(cc_a1_renorm).transpose() 
    cc_a2rd=pd.DataFrame(cc_a2_renorm).transpose() 
    cc_ad=pd.DataFrame(cc_a_sampled).transpose() 
    with pd.ExcelWriter(path + '/' + SYS + '_control_coordinates.xlsx') as writer: #writer object needed to have the lists written in different sheets
        #also need the data to be a dataframe to enable writing to XLSX file
        cc_a1d.to_excel(writer, sheet_name='cc_a1')
        cc_a2d.to_excel(writer, sheet_name='cc_a2')
        cc_a1rd.to_excel(writer, sheet_name='cc_a1_renorm')
        cc_a2rd.to_excel(writer, sheet_name='cc_a2_renorm')
        cc_ad.to_excel(writer, sheet_name='cc_a_sampled')
