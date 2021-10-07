# -*- coding: utf-8 -*-
"""
Created on Sun Apr  4 14:37:52 2021

@author: Aël
"""

# IMPORT MODULES
import reg 
import aimall_utils
import numpy as np
import sys
#import numpy.polynomial.polynomial as poly
import matplotlib.pyplot as plt
import pandas as pd
import reg_vis as rv

#from contextlib import redirect_stdout

"""
#original parameter input - you can copy/paste the content of a parameter file in there.
##############################    OPTIONS    ##################################
#SYSTEM NAME
SYS0 = 'BH2'

#DEFINE PATHS AND FILES FOR PATH 1: 
path1 = SYS0+"-out" #Path to the folder with .wfn 
F1=[i for i in range(0,155,5)] #the Forward SP/AIM file numbers ; don't forget to add 5 to the last one (e.g. F0-F70 -> 0,75)
R1=[i for i in range(5,70,5)] #the Reverse SP/AIM file numbers ; don't forget to add 5 to the last one

#DEFINE PATHS AND FILES FOR PATH 2:
path2 = SYS0+"-in" #Path to the folder with .wfn 
F2=[i for i in range(0,130,5)] #the Forward SP/AIM file numbers ; don't forget to add 5 to the last one
R2=[i for i in range(5,65,5)] #the Reverse SP/AIM file numbers ; don't forget to add 5 to the last one

#CALCULATION OPTIONS
inv_1=False #True if the file numbering is inverted for path 1 (Forward->Reactant and Reverse->Product)
inv_2=False #True if the file numbering is inverted for path 2 (Forward->Reactant and Reverse->Product)
groups=[['b10','h11','h12']]
POINTS = 4 #minimal number of points on each segment
INFLEX = True #if the inflexion points have to be considered as critical points or not

#OUTPUT OPTIONS
WRITE = True #if you want the script to write CSV and XLSX files with term values, REG values, energy values, control coordinate values
AUTO = True #Automatically divide the energy profile into segments, or not 
AUTO_fit = False #Automatically divide the energy difference profile into segments, or not 
FORCE = False #If you want the force, or not
COLOR_TERMS = False # If you want your term names to be colored in the segment analysis, or not
SIZE_TITLE = 18 #fontsize of plot titles
SIZE_LABEL = 16 #fontsize of plot axis and tick labels
SAVE_FIG = True #if you want to save plot figures to the disk
ANNOTATE = True #if you want to write the segment number on the complete energy profile plot or not
DETAILED_ANALYSIS = True #if you want to have the plots of significant terms for each segment or not
LABELS = True #if you want to have point labels on the complete energy profile plot or not
n_terms = 5 #number of terms to be displayed on each segment graph (=number of significant terms)
"""


def REGdiff(SYS0,path1,F1,R1,path2,F2,R2,inv_1,inv_2,groups,POINTS,INFLEX,WRITE,AUTO,AUTO_fit,FORCE,COLOR_TERMS,COLOR_SIGN,SIZE_TITLE,
            SIZE_LABEL,SAVE_FIG,ANNOTATE,DETAILED_ANALYSIS,LABELS,n_terms,lang="en"):
    path = SYS0+'-diff/'#test_C1+C2/'
    verbose = False #only for debugging purposes. Don't set to True unless you want to have a stream of blurb on your screen.
    if verbose == False:
        out = path+"out_noverbose.txt" 
    else:
        out = path+"out_verbose.txt"
        sys.stdout = open(out, "w")
       
    #with open(out, 'w') as f:
        #with redirect_stdout(f):
    print("====================== STARTING REGDIFF ======================")
    
    ###############################################################################
    #                                                                             #
    #                                 SETUP FILES                                 #
    #                                                                             #
    ###############################################################################

    print("====================== SETTING UP FILES ======================")
    
    #COMMON
    SYS1 = SYS0+'-out'
    SYS2 = SYS0+'-in'
    SYS_f1 = SYS0+'-fit-out'
    SYS_f2 = SYS0+'-fit-in'
    SYS = SYS0+'-diff (in-out)'
    
    model = 'RB3LYP'
    #path1 += '/test/'
    #path2 += '/test/'
    NOFIT = True # perform analysis of raw data
    FIT = False # perform analysis of fitted data
    DIFF = False # perform analysis of fitted, substracted data
    DOGROUP = False # perform analysis of grouped data (group atoms with "groups" variable)
    DOGROUP2 = False # perform analysis of grouped data (group electrostatic Vc and exchange Vx terms)
    
    #GET FILES AND CONTROL COORDINATES FOR PATH 1:
    res_get_list_wfx1=aimall_utils.get_list_wfx(path1,SYS0,"out",F1,R1,inv_1, verbose=verbose)
    cc_a1=res_get_list_wfx1[0]
    folders1=res_get_list_wfx1[1]
    folders_disp1=res_get_list_wfx1[2]
    wfx_files1=res_get_list_wfx1[3]
    
    res_get_list_wfx2=aimall_utils.get_list_wfx(path2,SYS0,"in",F2,R2,inv_2, verbose=verbose)
    cc_a2=res_get_list_wfx2[0]
    folders2=res_get_list_wfx2[1]
    folders_disp2=res_get_list_wfx2[2]
    wfx_files2=res_get_list_wfx2[3]
    
    #GET ATOM LIST FROM ANY .WFX FILE:
    atoms1 = aimall_utils.get_atom_list_wfx(wfx_files1[0], verbose=verbose)
    atoms2 = aimall_utils.get_atom_list_wfx(wfx_files2[0], verbose=verbose)
    
    #CHECK THE MODEL USED IN AIMALL
    aimall_utils.model_check(wfx_files1, model, folders1, atoms1, SYS1, verbose=verbose)
    aimall_utils.model_check(wfx_files2, model, folders2, atoms2, SYS2, verbose=verbose)
    
    ###############################################################################
    #                                                                             #
    #                           GET TERMS AND ENERGIES                            #
    #                                                                             #
    ###############################################################################
    
    print("================= READING TERMS AND ENERGIES =================")
    
    #DEFINE THE DESIRED TERMS:
    intra_prop = ['E_IQA_Intra(A)'] #intra atomic properties
    inter_prop = ['VC_IQA(A,B)', 'VX_IQA(A,B)'] #inter atomic properties
    
    #GET TOTAL ENERGY FROM THE .WFN FILES:
    total_energy_wfx1 = aimall_utils.get_aimall_wfx_energies(wfx_files1, verbose=verbose)
    total_energy_wfx2 = aimall_utils.get_aimall_wfx_energies(wfx_files2, verbose=verbose)

    #GET INTRA-ATOMIC TERMS:
    res_iqa_intra1 = aimall_utils.intra_property_from_int_file(folders1, intra_prop, atoms1, verbose=verbose)
    iqa_intra1 = res_iqa_intra1[0]
    iqa_intra_header1 = np.array(res_iqa_intra1[1]) #used for reference
    res_iqa_intra2 = aimall_utils.intra_property_from_int_file(folders2, intra_prop, atoms2, verbose=verbose)
    iqa_intra2 = res_iqa_intra2[0]
    iqa_intra_header2 = np.array(res_iqa_intra2[1]) #used for reference
        
    #change the term names to shorten them
    for i in range(len(iqa_intra_header1)):
            iqa_intra_header1[i] = iqa_intra_header1[i].replace("E_IQA_Intra(A)","E_Intra")
            iqa_intra_header2[i] = iqa_intra_header2[i].replace("E_IQA_Intra(A)","E_Intra")
    
    #GET INTER-ATOMIC TERMS:
    res_iqa_inter1 = aimall_utils.inter_property_from_int_file(folders1, inter_prop, atoms1, verbose=verbose)
    iqa_inter1 = res_iqa_inter1[0]
    iqa_inter_header1 = np.array(res_iqa_inter1[1]) #used for reference
    res_iqa_inter2 = aimall_utils.inter_property_from_int_file(folders2, inter_prop, atoms2, verbose=verbose)
    iqa_inter2 = res_iqa_inter2[0]
    iqa_inter_header2 = np.array(res_iqa_inter2[1]) #used for reference
    
    #change the term names to shorten them and change Vc to Vcl
    for i in range(len(iqa_inter_header1)):
            iqa_inter_header1[i] = iqa_inter_header1[i].replace("C_IQA(A,B)","cl")
            iqa_inter_header2[i] = iqa_inter_header2[i].replace("C_IQA(A,B)","cl")
            iqa_inter_header1[i] = iqa_inter_header1[i].replace("X_IQA(A,B)","X")
            iqa_inter_header2[i] = iqa_inter_header2[i].replace("X_IQA(A,B)","X")
    
    #GET TOTAL DISPERSION TERMS:
    res_iqa_disp1 = aimall_utils.disp_property(folders_disp1, verbose=verbose)
    iqa_disp1 = res_iqa_disp1[0]
    iqa_disp_header1 = np.array(res_iqa_disp1[1]) #used for reference
    res_iqa_disp2 = aimall_utils.disp_property(folders_disp2, verbose=verbose)
    iqa_disp2 = res_iqa_disp2[0]
    iqa_disp_header2 = np.array(res_iqa_disp2[1]) #used for reference
    
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
    
    ###############################################################################
    #                                                                             #
    #                               NORMALIZATION                                 #
    #                                                                             #
    ###############################################################################
    
    print("======================== NORMALIZATION =======================")
    
    #build list of critical points
    crit_list1 = reg.find_critical(total_energy_wfx1, cc_a1, use_inflex=True, verbose=verbose) 
    crit_list2 = reg.find_critical(total_energy_wfx2, cc_a2, use_inflex=True, verbose=verbose)
    
    #harmonise the lists in case the pathways don't have the same number of critical points
    crit_trim = reg.crit_trim(crit_list1, crit_list2, SYS1, SYS2)
    crit_list1 = crit_trim[0]
    crit_list2 = crit_trim[1]
    #crit_list1.pop(0)
    #crit_list2.pop(0)
    
    #NORMALIZE THE CONTROL COORDINATES           
    cc_a1_renorm = reg.normalize(cc_a1, crit_list1, verbose=verbose)
    cc_a2_renorm = reg.normalize(cc_a2, crit_list2, verbose=verbose)
    
    #some plots to show the normalisation results
    rv.plot_energy(total_energy_wfx1, total_energy_wfx2, cc_a1, cc_a2, crit_list1, crit_list2, SYS0, SIZE_TITLE, SIZE_LABEL, path, lang) #1 = OUT, 2 = IN
    rv.plot_renorm(total_energy_wfx1, total_energy_wfx2, cc_a1, cc_a2, cc_a1_renorm, cc_a2_renorm, crit_list1, SYS0, SIZE_TITLE, SIZE_LABEL,path,lang) #1 = OUT, 2 = IN
    
    #some more plots that I could probably get rid of
    plot0o = plt.figure(figsize=(10,7)) #ENERGY OUT ONLY
    plot1 = plot0o.add_subplot(111)
    plot1.plot(cc_a1,2625.5*(total_energy_wfx1-total_energy_wfx1[0]),'o--',label='$E_{mol}$', color='green')
    for i in crit_list1: #split the segments
        plot1.axvline(cc_a1[i], linestyle = 'dotted', color = 'blue')
    #plot1.set_title('$\Delta E$ ($kJ mol^{-1}$)',fontsize=18)
    plot1.set_xlabel(r'$\xi$ (a.m.u.)',fontsize=18)
    plot1.set_xlim(cc_a1[0],cc_a1[crit_list1[0]])
    plot1.set_ylabel('$\Delta E$ ($kJ.mol^{-1}$)',fontsize=18)
    plot1.set_ylim(-80,80)
    plot1.legend(fontsize=18,loc='best',labelcolor='linecolor')
    plot1.tick_params(axis='both',labelsize=18)
    plot0o.savefig('BH2-diff/'+SYS1+'_test_energy')
    plot1.plot(cc_a1,2625.5*(iqa_inter1[35]-iqa_inter1[35][0]),'o--',label='D1',color = 'blue')
    plot1.legend(fontsize=18,loc='best',labelcolor='linecolor')
    plot0o.savefig('BH2-diff/'+SYS1+'_test_D1')
    plot1.plot(cc_a1,2625.5*(iqa_inter1[18]-iqa_inter1[18][0]),'o--',label='D2',color = 'red')
    plot1.legend(fontsize=18,loc='best',labelcolor='linecolor')
    plot0o.savefig('BH2-diff/'+SYS1+'_test_D1+D2')
    
    plot0i = plt.figure(figsize=(10,7)) #ENERGY OUT ONLY
    plot2 = plot0i.add_subplot(111)
    plot2.plot(2625.5*(total_energy_wfx1[:crit_list1[0]+1]-total_energy_wfx1[0]),2625.5*(iqa_inter1[35][:crit_list1[0]+1]-iqa_inter1[35][0]),'o',label='D1',color = 'blue')
    plot2.plot(2625.5*(total_energy_wfx1[:crit_list1[0]+1]-total_energy_wfx1[0]),2625.5*(iqa_inter1[18][:crit_list1[0]+1]-iqa_inter1[18][0]),'o',label='D2',color = 'red')
    plot2.set_xlabel('$\Delta E_{mol}$ ($kJ.mol^{-1}$)',fontsize=18, color='green')
    plot2.set_ylabel('$D_{i}$ ($kJ.mol^{-1}$)',fontsize=18, color='purple')
    regc1=reg.regression((total_energy_wfx1[:crit_list1[0]+1]), iqa_inter1[18][:crit_list1[0]+1], "norm")[0]
    rege1=reg.regression((total_energy_wfx1[:crit_list1[0]+1]), iqa_inter1[18][:crit_list1[0]+1], "norm")[1]
    regr1=reg.regression((total_energy_wfx1[:crit_list1[0]+1]), iqa_inter1[18][:crit_list1[0]+1], "norm")[2]
    plot2.annotate("R² = "+str(regr1**2)[:6],xy=(52,-40),fontsize=18, color = 'red')
    plot2.annotate("$m_{REG} = $"+str(regc1)[:6],xy=(52,-30),fontsize=18, color = 'red')
    regc2=reg.regression((total_energy_wfx1[:crit_list1[0]+1]), iqa_inter1[35][:crit_list1[0]+1], "norm")[0]
    rege2=reg.regression((total_energy_wfx1[:crit_list1[0]+1]), iqa_inter1[35][:crit_list1[0]+1], "norm")[1]
    regr2=reg.regression((total_energy_wfx1[:crit_list1[0]+1]), iqa_inter1[35][:crit_list1[0]+1], "norm")[2]
    x1=np.linspace(0,2625.5*(total_energy_wfx1[crit_list1[0]+1]-total_energy_wfx1[0]),100)
    y1=regc1*x1+rege1
    plot2.plot(x1,y1,'-r')
    y2=regc2*x1+rege2
    plot2.plot(x1,y2,'-b')
    plot2.annotate("R² = "+str(regr2**2)[:6],xy=(52,10),fontsize=18, color = 'blue')
    plot2.annotate("$m_{REG} = $"+str(regc2)[:6],xy=(52,20),fontsize=18, color = 'blue')
    plot2.legend(fontsize=18,loc='best',labelcolor='linecolor')
    
    plot2.tick_params(axis='both',labelsize=18)
    plot0i.savefig('BH2-diff/'+SYS1+'_test2')
    """
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
    #GENERATE THE FULL CRITICAL POINT LIST (+EXTREMA AND SEGMENT LENGTHS)
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
    
    crit_list_fit_0N=[0]
    crit_list_fit_0N += crit_list2
    crit_list_fit_0N.append(len(cc_a2_renorm)-1)
    
    maxlong=max(maxlong1,maxlong2)
    
    #generate the critical point list of the normalized reaction coordinate; all segments have the same lengths which is "maxlong".
    crit_list_fit=[(maxlong-1)*i for i in range(1,len(crit_list1)+1)]
    print("Critical point list:",crit_list1)
    print("iqa_intra1 length",len(iqa_intra1))
    
    ###############################################################################
    #                                                                             #
    #                              GET TERM SIGNS                                 #
    #                                                                             #
    ###############################################################################
    
    print("======================= GET TERM SIGNS =======================")
    
    #Read term signs
    sign_seg_intra1 = reg.sign_segm(iqa_intra1, crit_list1, verbose=verbose)
    sign_seg_inter1 = reg.sign_segm(iqa_inter1, crit_list1, verbose=verbose)
    sign_seg_disp1 = reg.sign_segm(iqa_disp1, crit_list1, verbose=verbose)
    sign_seg_intra2 = reg.sign_segm(iqa_intra2, crit_list2, verbose=verbose)
    sign_seg_inter2 = reg.sign_segm(iqa_inter2, crit_list2, verbose=verbose)
    sign_seg_disp2 = reg.sign_segm(iqa_disp2, crit_list2, verbose=verbose)
    sign_seg1 = [sign_seg_intra1, sign_seg_inter1, sign_seg_disp1]
    sign_seg2 = [sign_seg_intra2, sign_seg_inter2, sign_seg_disp2]
    
    #some manipulations to get term signs for substracted terms (eg. to know if a term is positive on one pathway and negative on the other)
    sign_seg_intra1a = np.array(sign_seg_intra1)
    sign_seg_intra2a = np.array(sign_seg_intra2)
    sign_seg_intra = sign_seg_intra1a + sign_seg_intra2a
    sign_seg_intra = sign_seg_intra.tolist()
    sign_seg_inter1a = np.array(sign_seg_inter1)
    sign_seg_inter2a = np.array(sign_seg_inter2)
    sign_seg_inter = sign_seg_inter1a + sign_seg_inter2a
    sign_seg_inter = sign_seg_inter.tolist()
    sign_seg_disp1a = np.array(sign_seg_disp1)
    sign_seg_disp2a = np.array(sign_seg_disp2)
    sign_seg_disp = sign_seg_disp1a + sign_seg_disp2a
    sign_seg_disp = sign_seg_disp.tolist()
    sign_seg = [sign_seg_intra, sign_seg_inter, sign_seg_disp]
    #print(sign_seg)
    
    ###############################################################################
    #                                                                             #
    #                                   FITTING                                   #
    #                                                                             #
    ###############################################################################
    
    print("=========================== FITTING ==========================")
    
    #FITTING
    e1_cc = reg.fit(total_energy_wfx1, crit_list1_0N, maxlong, cc_a1_renorm, verbose=verbose)
    total_energy_wfx1_fit = e1_cc[0]
    cc_a_sampled = e1_cc[1]
    iqa_intra1_fit = reg.fit(iqa_intra1, crit_list1_0N, maxlong, cc_a1_renorm, verbose=verbose)[0]
    iqa_inter1_fit = reg.fit(iqa_inter1, crit_list1_0N, maxlong, cc_a1_renorm, verbose=verbose)[0]
    iqa_disp1_fit = reg.fit(iqa_disp1, crit_list1_0N, maxlong, cc_a1_renorm, verbose=verbose)[0]
    #print(iqa_intra1_fit)
    
    total_energy_wfx2_fit = reg.fit(total_energy_wfx2, crit_list2_0N, maxlong, cc_a2_renorm, verbose=verbose)[0]
    iqa_intra2_fit = reg.fit(iqa_intra2, crit_list2_0N, maxlong, cc_a2_renorm, verbose=verbose)[0]
    iqa_inter2_fit = reg.fit(iqa_inter2, crit_list2_0N, maxlong, cc_a2_renorm, verbose=verbose)[0]
    iqa_disp2_fit = reg.fit(iqa_disp2, crit_list2_0N, maxlong, cc_a2_renorm, verbose=verbose)[0]
    
    rv.plot_fit(total_energy_wfx1,total_energy_wfx2,total_energy_wfx1_fit,total_energy_wfx2_fit,cc_a1_renorm,cc_a2_renorm,
                cc_a_sampled,crit_list_fit,SYS0,SIZE_TITLE,SIZE_LABEL,path,lang)
    
    ###############################################################################
    #                                                                             #
    #                                SUBSTRACTION                                 #
    #                                                                             #
    ###############################################################################
    
    print("========================= SUBSTRACTION ========================")
    
    #SUBSTRACT
    total_energy_wfx = total_energy_wfx2_fit - total_energy_wfx1_fit
    iqa_inter = iqa_inter2_fit - iqa_inter1_fit
    iqa_intra = iqa_intra2_fit - iqa_intra1_fit
    iqa_disp = iqa_disp2_fit - iqa_disp1_fit
    #print(iqa_intra[0])  
    
    rv.plot_substract(total_energy_wfx1_fit,total_energy_wfx2_fit,total_energy_wfx,cc_a_sampled,crit_list_fit,
                      SYS0,SIZE_TITLE,SIZE_LABEL,path,lang)
    
    rv.plot_substract_term(iqa_intra1_fit,iqa_intra1,iqa_intra2_fit,iqa_intra2,iqa_intra,iqa_intra_header1,0,
                           cc_a1_renorm,cc_a2_renorm,cc_a_sampled,crit_list_fit,SYS0,SIZE_TITLE,SIZE_LABEL,path,lang)
    rv.plot_substract_term(iqa_inter1_fit,iqa_inter1,iqa_inter2_fit,iqa_inter2,iqa_inter,iqa_inter_header1,0,
                           cc_a1_renorm,cc_a2_renorm,cc_a_sampled,crit_list_fit,SYS0,SIZE_TITLE,SIZE_LABEL,path,lang) 
    
    plot_ = plt.figure(figsize=(10,7))
    plot2 = plot_.add_subplot(111)

    if lang=="fr":
        plot2.plot(cc_a1_renorm,2625.5*(total_energy_wfx1-total_energy_wfx1[0]),color='deepskyblue',linestyle='--',marker='o',label='Énergie (out)') #markeredgecolor='darkblue'
        plot2.plot(cc_a2_renorm,2625.5*(total_energy_wfx2-total_energy_wfx2[0]),color='orange',linestyle='--',marker='o',label='Énergie (in)') #,markeredgecolor='red'
    
    plot2.tick_params(axis='both',labelsize=15)
    plot2.set_xticks([i for i in range(0,5)])
    plot2.grid(b=True,axis='x')
    #plot_.savefig(path+'/'+SYS0+'_energy_fitting')
    plot2.plot(cc_a_sampled,2625.5*(total_energy_wfx1_fit-total_energy_wfx1_fit[0]),
               color='blue',linestyle='--',marker='.')
    plot2.plot(cc_a_sampled,2625.5*(total_energy_wfx2_fit-total_energy_wfx2_fit[0]),
               color='red',linestyle='--',marker='.')
    if lang == "fr":
        plot2.plot(cc_a_sampled[crit_list_fit[-1]:],2625.5*(total_energy_wfx1_fit[crit_list_fit[-1]:]-total_energy_wfx1_fit[0]),
                   color='blue',linestyle='-',label='Énergie rééch. (out)')
        plot2.plot(cc_a_sampled[crit_list_fit[-1]:],2625.5*(total_energy_wfx2_fit[crit_list_fit[-1]:]-total_energy_wfx2_fit[0]),
                   color='red',linestyle='-',label='Énergie rééch. (in)')
        plot2.set_title('Rééchantillonnage de l\'énergie ('+SYS0+')',fontsize=SIZE_TITLE)
        plot2.set_xlabel(r'$\tilde \xi$ (a.m.u.)',fontsize=SIZE_LABEL) #Coordonnée de réaction normalisée
        plot2.set_ylabel('$\Delta E_{SCF}$ ($kJ.mol^{-1}$)',fontsize=SIZE_LABEL)
        plot2.set_ylim(-100,175)
        plot2.plot(cc_a_sampled,2625.5*(total_energy_wfx-total_energy_wfx[0]),color='limegreen',linestyle='-',marker='.',label='Énergie (diff. in-out)')
    else:
        plot2.plot(cc_a_sampled,2625.5*(total_energy_wfx-total_energy_wfx[0]),color='limegreen',linestyle='-',marker='.',label='Energy (in-out diff.)')
    plot2.legend(fontsize=SIZE_LABEL-2)
    plot_.savefig(path+'/'+SYS0+'_energy_fitting_Nref_'+lang)
    
    """
    #more plots
    for i in range(len(iqa_intra1_fit)):
        rv.plot_substract_term(iqa_intra1_fit,iqa_intra1,iqa_intra2_fit,iqa_intra2,iqa_intra,iqa_intra_header1,i,
                           cc_a1_renorm,cc_a2_renorm,cc_a_sampled,crit_list_fit,SYS0,SIZE_TITLE,SIZE_LABEL,path,lang)
    for i in range(len(iqa_inter1_fit)):
        rv.plot_substract_term(iqa_inter1_fit,iqa_inter1,iqa_inter2_fit,iqa_inter2,iqa_inter,iqa_inter_header1,i,
                           cc_a1_renorm,cc_a2_renorm,cc_a_sampled,crit_list_fit,SYS0,SIZE_TITLE,SIZE_LABEL,path,lang) 
    """
     
    #CALCULATE THE IQA TOTAL ENERGY 
    total_energy_iqa1 = sum(iqa_inter1) + sum(iqa_intra1) + sum(iqa_disp1) # used to calculate the integration error.
    total_energy_iqa2 = sum(iqa_inter2) + sum(iqa_intra2) + sum(iqa_disp2) # used to calculate the integration error.
    total_energy_iqa1_fit = reg.fit(total_energy_iqa1, crit_list1_0N, maxlong, cc_a1_renorm, verbose=verbose)[0]
    total_energy_iqa2_fit = reg.fit(total_energy_iqa2, crit_list2_0N, maxlong, cc_a2_renorm, verbose=verbose)[0]
    total_energy_iqa1_fit1 = sum(iqa_inter1_fit) + sum(iqa_intra1_fit) + sum(iqa_disp1_fit) # used to calculate the integration error.
    total_energy_iqa2_fit1 = sum(iqa_inter2_fit) + sum(iqa_intra2_fit) + sum(iqa_disp2_fit) # used to calculate the integration error.
    
    diff_total_energy_iqa = sum(iqa_inter) + sum(iqa_intra) + sum(iqa_disp) # used to calculate the integration error.
    #print(total_energy_iqa1,total_energy_iqa2,diff_total_energy_iqa)
    
    """
    start
    This is a section in progress to automate a "resegmentation" when the difference of the energy profiles isn't monotonous, which leads
    to awful results. Well, we'd like to try to change that to see if the results can be improved.
    """
    
    """
    Some French-
    pour les cas problématiques il faudrait faire une détection automatique avec la dérivée, un message "segment non monotone" 
    et une figure qui montre à l'utilisateur où sont les minima/maxima détectés sur le segment.
    Pour le moment on va faire à la mano
    """
    CRIT_DIFF = False
    if CRIT_DIFF == True:
        ###############################################################################
        #                                                                             #
        #                          DIFFERENCE CRITICAL POINTS                         #
        #                                                                             #
        ###############################################################################
        
        print("=================== DIFFERENCE CRITICAL POINTS ================")
        
        crit_list_diff = reg.find_critical(total_energy_wfx, cc_a_sampled, use_inflex=True, verbose=verbose)
        print("crit_list_diff",crit_list_diff)
        crit_list_diff2 = reg.find_critical(total_energy_wfx, cc_a_sampled, use_inflex=False, verbose=verbose)
        print("crit_list_diff2",crit_list_diff2)
        
        if SYS0 == 'BCl2':
            crit_list_diff2.pop(2)
            #crit_list_diff2.pop(2)
            #crit_list_diff2.pop(2)
        
        if SYS0 == 'BH3-':
            crit_list_diff2.pop(0)
            crit_list_diff2.pop(0)
        #    crit_list_diff2.pop(-2)
        #    crit_list_diff2.pop(-2)
        
        if SYS0 == 'BMe2':
            crit_list_diff2.insert(1,49)
            #crit_list_diff2.pop(3)
            #crit_list_diff2.pop(1)
        #    crit_list_diff2.pop(-2)
        
        if SYS0 == 'CH3':
            crit_list_diff2.pop(0)
            crit_list_diff2.pop(0)
            crit_list_diff2.pop(2)
            crit_list_diff2.pop(2)
            crit_list_diff2.pop(3)
            crit_list_diff2.pop(-1)
        
        if SYS0 == 'CHO':
            crit_list_diff2.pop(3)
            crit_list_diff2.pop(0)
        
        #if SYS0 == 'Cl':
        #    crit_list_diff2.pop(-1)
        
        if SYS0 == 'CN':
            for i in range(0,5):
                crit_list_diff2.pop(0)    
            crit_list_diff2.pop(-2)
        
        if SYS0 == 'F':
            crit_list_diff2.pop(0)
        
        if SYS0 == 'NO2':
            crit_list_diff2.pop(6)
            crit_list_diff2.pop(5)
            crit_list_diff2.pop(4)
            crit_list_diff2.pop(3)
            crit_list_diff2.pop(1)
        
        if SYS0 == 'NH2':
            crit_list_diff2.pop(0)    
        
        if SYS0 == 'NH3+':
            crit_list_diff2.pop(-1)
        
        if SYS0 == 'PH2':
            crit_list_diff2.pop(5)
            crit_list_diff2.pop(5)
            crit_list_diff2.pop(5)
            crit_list_diff2.pop(2)
            crit_list_diff2.pop(0) 

        if SYS0 == 'SiF3':
            crit_list_diff2.pop(-1)
            crit_list_diff2.pop(-1)
            #crit_list_diff2.pop(-2)
            crit_list_diff2.pop(-1)
            crit_list_diff2.pop(0)
            crit_list_diff2.pop(0)
            crit_list_diff2.pop(0)
             
        print("crit_list_diff2",crit_list_diff2)
        #cc_a_renorm = reg.normalize(cc_a_sampled,crit_list_diff2, verbose=verbose)
    
        crd = reg.norm_diff(cc_a1, cc_a2, cc_a_sampled, crit_list_diff2, crit_list1_0N, crit_list2_0N, verbose=verbose) #WARNING : use crit_list_diff2
        crd1=crd[0]
        crd2=crd[1]
        print("crd0",crd[0])
        print("crd1",crd[1])
        
        crit_list1_0N = [0]
        crit_list1_0N += crd1
        crit_list1_0N.append(len(cc_a1)-1)
        seg_long1 = [0]
        for i in range(len(crit_list1_0N)-1):
            seg_long1.append(crit_list1_0N[i+1]-crit_list1_0N[i])
        seg_long1.sort()
        maxlong1 = seg_long1[-1]
    
        crit_list2_0N = [0]
        crit_list2_0N += crd2
        crit_list2_0N.append(len(cc_a2)-1)
        seg_long2 = [0]
        for i in range(len(crit_list2_0N)-1):
            seg_long2.append(crit_list2_0N[i+1]-crit_list2_0N[i])
        seg_long2.sort()
        maxlong2 = seg_long2[-1]
        
        cc_a1_renorm = reg.normalize(cc_a1,crd[0], verbose=verbose)
        cc_a2_renorm = reg.normalize(cc_a2,crd[1], verbose=verbose)
        
        ###############################################################################
        #                                                                             #
        #                                   FITTING (2)                               #
        #                                                                             #
        ###############################################################################
        
        print("========================== REFITTING =========================")
        
        #FITTING
        total_energy_wfx1_fit = reg.fit(total_energy_wfx1, crit_list1_0N, maxlong, cc_a1_renorm, verbose=verbose)[0]
        cc_a_sampled = reg.fit(total_energy_wfx1, crit_list1_0N, maxlong, cc_a1_renorm, verbose=verbose)[1]
        iqa_intra1_fit = reg.fit(iqa_intra1, crit_list1_0N, maxlong, cc_a1_renorm, verbose=verbose)[0]
        iqa_inter1_fit = reg.fit(iqa_inter1, crit_list1_0N, maxlong, cc_a1_renorm, verbose=verbose)[0]
        iqa_disp1_fit = reg.fit(iqa_disp1, crit_list1_0N, maxlong, cc_a1_renorm, verbose=verbose)[0]
        #print(iqa_intra1_fit)
        total_energy_wfx2_fit = reg.fit(total_energy_wfx2, crit_list2_0N, maxlong, cc_a2_renorm, verbose=verbose)[0]
        iqa_intra2_fit = reg.fit(iqa_intra2, crit_list2_0N, maxlong, cc_a2_renorm, verbose=verbose)[0]
        iqa_inter2_fit = reg.fit(iqa_inter2, crit_list2_0N, maxlong, cc_a2_renorm, verbose=verbose)[0]
        iqa_disp2_fit = reg.fit(iqa_disp2, crit_list2_0N, maxlong, cc_a2_renorm, verbose=verbose)[0]
        
        rv.plot_fit(total_energy_wfx1,total_energy_wfx2,total_energy_wfx1_fit,total_energy_wfx2_fit,cc_a1_renorm,cc_a2_renorm,
                    cc_a_sampled,crit_list_fit,SYS0,SIZE_TITLE,SIZE_LABEL,path,lang)
        
        ###############################################################################
        #                                                                             #
        #                                SUBSTRACTION (2)                             #
        #                                                                             #
        ###################################################################### #########
        
        print("======================= RESUBSTRACTING ======================")
    
        #SUBSTRACT
        total_energy_wfx = total_energy_wfx2_fit - total_energy_wfx1_fit
        iqa_inter = iqa_inter2_fit - iqa_inter1_fit
        iqa_intra = iqa_intra2_fit - iqa_intra1_fit
        iqa_disp = iqa_disp2_fit - iqa_disp1_fit
        #print(iqa_intra[0])
        
        rv.plot_substract(total_energy_wfx1_fit, total_energy_wfx2_fit, total_energy_wfx, cc_a_sampled, crit_list_fit,
                          SYS0, SIZE_TITLE, SIZE_LABEL, path, lang)
        
        rv.plot_substract_term(iqa_intra1_fit, iqa_intra1, iqa_intra2_fit, iqa_intra2, iqa_intra, iqa_intra_header1, 0,
                               cc_a1_renorm, cc_a2_renorm, cc_a_sampled, crit_list_fit, SYS0, SIZE_TITLE, SIZE_LABEL, path, lang)
        rv.plot_substract_term(iqa_inter1_fit, iqa_inter1, iqa_inter2_fit, iqa_inter2, iqa_inter, iqa_inter_header1, 0,
                               cc_a1_renorm, cc_a2_renorm, cc_a_sampled, crit_list_fit, SYS0, SIZE_TITLE, SIZE_LABEL, path, lang)
        
        """
        total_energy_wfx1 = 2625.5*total_energy_wfx1
        total_energy_wfx2 = 2625.5*total_energy_wfx2
        total_energy_wfx1_fit = 2625.5*total_energy_wfx1_fit
        total_energy_wfx2_fit = 2625.5*total_energy_wfx2_fit
        total_energy_wfx = 2625.5*total_energy_wfx
        """
        ###############################################################################
        #                                                                             #
        #                              GET TERM SIGNS                                 #
        #                                                                             #
        ###############################################################################
        
        print("======================= GET TERM SIGNS =======================")
    
        #Read term signs
        sign_seg_intra1 = reg.sign_segm(iqa_intra1, crit_list1, verbose=verbose)
        sign_seg_inter1 = reg.sign_segm(iqa_inter1, crit_list1, verbose=verbose)
        sign_seg_disp1 = reg.sign_segm(iqa_disp1, crit_list1, verbose=verbose)
        sign_seg_intra2 = reg.sign_segm(iqa_intra2, crit_list2, verbose=verbose)
        sign_seg_inter2 = reg.sign_segm(iqa_inter2, crit_list2, verbose=verbose)
        sign_seg_disp2 = reg.sign_segm(iqa_disp2, crit_list2, verbose=verbose)
        sign_seg1 = [sign_seg_intra1, sign_seg_inter1, sign_seg_disp1]
        sign_seg2 = [sign_seg_intra2, sign_seg_inter2, sign_seg_disp2]
        
        #some manipulations to get term signs for substracted terms (eg. to know if a term is positive on one pathway and negative on the other)
        sign_seg_intra1a = np.array(sign_seg_intra1)
        sign_seg_intra2a = np.array(sign_seg_intra2)
        sign_seg_intra = sign_seg_intra1a + sign_seg_intra2a
        sign_seg_intra = sign_seg_intra.tolist()
        sign_seg_inter1a = np.array(sign_seg_inter1)
        sign_seg_inter2a = np.array(sign_seg_inter2)
        sign_seg_inter = sign_seg_inter1a + sign_seg_inter2a
        sign_seg_inter = sign_seg_inter.tolist()
        sign_seg_disp1a = np.array(sign_seg_disp1)
        sign_seg_disp2a = np.array(sign_seg_disp2)
        sign_seg_disp = sign_seg_disp1a + sign_seg_disp2a
        sign_seg_disp = sign_seg_disp.tolist()
        sign_seg = [sign_seg_intra, sign_seg_inter, sign_seg_disp]
        #print(sign_seg)
    
    ###############################################################################
    #                                                                             #
    #                               REG ANALYSIS                                  #
    #                                                                             #
    ###############################################################################
    
    print("======================== REG ANALYSIS =======================")
    
    #REG ANALYSIS
    
    #INTRA ATOMIC CONTRIBUTION
    if NOFIT == True:
        reg_intra1 = reg.reg(total_energy_wfx1, cc_a1, iqa_intra1, np=POINTS, critical=AUTO_fit, inflex=INFLEX, 
                             critical_index=crit_list1, verbose=verbose)
        reg_intra2 = reg.reg(total_energy_wfx2, cc_a2, iqa_intra2, np=POINTS, critical=AUTO_fit, inflex=INFLEX, 
                             critical_index=crit_list2, verbose=verbose)
    if FIT == True:
        reg_intra1_fit = reg.reg(total_energy_wfx1_fit, cc_a_sampled, iqa_intra1_fit, np=POINTS, critical=AUTO_fit, 
                                 inflex=INFLEX, critical_index=crit_list_fit, verbose=verbose)
        reg_intra2_fit = reg.reg(total_energy_wfx2_fit, cc_a_sampled, iqa_intra2_fit, np=POINTS, critical=AUTO_fit, 
                                 inflex=INFLEX, critical_index=crit_list_fit, verbose=verbose)
    if DIFF == True:
        reg_intra = reg.reg(total_energy_wfx, cc_a_sampled, iqa_intra, np=POINTS, critical=AUTO_fit, inflex=INFLEX, 
                            critical_index=crit_list_fit, verbose=verbose)
    
    #INTER ATOMIC CONTRIBUTION
    
    if NOFIT == True:
        reg_inter1 = reg.reg(total_energy_wfx1, cc_a1, iqa_inter1, np=POINTS, critical=AUTO_fit, inflex=INFLEX, 
                             critical_index=crit_list1, verbose=verbose)
        reg_inter2 = reg.reg(total_energy_wfx2, cc_a2, iqa_inter2, np=POINTS, critical=AUTO_fit, inflex=INFLEX, 
                             critical_index=crit_list2, verbose=verbose)
    if FIT == True:
        reg_inter1_fit = reg.reg(total_energy_wfx1_fit, cc_a_sampled, iqa_inter1_fit, np=POINTS, critical=AUTO_fit, 
                                 inflex=INFLEX, critical_index=crit_list_fit, verbose=verbose)
        reg_inter2_fit = reg.reg(total_energy_wfx2_fit, cc_a_sampled, iqa_inter2_fit, np=POINTS, critical=AUTO_fit, 
                                 inflex=INFLEX, critical_index=crit_list_fit, verbose=verbose)
    if DIFF == True:
        reg_inter = reg.reg(total_energy_wfx, cc_a_sampled, iqa_inter, np=POINTS, critical=AUTO_fit, inflex=INFLEX, 
                            critical_index=crit_list_fit, verbose=verbose)
    
    #GLOBAL DISPERSION CONTRIBUTION
    
    if NOFIT == True:
        reg_disp1 = reg.reg(total_energy_wfx1, cc_a1, iqa_disp1, np=POINTS, critical=AUTO_fit, inflex=INFLEX, 
                            critical_index=crit_list1, verbose=verbose)
        reg_disp2 = reg.reg(total_energy_wfx2, cc_a2, iqa_disp2, np=POINTS, critical=AUTO_fit, inflex=INFLEX, 
                            critical_index=crit_list2, verbose=verbose)
    if FIT == True:
        reg_disp1_fit = reg.reg(total_energy_wfx1_fit, cc_a_sampled, iqa_disp1_fit, np=POINTS, critical=AUTO_fit, 
                                inflex=INFLEX, critical_index=crit_list_fit, verbose=verbose)
        reg_disp2_fit = reg.reg(total_energy_wfx2_fit, cc_a_sampled, iqa_disp2_fit, np=POINTS, critical=AUTO_fit, 
                                inflex=INFLEX, critical_index=crit_list_fit, verbose=verbose)
    if DIFF == True:
        reg_disp = reg.reg(total_energy_wfx, cc_a_sampled, iqa_disp, np=POINTS, critical=AUTO_fit, inflex=INFLEX, 
                           critical_index=crit_list_fit, verbose=verbose)
    
    #GROUP INTER ATOMIC TERMS:
    #iqa_disp_g = np.array(reg.group_inter(iqa_disp, iqa_disp_header, groups, verbose=verbose)[0])
    #iqa_disp_header_g = np.array(reg.group_inter(iqa_disp, iqa_disp_header, groups, verbose=verbose)[1])
    
    ###############################################################################
    #                                                                             #
    #                             INTEGRATION ERROR                               #
    #                                                                             #
    ###############################################################################
    
    print("====================== INTEGRATION ERROR ====================")
    
    #CALCULATE INTEGRATION (RECOVERY) ERROR:
    rmse_integration = reg.integration_error(total_energy_wfx1, total_energy_iqa1)
    print('Integration error [kJ/mol](RMSE) - Pathway 1')
    print("%.2f" % rmse_integration[1])
    rmse_integration = reg.integration_error(total_energy_wfx2, total_energy_iqa2)
    print('Integration error [kJ/mol](RMSE) - Pathway 2')
    print("%.2f" % rmse_integration[1])
    rmse_integration = reg.integration_error(total_energy_wfx1_fit, total_energy_iqa1_fit)
    print('Integration error [kJ/mol](RMSE) - Pathway 1 (fitted)')
    print("%.2f" % rmse_integration[1])
    rmse_integration = reg.integration_error(total_energy_wfx2_fit, total_energy_iqa2_fit)
    print('Integration error [kJ/mol](RMSE) - Pathway 2 (fitted)')
    print("%.2f" % rmse_integration[1])
    rmse_integration = reg.integration_error(total_energy_wfx1_fit, total_energy_iqa1_fit1)
    print('Integration error [kJ/mol](RMSE) - Pathway 1 (fitted)')
    print("%.2f" % rmse_integration[1])
    rmse_integration = reg.integration_error(total_energy_wfx2_fit, total_energy_iqa2_fit1)
    print('Integration error [kJ/mol](RMSE) - Pathway 2 (fitted)')
    print("%.2f" % rmse_integration[1])
    
    ###############################################################################
    #                                                                             #
    #                     WRITE CSV FILES (WITHOUT GROUPS)                        #
    #                                                                             #
    ###############################################################################
    
    print("===================== WRITING CSV FILES ===================")
    
    if WRITE == True: 
        
        #CSV files
        if NOFIT == True:
            dataframe_list1=reg.csv(reg_inter1, iqa_inter_header1, reg_intra1, iqa_intra_header1, reg_disp1, iqa_disp_header1, SYS1,
                                    path1, False, sign_seg1, verbose=verbose)#[1]
            dataframe_list2=reg.csv(reg_inter2, iqa_inter_header2, reg_intra2, iqa_intra_header2, reg_disp2, iqa_disp_header2, SYS2,
                                    path2, False, sign_seg2, verbose=verbose)#[1]
        if FIT == True:
            dataframe_list1_fit=reg.csv(reg_inter1_fit, iqa_inter_header1, reg_intra1_fit, iqa_intra_header1, reg_disp1_fit,
                                        iqa_disp_header1, SYS_f1, path1, False, sign_seg1, verbose=verbose)
            dataframe_list2_fit=reg.csv(reg_inter2_fit, iqa_inter_header1, reg_intra2_fit, iqa_intra_header1, reg_disp2_fit,
                                        iqa_disp_header2, SYS_f2, path2, False, sign_seg2, verbose=verbose)
        if DIFF == True:
            dataframe_list=reg.csv(reg_inter, iqa_inter_header1, reg_intra, iqa_intra_header1, reg_disp, iqa_disp_header1, SYS, path,
                                   False, sign_seg, verbose=verbose)
        
        #Excel files
        if NOFIT == True:
            reg.excel(iqa_inter1, iqa_inter_header1, iqa_intra1, iqa_intra_header1, iqa_disp1, iqa_disp_header1, SYS1, path1, False)
            reg.excel(iqa_inter2, iqa_inter_header2, iqa_intra2, iqa_intra_header2, iqa_disp2, iqa_disp_header2, SYS2, path2, False)
        if FIT == True:
            reg.excel(iqa_inter1_fit, iqa_inter_header1, iqa_intra1_fit, iqa_intra_header1, iqa_disp1_fit, iqa_disp_header1, SYS_f1, path1, False)
            reg.excel(iqa_inter2_fit, iqa_inter_header2, iqa_intra2_fit, iqa_intra_header2, iqa_disp2_fit, iqa_disp_header2, SYS_f2, path2, False)
        if DIFF == True:
            reg.excel(iqa_inter, iqa_inter_header1, iqa_intra, iqa_intra_header1, iqa_disp, iqa_disp_header1, SYS, path, False)
        
        #print("DATA",dataframe_list1[0]['R2'])
        
        #GENERATE ENERGY AND CONTROL COORDINATE XLSX FILES
        if NOFIT == True:
            tot_en_wfx1=pd.DataFrame(total_energy_wfx1).transpose()
            tot_en_wfx2=pd.DataFrame(total_energy_wfx2).transpose()
            cc_a1d=pd.DataFrame(cc_a1).transpose() 
            cc_a2d=pd.DataFrame(cc_a2).transpose()
            with pd.ExcelWriter(path + '/' + SYS + '_energy.xlsx') as writer: #writer object needed to have the lists written in different sheets
            #also need the data to be a dataframe to enable writing to XLSX file
                tot_en_wfx1.to_excel(writer, sheet_name='E_WFX_OUT')
                tot_en_wfx2.to_excel(writer, sheet_name='E_WFX_IN')
            with pd.ExcelWriter(path + '/' + SYS + '_control_coordinates.xlsx') as writer: #writer object needed to have the lists written in different sheets
            #also need the data to be a dataframe to enable writing to XLSX file
                cc_a1d.to_excel(writer, sheet_name='cc_a1')
                cc_a2d.to_excel(writer, sheet_name='cc_a2')
        if FIT == True:
            tot_en_wfx1_fit=pd.DataFrame(total_energy_wfx1_fit).transpose()
            tot_en_wfx2_fit=pd.DataFrame(total_energy_wfx2_fit).transpose()
            cc_a1rd=pd.DataFrame(cc_a1_renorm).transpose() 
            cc_a2rd=pd.DataFrame(cc_a2_renorm).transpose()
            with pd.ExcelWriter(path + '/' + SYS + '_energy.xlsx') as writer: #writer object needed to have the lists written in different sheets
            #also need the data to be a dataframe to enable writing to XLSX file
                tot_en_wfx1_fit.to_excel(writer, sheet_name='E_WFX_FIT_OUT')
                tot_en_wfx2_fit.to_excel(writer, sheet_name='E_WFX_FIT_OUT')
            with pd.ExcelWriter(path + '/' + SYS + '_control_coordinates.xlsx') as writer: #writer object needed to have the lists written in different sheets
            #also need the data to be a dataframe to enable writing to XLSX file
                cc_a1rd.to_excel(writer, sheet_name='cc_a1_renorm')
                cc_a2rd.to_excel(writer, sheet_name='cc_a2_renorm')
        if DIFF == True:
            tot_en_wfx=pd.DataFrame(total_energy_wfx).transpose()
            cc_ad=pd.DataFrame(cc_a_sampled).transpose()
            with pd.ExcelWriter(path + '/' + SYS + '_energy.xlsx') as writer: #writer object needed to have the lists written in different sheets
            #also need the data to be a dataframe to enable writing to XLSX file
                tot_en_wfx.to_excel(writer, sheet_name='E_WFX_DIFF')
            with pd.ExcelWriter(path + '/' + SYS + '_control_coordinates.xlsx') as writer: #writer object needed to have the lists written in different sheets
            #also need the data to be a dataframe to enable writing to XLSX file
                cc_ad.to_excel(writer, sheet_name='cc_a_sampled')

    print("===================== PLOTTING SEGMENTS ===================")

    #if AUTO == True:
    #    critical_points = reg.find_critical(total_energy_wfx, cc_a, min_points=POINTS, use_inflex=INFLEX)
    #elif AUTO == False:
    #    critical_points = tp
    
    #Violin plots
    if NOFIT == True:
        rv.plot_violin([dataframe_list1[i]['R2'] for i in range(len(reg_inter1[0]))], save=SAVE_FIG, file_name = path1 +'/'+SYS1+'_violin.png',lang=lang)
        rv.plot_violin([dataframe_list2[i]['R2'] for i in range(len(reg_inter2[0]))], save=SAVE_FIG, file_name = path2 +'/'+SYS2+'_violin.png',lang=lang)
    if FIT == True:
        rv.plot_violin([dataframe_list1_fit[i]['R2'] for i in range(len(reg_inter1_fit[0]))], save=SAVE_FIG, file_name = path1 +'/'+SYS_f1+'_violin.png',lang=lang)
        rv.plot_violin([dataframe_list2_fit[i]['R2'] for i in range(len(reg_inter2_fit[0]))], save=SAVE_FIG, file_name = path2 +'/'+SYS_f2+'_violin.png',lang=lang)
    if DIFF == True:
        rv.plot_violin([dataframe_list[i]['R2'] for i in range(len(reg_inter[0]))], save=SAVE_FIG, file_name = path1 +'/'+SYS+'_violin.png',lang=lang)

    xlab=r'Reaction coordinate (a.m.u)'
    ylab=r'Relative Energy [$kJ$ $mol^{-1}$]'
    if NOFIT == True:
        rv.plot_segment(cc_a1, 2625.50*(total_energy_wfx1-total_energy_wfx1[0]), crit_list1, annotate=ANNOTATE, label=LABELS,
                    y_label = ylab, x_label = xlab, title = SYS1, size_title=SIZE_TITLE, size_label=SIZE_LABEL, save = SAVE_FIG, 
                    file_name=path1+'/'+SYS1+'REG_analysis_'+lang+'.png', force = FORCE, diff=False, lang=lang, ylim=(-90,120), verbose=verbose)   
        rv.plot_segment(cc_a2, 2625.50*(total_energy_wfx2-total_energy_wfx2[0]), crit_list2, annotate=ANNOTATE, label=LABELS,
                    y_label = ylab, x_label = xlab, title = SYS2, size_title=SIZE_TITLE, size_label=SIZE_LABEL, save = SAVE_FIG, 
                    file_name=path2+'/'+SYS2+'REG_analysis_'+lang+'.png', force = FORCE, diff=False, lang=lang, ylim=(-80,60),verbose=verbose)   
    if FIT == True:
        rv.plot_segment(cc_a_sampled, 2625.50*(total_energy_wfx1_fit-total_energy_wfx1_fit[0]), crit_list_fit, annotate=ANNOTATE, label=LABELS,
                    y_label = ylab, x_label = xlab, title = SYS_f1, size_title=SIZE_TITLE, size_label=SIZE_LABEL, save = SAVE_FIG, 
                    file_name=path+'/'+SYS+'REG_analysis_'+lang+'.png', force = FORCE, diff=False, lang=lang, verbose=verbose)   
        rv.plot_segment(cc_a_sampled, 2625.50*(total_energy_wfx2_fit-total_energy_wfx2_fit[0]), crit_list_fit, annotate=ANNOTATE, label=LABELS,
                    y_label = ylab, x_label = xlab, title = SYS_f2, size_title=SIZE_TITLE, size_label=SIZE_LABEL, save = SAVE_FIG, 
                    file_name=path+'/'+SYS+'REG_analysis_'+lang+'.png', force = FORCE, diff=False, lang=lang, verbose=verbose)   
    if DIFF == True:
        rv.plot_segment(cc_a_sampled, 2625.50*(total_energy_wfx-total_energy_wfx[0]), crit_list_fit, annotate=ANNOTATE, label=LABELS,
                    y_label = ylab, x_label = xlab, title = SYS, size_title=SIZE_TITLE, size_label=SIZE_LABEL, save = SAVE_FIG, 
                    file_name=path+'/'+SYS+'REG_analysis_'+lang+'.png', force = FORCE, diff=True, lang=lang, verbose=verbose)   
    #AC:Changed the energy plot by substracting the energy of the first point
    
    ###############################################################################
    #                                                                             #
    #                               GROUP TERMS (ATOMS and Vcl+Vx)                #
    #                                                                             #
    ###############################################################################
    
    print("======================= GROUPING TERMS =====================")

    if DOGROUP == True:
        #GROUP INTRA ATOMIC TERMS:     
        iqa_intra1_g = np.array(reg.group_intra(iqa_intra1, iqa_intra_header1, groups, verbose=verbose)[0])          #path1
        iqa_intra_header1_g = np.array(reg.group_intra(iqa_intra1, iqa_intra_header1, groups, verbose=verbose)[1])   #path1
        iqa_intra2_g = np.array(reg.group_intra(iqa_intra2, iqa_intra_header2, groups, verbose=verbose)[0])          #path2
        iqa_intra_header2_g = np.array(reg.group_intra(iqa_intra2, iqa_intra_header2, groups, verbose=verbose)[1])   #path2
            
        if DIFF == True:
            iqa_intra_g = np.array(reg.group_intra(iqa_intra, iqa_intra_header1, groups, verbose=verbose)[0])            #diff
            iqa_intra_header_g = np.array(reg.group_intra(iqa_intra, iqa_intra_header1, groups, verbose=verbose)[1])     #diff
            #INTRA ATOMIC CONTRIBUTION
            reg_intra_g = reg.reg(total_energy_wfx, cc_a_sampled, iqa_intra_g, np=POINTS, critical=AUTO_fit, 
                                  inflex=INFLEX, critical_index=crit_list_fit, verbose=verbose)
    
        #GROUP INTER ATOMIC TERMS: 
        iqa_inter1_g = np.array(reg.group_inter(iqa_inter1, iqa_inter_header1, groups, verbose=verbose)[0])          #path1
        iqa_inter_header1_g = np.array(reg.group_inter(iqa_inter1, iqa_inter_header1, groups, verbose=verbose)[1])   #path1
        iqa_inter2_g = np.array(reg.group_inter(iqa_inter2, iqa_inter_header2, groups, verbose=verbose)[0])          #path2
        iqa_inter_header2_g = np.array(reg.group_inter(iqa_inter2, iqa_inter_header2, groups, verbose=verbose)[1])   #path2
        
        #get signs of grouped terms
        sign_seg_intra1_g = reg.sign_segm(iqa_intra1_g,crit_list1, verbose=verbose)
        sign_seg_inter1_g = reg.sign_segm(iqa_inter1_g,crit_list1, verbose=verbose)
        sign_seg_disp1_g = reg.sign_segm(iqa_disp1,crit_list1, verbose=verbose)
        sign_seg1_g = [sign_seg_intra1_g,sign_seg_inter1_g,sign_seg_disp1_g]
        sign_seg_intra2_g = reg.sign_segm(iqa_intra2_g,crit_list2, verbose=verbose)
        sign_seg_inter2_g = reg.sign_segm(iqa_inter2_g,crit_list2, verbose=verbose)
        sign_seg_disp2_g = reg.sign_segm(iqa_disp2,crit_list2, verbose=verbose)
        sign_seg2_g = [sign_seg_intra2_g,sign_seg_inter2_g,sign_seg_disp2_g]
        sign_seg_intra1_ga = np.array(sign_seg_intra1_g)
        sign_seg_intra2_ga = np.array(sign_seg_intra2_g)
        sign_seg_intra_g = sign_seg_intra1_ga + sign_seg_intra2_ga
        sign_seg_intra_g = sign_seg_intra_g.tolist()
        sign_seg_inter1_ga = np.array(sign_seg_inter1_g)
        sign_seg_inter2_ga = np.array(sign_seg_inter2_g)
        sign_seg_inter_g = sign_seg_inter1_ga + sign_seg_inter2_ga
        sign_seg_inter_g = sign_seg_inter_g.tolist()
        sign_seg_disp1_ga = np.array(sign_seg_disp1_g)
        sign_seg_disp2_ga = np.array(sign_seg_disp2_g)
        sign_seg_disp_g = sign_seg_disp1_ga + sign_seg_disp2_ga
        sign_seg_disp_g = sign_seg_disp_g.tolist()
        
        sign_seg_g = [sign_seg_intra_g,sign_seg_inter_g,sign_seg_disp_g]
        if DIFF == True:
            iqa_inter_g = np.array(reg.group_inter(iqa_inter, iqa_inter_header1, groups, verbose=verbose)[0])            #diff
            iqa_inter_header_g = np.array(reg.group_inter(iqa_inter, iqa_inter_header1, groups, verbose=verbose)[1])     #diff
            #INTER ATOMIC CONTRIBUTION
            reg_inter_g = reg.reg(total_energy_wfx, cc_a_sampled, iqa_inter_g, np=POINTS, critical=AUTO_fit, inflex=INFLEX, 
                                  critical_index=crit_list_fit, verbose=verbose)
            dataframe_list_g=reg.csv(reg_inter_g, iqa_inter_header_g, reg_intra_g, iqa_intra_header_g, reg_disp, iqa_disp_header1,
                                     SYS, path, True, sign_seg_g, verbose=verbose)
            reg.excel(iqa_inter_g, iqa_inter_header_g, iqa_intra_g, iqa_intra_header_g, iqa_disp, iqa_disp_header1, SYS, path, True)
            rv.plot_violin([dataframe_list_g[i]['R2'] for i in range(len(reg_inter_g[0]))], save=SAVE_FIG, file_name = path +'/'+SYS+'_violin_g_'+lang+'.png', lang=lang)
    
    if DOGROUP2 == True:
        #GROUP VC AND VX
        groupedVcVx_1 = reg.group_vcvx(iqa_inter1_g, iqa_inter_header1_g, verbose=verbose)
        iqa_inter1_g2 = groupedVcVx_1[0]
        iqa_inter_header1_g2 = groupedVcVx_1[1]
        groupedVcVx_2 = reg.group_vcvx(iqa_inter2_g, iqa_inter_header2_g, verbose=verbose)
        iqa_inter2_g2 = groupedVcVx_2[0]
        iqa_inter_header2_g2 = groupedVcVx_2[1]
        sign_seg_inter1_g2=reg.sign_segm(iqa_inter1_g2, crit_list1, verbose=verbose)
        sign_seg1_g2 = [sign_seg_intra1_g, sign_seg_inter1_g2, sign_seg_disp1_g]
        sign_seg_inter2_g2=reg.sign_segm(iqa_inter2_g2, crit_list2, verbose=verbose)
        sign_seg2_g2 = [sign_seg_intra2_g, sign_seg_inter2_g2, sign_seg_disp2_g]
        if DIFF == True:
            groupedVcVx_g = reg.group_vcvx(iqa_inter_g, iqa_inter_header_g, verbose=verbose)
            iqa_inter_g2 = groupedVcVx_g[0]
            iqa_inter_header_g2 = groupedVcVx_g[1]
            sign_seg_inter1_g2a = np.array(sign_seg_inter1_g2)
            sign_seg_inter2_g2a = np.array(sign_seg_inter2_g2)
            sign_seg_inter_g2 = sign_seg_inter1_g2a + sign_seg_inter2_g2a
            sign_seg_inter_g2 = sign_seg_inter_g2.tolist()
            global sign_seg_g2
            sign_seg_g2 = [sign_seg_intra_g,sign_seg_inter_g2,sign_seg_disp_g]
            #INTER ATOMIC CONTRIBUTION
            reg_inter_g2 = reg.reg(total_energy_wfx, cc_a_sampled, iqa_inter_g2, np=POINTS, critical=AUTO_fit, inflex=INFLEX, 
                                   critical_index=crit_list_fit, verbose=verbose)
            pathG2=path+'/Vc+Vx_09'
            dataframe_list_g2=reg.csv(reg_inter_g2, iqa_inter_header_g2, reg_intra_g, iqa_intra_header_g, reg_disp, 
                                      iqa_disp_header1, SYS, pathG2, True, sign_seg_g2, verbose=verbose)
            reg.excel(iqa_inter_g2, iqa_inter_header_g2, iqa_intra_g, iqa_intra_header_g, iqa_disp, iqa_disp_header1, SYS, pathG2, True)
            rv.plot_violin([dataframe_list_g[i]['R2'] for i in range(len(reg_inter_g[0]))], save=SAVE_FIG, file_name = pathG2 +'/'+SYS+'_violin_g2_'+lang+'.png', lang=lang)

    ###############################################################################
    #                                                                             #
    #                           WRITE CSV FILES (GROUPED)                         #
    #                                                                             #
    ###############################################################################
    
    print("===================== WRITING CSV FILES ===================")
    
    if WRITE == True:
        if DIFF == True and DOGROUP == True:
            dataframe_list_g=reg.csv(reg_inter_g, iqa_inter_header_g, reg_intra_g, iqa_intra_header_g, reg_disp, iqa_disp_header1,
                                     SYS, path, True, sign_seg_g, verbose=verbose)
            reg.excel(iqa_inter_g,iqa_inter_header_g,iqa_intra_g,iqa_intra_header_g,iqa_disp,iqa_disp_header1,SYS,path,True)
            rv.plot_violin([dataframe_list_g[i]['R2'] for i in range(len(reg_inter_g[0]))], save=SAVE_FIG, file_name = path +'/'+SYS+'_violin_g.png', lang=lang)
        
        if DIFF == True and DOGROUP2 == True:
            pathG2=path+'/Vc+Vx_09'
            dataframe_list_g2=reg.csv(reg_inter_g2, iqa_inter_header_g2, reg_intra_g, iqa_intra_header_g, reg_disp, iqa_disp_header1,
                                      SYS, pathG2, True, sign_seg_g2, verbose=verbose)
            reg.excel(iqa_inter_g2,iqa_inter_header_g2,iqa_intra_g,iqa_intra_header_g,iqa_disp,iqa_disp_header1,SYS,pathG2,True)
            rv.plot_violin([dataframe_list_g[i]['R2'] for i in range(len(reg_inter_g[0]))], save=SAVE_FIG, file_name = pathG2 +'/'+SYS+'_violin_g2.png', lang=lang)
    
    
    ###############################################################################
    #                                                                             #
    #                           PLOT SEGMENTS AND TERMS                           #
    #                                                                             #
    ###############################################################################
    
    print("=================== PLOTTING SEGMENTS (2) ==================")

    if DETAILED_ANALYSIS ==True:
    
        if NOFIT == True: #original REG analysis for each pathway, no substraction or fitting
            for i in range(len(reg_inter1[0])): #AC change: R to R2.
                rv.generate_data_vis(dataframe_list1[i], [dataframe_list1[i]['R2'] for i in range(len(reg_inter1[0]))], 
                                                    n_terms, save=SAVE_FIG, file_name=path1+'/'+SYS1+'_detailed_seg_'+ str(i+1) +'_'+lang+'.png', 
                                                    title = SYS1 +' seg. '+str(i+1), size_title=SIZE_TITLE, size_label=SIZE_LABEL, color = COLOR_TERMS,
                                                    color_sign = COLOR_SIGN, sign_term = [sign_seg1[0][i],sign_seg1[1][i],sign_seg1[2][i]], lang=lang, verbose=verbose)
            for i in range(len(reg_inter2[0])): #AC change: R to R2.
                rv.generate_data_vis(dataframe_list2[i], [dataframe_list2[i]['R2'] for i in range(len(reg_inter2[0]))], 
                                                    n_terms, save=SAVE_FIG, file_name=path2+'/'+SYS2+'_detailed_seg_'+ str(i+1) +'_'+lang+'.png', 
                                                    title = SYS2 +' seg. '+str(i+1), size_title=SIZE_TITLE, size_label=SIZE_LABEL, color = COLOR_TERMS,
                                                    color_sign = COLOR_SIGN, sign_term = [sign_seg2[0][i],sign_seg2[1][i],sign_seg2[2][i]], lang=lang, verbose=verbose)
        if FIT == True: #REG analysis for fitted pathways, no subtraction
            for i in range(len(reg_inter1_fit[0])): #AC change: R to R2.
                rv.generate_data_vis(dataframe_list1_fit[i], [dataframe_list1_fit[i]['R2'] for i in range(len(reg_inter1_fit[0]))], 
                                                    n_terms, save=SAVE_FIG, file_name=path1+'/'+SYS_f1+'_detailed_seg_'+ str(i+1) +'_'+lang+'.png', 
                                                    title = SYS_f1 +' seg. '+str(i+1), size_title=SIZE_TITLE, size_label=SIZE_LABEL, color = COLOR_TERMS, 
                                                    color_sign = COLOR_SIGN, sign_term = [sign_seg1[0][i],sign_seg1[1][i],sign_seg1[2][i]], lang=lang, verbose=verbose)
            for i in range(len(reg_inter2_fit[0])): #AC change: R to R2.
                rv.generate_data_vis(dataframe_list2_fit[i], [dataframe_list2_fit[i]['R2'] for i in range(len(reg_inter2_fit[0]))], 
                                                    n_terms, save=SAVE_FIG, file_name=path2+'/'+SYS_f2+'_detailed_seg_'+ str(i+1) +'_'+lang+'.png', 
                                                    title = SYS_f2 +' seg. '+str(i+1), size_title=SIZE_TITLE, size_label=SIZE_LABEL, color = COLOR_TERMS, 
                                                    color_sign = COLOR_SIGN, sign_term = [sign_seg2[0][i],sign_seg2[1][i],sign_seg2[2][i]], lang=lang, verbose=verbose)
        if DIFF == True: #REG analysis for substracted profiles/terms
            for i in range(len(reg_inter[0])): #AC change: R to R2.
                rv.generate_data_vis(dataframe_list[i], [dataframe_list[i]['R2'] for i in range(len(reg_inter[0]))], 
                                                    n_terms, save=SAVE_FIG, file_name=path+'/'+SYS+'_detailed_seg_'+ str(i+1) +'_'+lang+'.png', 
                                                    title = SYS +' seg. '+str(i+1), size_title=SIZE_TITLE, size_label=SIZE_LABEL, color = COLOR_TERMS, 
                                                    color_sign = COLOR_SIGN, sign_term = [sign_seg[0][i],sign_seg[1][i],sign_seg[2][i]], lang=lang, verbose=verbose)
            if DOGROUP == True: #same but with grouped atoms
                if lang == "fr":
                    gr = "groupé"
                else:
                    gr = "grouped"
                for i in range(len(reg_inter_g[0])): #AC change: R to R2.
                    if i == 0 and len(crit_list1) == 2:    
                        rv.generate_data_vis(dataframe_list_g[i], [dataframe_list_g[i]['R2'] for i in range(len(reg_inter_g[0]))], 
                                                            n_terms, save=SAVE_FIG, file_name=path+'/'+SYS+'_detailed_seg_1+2_g_'+lang+'_3seg.png', 
                                                            title = SYS +' '+gr+' seg. '+str(i+1), size_title=SIZE_TITLE, size_label=SIZE_LABEL, color = COLOR_TERMS, 
                                                            color_sign = COLOR_SIGN, sign_term = [sign_seg_g[0][i],sign_seg_g[1][i],sign_seg_g[2][i]], lang=lang, verbose=verbose)
                    elif len(crit_list1) == 2:
                        rv.generate_data_vis(dataframe_list_g[i], [dataframe_list_g[i]['R2'] for i in range(len(reg_inter_g[0]))], 
                                                            n_terms, save=SAVE_FIG, file_name=path+'/'+SYS+'_detailed_seg_'+ str(i+2) +'_g_'+lang+'_3seg.png', 
                                                            title = SYS +' '+gr+' seg. '+str(i+1), size_title=SIZE_TITLE, size_label=SIZE_LABEL, color = COLOR_TERMS, 
                                                            color_sign = COLOR_SIGN, sign_term = [sign_seg_g[0][i],sign_seg_g[1][i],sign_seg_g[2][i]], lang=lang, verbose=verbose)
                    else:
                        rv.generate_data_vis(dataframe_list_g[i], [dataframe_list_g[i]['R2'] for i in range(len(reg_inter_g[0]))], 
                                                            n_terms, save=SAVE_FIG, file_name=path+'/'+SYS+'_detailed_seg_'+ str(i+1) +'_g_'+lang+'.png', 
                                                            title = SYS +' '+gr+' seg. '+str(i+1), size_title=SIZE_TITLE, size_label=SIZE_LABEL, color = COLOR_TERMS, 
                                                            color_sign = COLOR_SIGN, sign_term = [sign_seg_g[0][i],sign_seg_g[1][i],sign_seg_g[2][i]], lang=lang, verbose=verbose)
            if DOGROUP2 == True: #same but with grouped atoms and Vcl+Vx grouped into Inter terms
                
                    for i in range(len(reg_inter_g2[0])): #AC change: R to R2.
                        if len(crit_list1) == 2 and i == 0:    
                            rv.generate_data_vis(dataframe_list_g2[i], [dataframe_list_g2[i]['R2'] for i in range(len(reg_inter_g2[0]))], 
                                                            n_terms, save=SAVE_FIG, file_name=pathG2+'/'+SYS+'_detailed_seg_1+2_g2_'+lang+'_3seg.png', 
                                                            title = SYS +' '+gr+' seg. 1+2', size_title=SIZE_TITLE, size_label=SIZE_LABEL, color = COLOR_TERMS, 
                                                            color_sign = COLOR_SIGN, sign_term = [sign_seg_g2[0][i],sign_seg_g2[1][i],sign_seg_g2[2][i]], lang=lang, verbose=verbose)
                        elif len(crit_list1) == 2:
                            rv.generate_data_vis(dataframe_list_g2[i], [dataframe_list_g2[i]['R2'] for i in range(len(reg_inter_g2[0]))], 
                                                            n_terms, save=SAVE_FIG, file_name=pathG2+'/'+SYS+'_detailed_seg_'+ str(i+2) +'_g2_'+lang+'_3seg.png', 
                                                            title = SYS +' '+gr+' seg. '+str(i+2), size_title=SIZE_TITLE, size_label=SIZE_LABEL, color = COLOR_TERMS, 
                                                            color_sign = COLOR_SIGN, sign_term = [sign_seg_g2[0][i],sign_seg_g2[1][i],sign_seg_g2[2][i]], lang=lang, verbose=verbose)
                        else:
                            rv.generate_data_vis(dataframe_list_g2[i], [dataframe_list_g2[i]['R2'] for i in range(len(reg_inter_g2[0]))], 
                                                            n_terms, save=SAVE_FIG, file_name=pathG2+'/'+SYS+'_detailed_seg_'+ str(i+1) +'_g2_'+lang+'.png', 
                                                            title = SYS +' '+gr+' seg. '+str(i+1), size_title=SIZE_TITLE, size_label=SIZE_LABEL, color = COLOR_TERMS, 
                                                            color_sign = COLOR_SIGN, sign_term = [sign_seg_g2[0][i],sign_seg_g2[1][i],sign_seg_g2[2][i]], lang=lang, verbose=verbose)                                                 
                            
    print("=================== REGDIFF COMPLETED ==================")
    if verbose == True:
        sys.stdout.close()
