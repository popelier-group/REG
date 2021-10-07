# -*- coding: utf-8 -*-
"""
Created on Sat May 29 11:28:16 2021

@author: Ezan

REG-NBO running script accounting for G09/NBO3 (E(2) terms)
with procedures to recalculate the critical points of the compared profiles
according to the critical points of the difference of the profiles.

For a more recent version see V3
"""

# IMPORT MODULES
import reg
import aimall_utils
import nbo_utils
import numpy as np
# import numpy.polynomial.polynomial as poly
import matplotlib.pyplot as plt
import pandas as pd
import reg_vis as rv
from adjustText import adjust_text
import sys

#This code is under heavy rework so there are lots of commented things (sorry) - it should work fine as is
#it's basically a copy of the REGdiff script but with REGdiff-specific parts commented out
"""
##############################    OPTIONS    ##################################
# SYSTEM NAME
SYS0 = 'BH2'

# DEFINE PATHS AND FILES FOR PATH 1:
path1 = SYS0+"-out"  # Path to the folder with .wfn
# 155 #the Forward SP/AIM file numbers ; don't forget to add 5 to the last one (e.g. F0-F70 -> 0,75)
F1 = [i for i in range(0, 155, 1)]
# 70 #the Reverse SP/AIM file numbers ; don't forget to add 5 to the last one
R1 = [i for i in range(1, 67, 1)]

# DEFINE PATHS AND FILES FOR PATH 2:
path2 = SYS0+"-in"  # Path to the folder with .wfn
# the Forward SP/AIM file numbers ; don't forget to add 5 to the last one
F2 = [i for i in range(0, 129, 1)]
# the Reverse SP/AIM file numbers ; don't forget to add 5 to the last one
R2 = [i for i in range(1, 64, 1)]

# CALCULATION OPTIONS
# True if the file numbering is inverted for path 1 (Forward->Reactant and Reverse->Product)
inv_1 = False
# True if the file numbering is inverted for path 2 (Forward->Reactant and Reverse->Product)
inv_2 = False
groups = [['b10', 'h11', 'h12']]
POINTS = 4  # minimal number of points on each segment
INFLEX = True  # if the inflexion points have to be considered as critical points or not

# OUTPUT OPTIONS
WRITE = True  # if you want the script to write CSV and XLSX files with term values, REG values, energy values, control coordinate values
AUTO = True  # Automatically divide the energy profile into segments, or not
AUTO_fit = False  # Automatically divide the energy difference profile into segments, or not
FORCE = False  # If you want the force, or not
# If you want your term names to be colored in the segment analysis, or not
COLOR_TERMS = False
COLOR_SIGN = True  # If you want your term names to be colored according to the sign of the terms in the segment analysis, or not
SIZE_TITLE = 18  # fontsize of plot titles
SIZE_LABEL = 16  # fontsize of plot axis and tick labels
SAVE_FIG = True  # if you want to save plot figures to the disk
ANNOTATE = True  # if you want to write the segment number on the complete energy profile plot or not
# if you want to have the plots of significant terms for each segment or not
DETAILED_ANALYSIS = True
LABELS = True  # if you want to have point labels on the complete energy profile plot or not
# number of terms to be displayed on each segment graph (=number of significant terms)
n_terms = 5

# rd !!!!
# rdrm.REGdiff
# nb.REGNBO(SYS0,path1,F1,R1,path2,F2,R2,inv_1,inv_2,groups,POINTS,INFLEX,WRITE,AUTO,AUTO_fit,FORCE,COLOR_TERMS,COLOR_SIGN,SIZE_TITLE,
#            SIZE_LABEL,SAVE_FIG,ANNOTATE,DETAILED_ANALYSIS,LABELS,n_terms)
"""

def REGNBO(SYS0, path1, F1, R1, path2, F2, R2, inv_1, inv_2, groups, POINTS, INFLEX, WRITE, AUTO, AUTO_fit, FORCE, COLOR_TERMS, COLOR_SIGN, SIZE_TITLE,
            SIZE_LABEL, SAVE_FIG, ANNOTATE, DETAILED_ANALYSIS, LABELS, n_terms, lang="en"):
    path = SYS0+'-diff/'#test_C1+C2/'
    verbose = False #only for debugging purposes. Don't set to True unless you want to have a stream of blurb on your screen.
    if verbose == False:
        out = path+"out_NBO_noverbose.txt" 
    else:
        out = path+"out_NBO_verbose.txt"
        #sys.stdout = open(out, "w")
       
    #with open(out, 'w') as f:
        #with redirect_stdout(f):
    print("===================== STARTING REG-NBO V2 =====================")
    
    ###############################################################################
    #                                                                             #
    #                                 SETUP FILES                                 #
    #                                                                             #
    ###############################################################################

    print("====================== SETTING UP FILES ======================")
    
    # COMMON
    SYS1 = SYS0+'-out'
    SYS2 = SYS0+'-in'
    #SYS_f1 = SYS0+'-fit-out'
    #SYS_f2 = SYS0+'-fit-in'
    SYS = SYS0+'-diff (in-out)'
    path = SYS0+'-NBO/V2'

    # GET FILES AND CONTROL COORDINATES FOR PATH 1:
    res_get_list_wfx1 = aimall_utils.get_list_wfx(
        path1, SYS0, "out", F1, R1, inv_1, verbose=verbose)
    cc_a1 = res_get_list_wfx1[0]
    wfx_files1 = res_get_list_wfx1[3]
    out_files1 = res_get_list_wfx1[4]

    res_get_list_wfx2 = aimall_utils.get_list_wfx(
        path2, SYS0, "in", F2, R2, inv_2, verbose=verbose)
    cc_a2 = res_get_list_wfx2[0]
    wfx_files2 = res_get_list_wfx2[3]
    out_files2 = res_get_list_wfx2[4]

    # GET ATOM LIST FROM ANY .WFX FILE:
    NOFIT = True
    #FIT = True
    DIFF = True
    #DOGROUP = False
    #DOGROUP2 = False

    ###############################################################################
    #                                                                             #
    #                           GET TERMS AND ENERGIES                            #
    #                                                                             #
    ###############################################################################

    print("================= READING TERMS AND ENERGIES =================")
        
    # GET TOTAL ENERGY FROM THE .WFN FILES:
    total_energy_wfx1 = aimall_utils.get_aimall_wfx_energies(wfx_files1, verbose=verbose)
    total_energy_wfx2 = aimall_utils.get_aimall_wfx_energies(wfx_files2, verbose=verbose)
    
    nbo_prop1 = nbo_utils.nbo_property_from_g09_file(out_files1, verbose=verbose)
    nbo_prop2 = nbo_utils.nbo_property_from_g09_file(out_files2, verbose=verbose)
    #print(nbo_prop1[0][0])
    #print(nbo_prop1[1][0])
    if SYS0 == "BH2":
        print(nbo_prop1[0][0].index('5.BD(1)C2-C3/19.LP*(1)B10'))
        print(nbo_prop1[0][1].index('5.BD(1)C2-C3/19.LP*(1)B10'))
        print(nbo_prop2[0][0].index('5.BD(1)C2-C3/19.LP*(1)B10'))
        print(nbo_prop2[0][1].index('5.BD(1)C2-C3/19.LP*(1)B10'))
    if SYS0 == "NH2":
        print(nbo_prop1[0][0].index('19.LP(1)N10/197.BD*(1)C2-C3')) 
        print(nbo_prop1[0][1].index('19.LP(1)N10/197.BD*(1)C2-C3'))
        print(nbo_prop2[0][0].index('19.LP(1)N10/197.BD*(1)C2-C3'))
        print(nbo_prop2[0][1].index('19.LP(1)N10/197.BD*(1)C2-C3'))        
    # CONVERT LISTS TO ARRAYS.
    total_energy_wfx1 = np.array(total_energy_wfx1)

    # print(iqa_intra1[0])
    total_energy_wfx2 = np.array(total_energy_wfx2)

    ###############################################################################
    #                                                                             #
    #                               NORMALIZATION                                 #
    #                                                                             #
    ###############################################################################

    print("======================== NORMALIZATION =======================")

    # NORMALIZE THE CONTROL COORDINATES
    crit_list1 = reg.find_critical(total_energy_wfx1, cc_a1, use_inflex=True, verbose=verbose)
    crit_list2 = reg.find_critical(total_energy_wfx2, cc_a2, use_inflex=True, verbose=verbose)
    #crit_list1.pop(0)
    #crit_list2.pop(0)
    
    crit_trim = reg.crit_trim(crit_list1, crit_list2, SYS1, SYS2)
    crit_list1 = crit_trim[0]
    crit_list2 = crit_trim[1]

    cc_a1_renorm = reg.normalize(cc_a1,crit_list1)
    cc_a2_renorm = reg.normalize(cc_a2,crit_list2)

    rv.plot_energy(total_energy_wfx1,total_energy_wfx2,cc_a1,cc_a2,crit_list1,
                   crit_list2,SYS0,SIZE_TITLE,SIZE_LABEL,path,lang) #1 = OUT, 2 = IN
    rv.plot_renorm(total_energy_wfx1,total_energy_wfx2,cc_a1,cc_a2,cc_a1_renorm,
                   cc_a2_renorm,crit_list1,SYS0,SIZE_TITLE,SIZE_LABEL,path,lang) #1 = OUT, 2 = IN


    # GENERATE THE FULL CRITICAL POINT LIST (+EXTREMA)
    crit_list1_0N = [0]
    crit_list1_0N += crit_list1
    crit_list1_0N.append(len(cc_a1)-1)
    seg_long1 = [0]
    for i in range(len(crit_list1_0N)-1):
        seg_long1.append(crit_list1_0N[i+1]-crit_list1_0N[i])
    seg_long1.sort()
    maxlong1 = seg_long1[-1]

    crit_list2_0N = [0]
    crit_list2_0N += crit_list2
    crit_list2_0N.append(len(cc_a2)-1)
    seg_long2 = [0]
    for i in range(len(crit_list2_0N)-1):
        seg_long2.append(crit_list2_0N[i+1]-crit_list2_0N[i])
    seg_long2.sort()
    maxlong2 = seg_long2[-1]
    
    """
    if SYS0 == "BH2":
        terms = ['5.BD(1)C2-C3/19.LP*(1)B10'] 
        titles = ['BD(C2-C3) → LP*(B10)']
        file_name_terms = ['C2-C3_B10']
        
        comp=['5.BD(1)C2-C3/19.LP*(1)B10','9.BD(1)C3-H8/19.LP*(1)B10']
        plot0t2 = plt.figure(figsize=(10, 7))  # ENERGY OUT/IN ONLY
        plot1 = plot0t2.add_subplot(111)
        plot1.plot(nbo_int1[0][nbo_header1[0].index(comp[0])], nbo_int1[0][nbo_header1[0].index(comp[1])], 
                   'ro--', label=comp[0]+" "+comp[1]+' (out)')
        regr=reg.regression(nbo_int1[0][nbo_header1[0].index(comp[0])], nbo_int1[0][nbo_header1[0].index(comp[1])], "norm")[2]
        plot1.annotate("R² = "+str(regr**2)[:6],xy=(0.85,0.85),xycoords="axes fraction",fontsize=14)
        #plot1.plot(nbo_int2[0][nbo_prop2[0][0].index(comp[0])], nbo_int2[0][nbo_prop2[0][0].index(comp[1])], 
        #           color='orange', marker='o', linestyle='--', label=comp[0]+" "+comp[1]+' (in)')    
        #regr=reg.regression(nbo_int2[0][nbo_prop2[0][0].index(comp[0])], nbo_int2[0][nbo_prop2[0][0].index(comp[1])], "norm")[2]
        #plot1.annotate("R² = "+str(regr**2)[:6]+" (in)",xy=(0.85,0.80),xycoords="axes fraction",fontsize=14)
        #plot1.annotate(label+" out",xy=(nbo_int1[0][nbo_prop1[0].index(term)][-1],
        #                             total_energy_wfx1[crit_list1[0]+1]-total_energy_wfx1[0]))
        #plot1.annotate(label+" in",xy=(nbo_int2[0][nbo_prop2[0].index(term)][-1],
        #                            total_energy_wfx2[crit_list2[0]+1]-total_energy_wfx2[0]))
        plot1.legend(fontsize=SIZE_LABEL)
        plot1.set_title(comp[0]+"\n"+comp[1], fontsize=18)
        plot1.set_xlabel(comp[0]+' ($kJ$ $mol^{-1}$)', fontsize=16)
        plot1.margins(0.1,0.1)
        plot1.set_ylabel(comp[1]+' ($kJ$ $mol^{-1}$)', fontsize=16)
        plot1.tick_params(axis='both', labelsize=15)
        plot0t2.savefig(path+'/'+SYS0+'_comps')
        
    elif SYS0 == "NH2":
        terms = ['19.LP(1)N10/197.BD*(1)C2-C3','8.BD(1)C3-C4/196.BD*(1)C1-H5','1.BD(1)C1-C2/203.BD*(1)C4-H9']
        titles = ['LP(N10) → BD*(C2-C3)','BD(C3-C4) → BD*(C1-H5)','BD(C1-C2) → BD*(C4-H9)']
        file_name_terms = ['N10_C2-C3','C3-C4_C1-H5', 'C1-C2_C4-H9']
    elif SYS0 == "NO":
        terms = ['5.BD(1)C2-C3/211.BD*(1)N10-O11','5.BD(1)C2-C3/212.BD*(2)N10-O11']
        titles = ['BD(C2-C3) → BD*(N10-O11)']
        file_name_terms = ['C2-C3_N10-O11']
    if len(crit_list1) == 2 and len(crit_list2) == 2:
        label1 = "TS"
        for i in range(len(file_name_terms)):
            file_name_terms[i] += '_3seg'
            titles[i] += ' (seg. 1+2)'
    else:
        label1 = "Inflexion 1"
    i=0
    for term in terms:
        title = titles[i]
        file_name_term = file_name_terms[i]
        plot0t = plt.figure(figsize=(10, 7))  # OUT/IN 
        plot1 = plot0t.add_subplot(111)
        plot1.plot(cc_a1[:crit_list1[0]+1], (total_energy_wfx1[:crit_list1[0]+1] - total_energy_wfx1[0]), 'bo--', label='$\Delta E_{SCF}$ (out)')
        plot1.plot(cc_a1[:crit_list1[0]+1], nbo_int1[0][nbo_header1[0].index(term)], 'ro--', label='E(2) (out)')
        plot1.plot(cc_a2[:crit_list2[0]+1], (total_energy_wfx2[:crit_list2[0]+1]-total_energy_wfx2[0]), 
                   color='lightblue', marker='o', linestyle='--', label='$\Delta E_{SCF}$ (in)')
        plot1.plot(cc_a2[:crit_list2[0]+1], nbo_int2[0][nbo_header2[0].index(term)], 
                   color='orange', marker='o', linestyle='--', label='E(2) (in)')
        plot1.legend(fontsize=SIZE_LABEL)
        plot1.set_title(title, fontsize=18)
        plot1.set_xlabel('Reaction coordinate', fontsize=16)
        plot1.set_xlim(min(cc_a1[0],cc_a2[0])-0.1,max(cc_a1[crit_list1[0]+1],cc_a2[crit_list2[0]+1])+0.1) #-5.2, 0.1)
        plot1.set_ylabel('Energy ($kJ$ $mol^{-1}$)', fontsize=16)
        #plot1.set_ylim(-700, 150)
        plot1.tick_params(axis='both', labelsize=15)
        plot0t.savefig(path+'/'+SYS0+'_'+file_name_term)
        
        plot0t = plt.figure(figsize=(10, 7))  # OUT/IN RENORM
        plot1 = plot0t.add_subplot(111)
        plot1.plot(cc_a1_renorm[:crit_list1[0]+1], (total_energy_wfx1[:crit_list1[0]+1] - total_energy_wfx1[0]), 'bo--', label='$\Delta E_{SCF}$ (out)')
        plot1.plot(cc_a1_renorm[:crit_list1[0]+1], nbo_int1[0][nbo_header1[0].index(term)], 'ro--', label='E(2) (out)')
        plot1.plot(cc_a2_renorm[:crit_list2[0]+1], (total_energy_wfx2[:crit_list2[0]+1]-total_energy_wfx2[0]), 
                   color='lightblue', marker='o', linestyle='--', label='$\Delta E_{SCF}$ (in)')
        plot1.plot(cc_a2_renorm[:crit_list2[0]+1], nbo_int2[0][nbo_header2[0].index(term)], 
                   color='orange', marker='o', linestyle='--', label='E(2) (in)')
        plot1.legend(fontsize=SIZE_LABEL)
        plot1.set_title(title, fontsize=18)
        plot1.set_xlabel('Normalized reaction coordinate', fontsize=16)
        plot1.set_xlim(-0.1,1.1) #-5.2, 0.1)
        plot1.set_ylabel('Energy ($kJ$ $mol^{-1}$)', fontsize=16)
        #plot1.set_ylim(-700, 150)
        plot1.tick_params(axis='both', labelsize=15)
        plot0t.savefig(path+'/'+SYS0+'_'+file_name_term+'_renorm')
        
        if len(crit_list1) == 3 and len(crit_list2) == 3:
            plot0t = plt.figure(figsize=(10, 7))  # OUT/IN RENORM (segs 1 and 2)
            plot1 = plot0t.add_subplot(111)
            plot1.plot(cc_a1_renorm[:crit_list1[1]+1], (total_energy_wfx1[:crit_list1[1]+1] - total_energy_wfx1[0]), 'bo--', label='$\Delta E_{SCF}$ (out)')
            plot1.plot(cc_a1_renorm[:crit_list1[0]+1], nbo_int1[0][nbo_header1[0].index(term)], 'ro--')
            plot1.plot(cc_a1_renorm[crit_list1[0]:crit_list1[1]+1], nbo_int1[1][nbo_header1[1].index(term)], 'ro--', label='E(2) (out)')
            plot1.plot(cc_a2_renorm[:crit_list2[1]+1], (total_energy_wfx2[:crit_list2[1]+1]-total_energy_wfx2[0]), 
                       color='lightblue', marker='o', linestyle='--', label='$\Delta E_{SCF}$ (in)')
            plot1.plot(cc_a2_renorm[:crit_list2[0]+1], nbo_int2[0][nbo_header2[0].index(term)], 
                       color='orange', marker='o', linestyle='--')
            plot1.plot(cc_a2_renorm[crit_list2[0]:crit_list2[1]+1], nbo_int2[1][nbo_header2[1].index(term)], 
                       color='orange', marker='o', linestyle='--', label='E(2) (in)')
            plot1.legend(fontsize=SIZE_LABEL)
            plot1.set_title(title, fontsize=18)
            plot1.set_xlabel('Normalized reaction coordinate', fontsize=16)
            plot1.set_xlim(-0.1,2.1) #-5.2, 0.1)
            plot1.set_ylabel('Energy ($kJ$ $mol^{-1}$)', fontsize=16)
            #plot1.set_ylim(-700, 150)
            plot1.tick_params(axis='both', labelsize=15)
            plot0t.savefig(path+'/'+SYS0+'_'+file_name_term+'_renorm_0-2')
    """
    """
        plot0o = plt.figure(figsize=(10, 7))  # ENERGY OUT ONLY
        plot1 = plot0o.add_subplot(111)
        plot1.plot(cc_a1[:crit_list1[0]+1], (total_energy_wfx1[:crit_list1[0]+1] - total_energy_wfx1[0]), 'bo--', label='Energy (out)')
        plot1.plot(cc_a1[:crit_list1[0]+1], nbo_int1[0][nbo_header1[0].index(term)], 'ro--', label=term+' (out)')
        plot1.legend()
        plot1.set_title(title, fontsize=18)
        plot1.set_xlabel('Reaction coordinate', fontsize=16)
        plot1.set_xlim(cc_a1[0]-0.1,cc_a1[crit_list1[0]+1]+0.1)
        plot1.set_ylabel('$\Delta E_{SCF}$ ($kJ$ $mol^{-1}$)', fontsize=16)
        #plot1.set_ylim(-700, 150)
        plot1.tick_params(axis='both', labelsize=15)
        plot0o.savefig(SYS1+'/'+SYS1+'_'+file_name_term)
    
        plot0i = plt.figure(figsize=(10, 7))  # ENERGY IN ONLY
        plot1 = plot0i.add_subplot(111)
        plot1.plot(cc_a2[:crit_list2[0]+1], (total_energy_wfx2[:crit_list2[0]+1]-total_energy_wfx2[0]),
                   color='lightblue', marker='o', linestyle='--', label='Energy (in)')
        plot1.plot(cc_a2[:crit_list2[0]+1], nbo_int2[0][nbo_header2[0].index(term)], 
                   color='orange', marker='o', linestyle='--', label=term+' (in)')
        plot1.legend()
        plot1.set_title(title, fontsize=18)
        plot1.set_xlabel('Reaction coordinate', fontsize=16)
        plot1.set_xlim(cc_a2[0]-0.1,cc_a2[crit_list1[0]+1]+0.1)
        plot1.set_ylabel('$\Delta E_{SCF}$ ($kJ$ $mol^{-1}$)', fontsize=16)
        #plot1.set_ylim(-700, 150)
        plot1.tick_params(axis='both', labelsize=15)
        plot0i.savefig(SYS2+'/'+SYS2+'_'+file_name_term)
    """
    """
        if i == 0:
            plot0x = plt.figure(figsize=(10, 7))  # ENERGY OUT/IN ONLY
            plot1 = plot0x.add_subplot(111)
            plot1.plot((total_energy_wfx1[:crit_list1[0]+1]-total_energy_wfx1[0]), nbo_int1[0][nbo_header1[0].index(term)], 
                       'ro--', label='E(2) (out)')
            plot1.plot((total_energy_wfx2[:crit_list2[0]+1]-total_energy_wfx2[0]), nbo_int2[0][nbo_header2[0].index(term)],
                       color='orange', marker='o', linestyle='--', label='E(2) (in)')    
            
            if (nbo_int1[0][nbo_header1[0].index(term)][-1]-nbo_int1[0][nbo_header1[0].index(term)][0])/(total_energy_wfx1[crit_list1[0]+1]-total_energy_wfx1[0]) > 0: 
                offset1=1.1
            elif (nbo_int1[0][nbo_header1[0].index(term)][-1])/(total_energy_wfx1[crit_list1[0]+1]-total_energy_wfx1[0]) < 0: 
                offset1=0.9
            if (nbo_int2[0][nbo_header2[0].index(term)][-1]-nbo_int2[0][nbo_header2[0].index(term)][0])/(total_energy_wfx2[crit_list2[0]+1]-total_energy_wfx2[0]) > 0: 
                offset2=1.1
            elif (nbo_int2[0][nbo_header2[0].index(term)][-1])/(total_energy_wfx2[crit_list2[0]+1]-total_energy_wfx2[0]) < 0: 
                offset2=0.9
            
            #offset1,offset2 = 1, 1
            start1y = nbo_int1[0][nbo_header1[0].index(term)][0]
            end1x = offset1*(total_energy_wfx1[crit_list1[0]+1]-total_energy_wfx1[0])
            end1y = offset1*nbo_int1[0][nbo_header1[0].index(term)][-1]
            mid1x = end1x/2
            mid1y = (start1y+end1y)/2
            start2y = nbo_int2[0][nbo_header2[0].index(term)][0]
            end2x = offset2*(total_energy_wfx2[crit_list2[0]+1]-total_energy_wfx2[0])
            end2y = offset2*nbo_int2[0][nbo_header2[0].index(term)][-1]
            mid2x = end2x/2
            mid2y = (start2y+end2y)/2
            plot1.annotate(label1+" out",xy=(end1x,end1y), fontsize=16)
            plot1.annotate("Reactant out",xy=(0, start1y), fontsize=16)
            plot1.annotate(label1+" in",xy=(end2x,end2y), fontsize=16)
            plot1.annotate("Reactant in",xy=(0, start2y), fontsize=16)
            regr1=reg.regression((total_energy_wfx1[:crit_list1[0]+1]-total_energy_wfx1[0]), nbo_int1[0][nbo_header1[0].index(term)], "norm")[2]
            plot1.annotate("R² = "+str(regr1**2)[:6],xy=(mid1x,mid1y),fontsize=14)
            regr2=reg.regression((total_energy_wfx2[:crit_list2[0]+1]-total_energy_wfx2[0]), nbo_int2[0][nbo_header2[0].index(term)], "norm")[2]
            plot1.annotate("R² = "+str(regr2**2)[:6],xy=(mid2x,mid2y),fontsize=14)
            plot1.legend(fontsize=14,loc='best')
            plot1.set_title(title, fontsize=18)
            plot1.set_xlabel('$\Delta E_{SCF}$ ($kJ$ $mol^{-1}$)', fontsize=16)
            plot1.margins(0.2,0.2)
            plot1.set_ylabel('NBO E(2) ($kJ$ $mol^{-1}$)', fontsize=16)
            plot1.tick_params(axis='both', labelsize=15)
            plot0x.savefig(path+'/'+SYS0+'_'+file_name_term)
        i+=1
    """
    """
    for i in range(len(nbo_header1)):
        for j in range(len(nbo_header1[i])):
            #for k in range(len(nbo_header1[i][j])):
            for x in nbo_header1[i][j]:#[k]:
                if x == "/" :
                    rep=nbo_header1[i][j].split("/") 
                    nbo_header1[i][j] = rep[0] + "\n" + rep[1]
    for i in range(len(nbo_header2)):
        for j in range(len(nbo_header2[i])):
            for x in nbo_header2[i][j]:
                if x == "/" :
                    rep=nbo_header2[i][j].split("/") 
                    nbo_header2[i][j] = rep[0] + "\n" + rep[1]
    """
    ###############################################################################
    #                                                                             #
    #                               REG ANALYSIS                                  #
    #                                                                             #
    ###############################################################################
    
    #now with failsafes if REG values are null
    
    """
    if NOFIT == True:
        reg_nbo1 = reg.regnbo(total_energy_wfx1, cc_a1, nbo_int1, np=POINTS,
                          critical=AUTO_fit, inflex=INFLEX, critical_index=crit_list1)
        reg_nbo2 = reg.regnbo(total_energy_wfx2, cc_a2, nbo_int2, np=POINTS,
                          critical=AUTO_fit, inflex=INFLEX, critical_index=crit_list2)
        reg_trim1 = reg.reg_zeroes(reg_nbo1, [nbo_header1, nbo_int1, nbo_diff1, nbo_F1], verbose=verbose)
        reg_nbo1    = reg_trim1[0]
        nbo_header1 = reg_trim1[1][0]
        nbo_int1    = reg_trim1[1][1]
        nbo_diff1   = reg_trim1[1][2]
        nbo_F1      = reg_trim1[1][3]
        
        reg_trim2 = reg.reg_zeroes(reg_nbo2, [nbo_header2, nbo_int2, nbo_diff2, nbo_F2], verbose=verbose)
        reg_nbo2    = reg_trim2[0]
        nbo_header2 = reg_trim2[1][0]
        nbo_int2    = reg_trim2[1][1]
        nbo_diff2   = reg_trim2[1][2]
        nbo_F2      = reg_trim2[1][3]
    """   

    ###############################################################################
    #                                                                             #
    #                                   FITTING                                   #
    #                                                                             #
    ###############################################################################
    
    print("=========================== FITTING ==========================")
    
    #Critical point list for the fitting
    crit_list_fit_0N=[0]
    crit_list_fit_0N += crit_list2
    crit_list_fit_0N.append(len(cc_a2_renorm)-1)

    maxlong=max(maxlong1,maxlong2)

    crit_list_fit=[(maxlong-1)*i for i in range(1,len(crit_list1)+1)]
    print("crit_list1", crit_list1)
    #print(len(nbo_int1))
    
    #FITTING
    e1_cc = reg.fit(total_energy_wfx1, crit_list1_0N, maxlong, cc_a1_renorm, verbose=verbose)
    total_energy_wfx1_fit = e1_cc[0]
    cc_a_sampled = e1_cc[1]
    total_energy_wfx2_fit = reg.fit(total_energy_wfx2, crit_list2_0N, maxlong, cc_a2_renorm, verbose=verbose)[0]
    
    #print(len(total_energy_wfx1_fit),len(cc_a_sampled),len(nbo_int1_fit),len(nbo_int1_fit[0]),nbo_int1_fit)
    #print(len(total_energy_wfx2_fit),len(cc_a_sampled),len(nbo_int2_fit),len(nbo_int2_fit[0]),nbo_int2_fit)
    
    rv.plot_fit(total_energy_wfx1,total_energy_wfx2,total_energy_wfx1_fit,total_energy_wfx2_fit,cc_a1_renorm,cc_a2_renorm,
                cc_a_sampled,crit_list_fit,SYS0,SIZE_TITLE,SIZE_LABEL,path,lang)
    
    """
    if FIT == True:
        reg_nbo1_fit = reg.regnbo(total_energy_wfx1_fit, cc_a_sampled, nbo_int1_fit,
                                 np=POINTS, critical=AUTO_fit, inflex=INFLEX, critical_index=crit_list_fit)
        reg_nbo2_fit = reg.regnbo(total_energy_wfx2_fit, cc_a_sampled, nbo_int2_fit,
                                 np=POINTS, critical=AUTO_fit, inflex=INFLEX, critical_index=crit_list_fit)
        
        reg_trim1_fit = reg.reg_zeroes(reg_nbo1_fit, [nbo_header1, nbo_int1_fit], verbose=verbose)
        reg_nbo1_fit    = reg_trim1_fit[0]
        if NOFIT == False:
            nbo_header1 = reg_trim1_fit[1][0]
        nbo_int1_fit    = reg_trim1_fit[1][1]
        
        reg_trim2_fit = reg.reg_zeroes(reg_nbo2_fit, [nbo_header2, nbo_int2_fit], verbose=verbose)
        reg_nbo2_fit    = reg_trim2_fit[0]
        if NOFIT == False:
            nbo_header2 = reg_trim2_fit[1][0]
        nbo_int2_fit    = reg_trim2_fit[1][1]
    """
    ###############################################################################
    #                                                                             #
    #                                SUBSTRACTION (1)                             #
    #                                                                             #
    ###############################################################################
    
    print("========================= SUBSTRACTION ========================")
    
    #SUBSTRACT
    total_energy_wfx = total_energy_wfx2_fit - total_energy_wfx1_fit
    
    rv.plot_substract(total_energy_wfx1_fit,total_energy_wfx2_fit,total_energy_wfx,cc_a_sampled,crit_list_fit,
                      SYS0,SIZE_TITLE,SIZE_LABEL,path,lang)
    
    ###############################################################################
    #                                                                             #
    #                          DIFFERENCE CRITICAL POINTS                         #
    #                                                                             #
    ###############################################################################
    
    print("=================== DIFFERENCE CRITICAL POINTS ================")
        
    crit_list_diff = reg.find_critical(total_energy_wfx, cc_a_sampled, use_inflex=True, verbose=verbose)
    crit_list_diff2 = reg.find_critical(total_energy_wfx, cc_a_sampled, use_inflex=False, verbose=verbose)
    if verbose==True:
        print("crit_list_diff (with inflexion points)",crit_list_diff)
        print("crit_list_diff2 (without inflexion points)",crit_list_diff2)
    
    #Removing extra points from critical point lists
    if SYS0 == 'NH2':
        crit_list_diff2.pop(0)

    if SYS0 == 'CHO':
        crit_list_diff2.pop(3)
        crit_list_diff2.pop(0)
        
    if SYS0 == 'SiF3':
        crit_list_diff2.pop(-1)
        crit_list_diff2.pop(-1)
        crit_list_diff2.pop(0)
    
    if verbose==True:
        print("crit_list_diff2 (without inflexion points, trimmed)",crit_list_diff2)
       
    #cc_a_renorm = reg.normalize(cc_a_sampled,crit_list_diff2)

    crd = reg.norm_diff(cc_a1, cc_a2, cc_a_sampled, crit_list_diff2, crit_list1_0N, crit_list2_0N, verbose=verbose)
    crd1=crd[0]
    crd2=crd[1]
    if verbose==True:
        print("crd1",crd[0])
        print("crd2",crd[1])
    
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
    
    cc_a1_renorm = reg.normalize(cc_a1, crd[0], verbose=verbose)
    cc_a2_renorm = reg.normalize(cc_a2, crd[1], verbose=verbose)
    
    ###############################################################################
    #                                                                             #
    #                                 TRIMMING (1)                                #
    #                                                                             #
    ###############################################################################
    
    print("======================= TRIMMING TERMS ======================")
    
    nbo_trim1 = nbo_utils.nbo_trim(nbo_prop1, crit_list1_0N, verbose=verbose)
    nbo_header1 = nbo_trim1[0]
    nbo_int1 = nbo_trim1[1]
    
    if SYS0 == 'BH2' and verbose==True:
        print("Indices of term: 5.BD(1)C2-C3/19.LP*(1)B10")
        print("Path 1, segment 1",nbo_int1[0][nbo_header1[0].index('5.BD(1)C2-C3/19.LP*(1)B10')])
    if SYS0 == 'NH2' and verbose==True:
        print("Indices of term: 19.LP(1)N10/197.BD*(1)C2-C3")
        print("Path 1, segment 1",nbo_int1[0][nbo_header1[0].index('19.LP(1)N10/197.BD*(1)C2-C3')])
    nbo_diff1 = nbo_trim1[2]
    nbo_F1 = nbo_trim1[3]
    
    nbo_trim2 = nbo_utils.nbo_trim(nbo_prop2, crit_list2_0N, verbose=verbose) 
    nbo_header2 = nbo_trim2[0]
    nbo_int2 = nbo_trim2[1]
    
    if SYS0 == 'BH2' and verbose==True:
        print("Indices of term: 5.BD(1)C2-C3/19.LP*(1)B10")
        print("Path 2, segment 1",nbo_int2[0][nbo_header2[0].index('5.BD(1)C2-C3/19.LP*(1)B10')])
    if SYS0 == 'NH2' and verbose==True:
        print("Indices of term: 19.LP(1)N10/197.BD*(1)C2-C3")
        print("Path 2, segment 1",nbo_int2[0][nbo_header2[0].index('19.LP(1)N10/197.BD*(1)C2-C3')])
    nbo_diff2 = nbo_trim2[2]
    nbo_F2 = nbo_trim2[3]
    
    ###############################################################################
    #                                                                             #
    #                                 TRIMMING (2)                                #
    #                                                                             #
    ###############################################################################
    
    print("===================== TRIMMING TERMS (2) ====================")
    
    nbotrim = nbo_utils.nbo_trimdiff(nbo_header1, nbo_int1, nbo_diff1, nbo_F1, nbo_header2, nbo_int2, 
                                     nbo_diff2, nbo_F2, verbose=verbose)
    nbo_header_t = nbotrim[0]
    nbo_int1_t = nbotrim[1]
    nbo_int2_t = nbotrim[5]
    
    ###############################################################################
    #                                                                             #
    #                                   FITTING                                   #
    #                                                                             #
    ###############################################################################

    print("=========================== FITTING ==========================")
    
    maxlong=max(maxlong1,maxlong2)

    crit_list_fit=[(maxlong-1)*i for i in range(1,len(crit_list1_0N)-1)]
    if verbose==True:
        print("crit_list_fit",crit_list_fit)
    #print(len(nbo_int1))
    
    #FITTING
    e1_cc = reg.fit(total_energy_wfx1, crit_list1_0N, maxlong, cc_a1_renorm, verbose=verbose)
    total_energy_wfx1_fit = e1_cc[0]
    cc_a_sampled = e1_cc[1]
    total_energy_wfx2_fit = reg.fit(total_energy_wfx2, crit_list2_0N, maxlong, cc_a2_renorm, verbose=verbose)[0]
    
    #print(len(total_energy_wfx1_fit),len(cc_a_sampled),len(nbo_int1_fit),len(nbo_int1_fit[0]),nbo_int1_fit)
    #print(len(total_energy_wfx2_fit),len(cc_a_sampled),len(nbo_int2_fit),len(nbo_int2_fit[0]),nbo_int2_fit)
    
    rv.plot_fit(total_energy_wfx1,total_energy_wfx2,total_energy_wfx1_fit,total_energy_wfx2_fit,cc_a1_renorm,cc_a2_renorm,
                cc_a_sampled,crit_list_fit,SYS0,SIZE_TITLE,SIZE_LABEL,path,lang)
    
    nbo_int1_fit = reg.fit(nbo_int1_t, crit_list1_0N, maxlong, cc_a1_renorm, verbose=verbose)[0]
    nbo_int2_fit = reg.fit(nbo_int2_t, crit_list2_0N, maxlong, cc_a2_renorm, verbose=verbose)[0]
    
    ###############################################################################
    #                                                                             #
    #                                SUBSTRACTION                                 #
    #                                                                             #
    ###############################################################################
    
    print("========================= SUBSTRACTION ========================")
    
    #SUBSTRACT
    total_energy_wfx = total_energy_wfx2_fit - total_energy_wfx1_fit
    nbo_int = nbo_int2_fit - nbo_int1_fit
    
    #print(nbo_int[0][0])
    #print(iqa_intra[0])
    
    if verbose==True:
        print("check",len(total_energy_wfx),len(cc_a_sampled),len(nbo_int))
        for i in range(len(nbo_int)):
            print("Length of nbo_int list number",i,":",len(nbo_int[i]))
    
    rv.plot_substract(total_energy_wfx1_fit,total_energy_wfx2_fit,total_energy_wfx,cc_a_sampled,crit_list_fit,
                      SYS0,SIZE_TITLE,SIZE_LABEL,path,lang)
    
    total_energy_wfx1 = 2625.5*total_energy_wfx1
    total_energy_wfx2 = 2625.5*total_energy_wfx2
    total_energy_wfx1_fit = 2625.5*total_energy_wfx1_fit
    total_energy_wfx2_fit = 2625.5*total_energy_wfx2_fit
    total_energy_wfx = 2625.5*total_energy_wfx
    """
    if SYS0 == "BH2":
        term = '5.BD(1)C2-C3\n19.LP*(1)B10' 
    elif SYS0 == "NH2":
        term = '19.LP(1)N10\n197.BD*(1)C2-C3' 
    if len(crit_list1) == 2 :
        label = "TS"
    else:
        label = "Inflexion 1"
    plot0y = plt.figure(figsize=(10, 7))  # ENERGY OUT/IN ONLY
    plot1 = plot0y.add_subplot(111)
    plot1.plot((total_energy_wfx[:crit_list_fit[0]+1]-total_energy_wfx[0]), nbo_int[0][nbo_header_t[0].index(term)], 
               'ro--', label='E(2) (diff)')
    plot1.annotate(label,xy=(0.9*(total_energy_wfx[crit_list_fit[0]+1]-total_energy_wfx1[0]), 
                                    0.9*nbo_int[0][nbo_header_t[0].index(term)][-1]), fontsize=16)
    plot1.annotate("Reactant",xy=(0, nbo_int[0][nbo_header_t[0].index(term)][0]), fontsize=16)
    regr=reg.regression((total_energy_wfx[:crit_list_fit[0]+1]-total_energy_wfx[0]), nbo_int[0][nbo_header_t[0].index(term)], "norm")[2]
    plot1.annotate("R² = "+str(regr**2)[:6],xy=(0.05,0.80),xycoords="axes fraction",fontsize=14)
    plot1.legend(fontsize=SIZE_LABEL)
    plot1.set_title(titles[0], fontsize=18)
    plot1.set_xlabel('$\Delta E_{SCF}$ ($kJ$ $mol^{-1}$)', fontsize=16)
    plot1.margins(0.2,0.2)
    plot1.set_ylabel('NBO E(2) ($kJ$ $mol^{-1}$)', fontsize=16)
    plot1.tick_params(axis='both', labelsize=15)
    plot0y.savefig(path+'/'+SYS+'_'+file_name_term[0])
    
    plot0t = plt.figure(figsize=(10, 7))  # OUT/IN RENORM
    plot1 = plot0t.add_subplot(111)
    plot1.plot(cc_a_sampled[:crit_list_fit[0]+1], (total_energy_wfx[:crit_list_fit[0]+1]-total_energy_wfx[0]), 'bo--', label='$\Delta E_{SCF}$ (diff)')
    plot1.plot(cc_a_sampled[:crit_list_fit[0]+1], nbo_int[0][nbo_header_t[0].index(term)], 'ro--', label='E(2) (diff)')
    plot1.legend(fontsize=SIZE_LABEL)
    plot1.set_title(titles[0], fontsize=18)
    plot1.set_xlabel('Normalized reaction coordinate', fontsize=16)
    plot1.set_xlim(-0.1,1.1) #-5.2, 0.1)
    plot1.set_ylabel('Energy ($kJ$ $mol^{-1}$)', fontsize=16)
    #plot1.set_ylim(-700, 150)
    plot1.tick_params(axis='both', labelsize=15)
    plot0t.savefig(path+'/'+SYS+'_'+file_name_term[0]+'_renorm')
    
    """
    """
    rv.plot_substract_term(nbo_int1_fit,nbo_int1,nbo_int2_fit,nbo_int2,nbo_int,nbo_header1[0],0, 
                           cc_a1_renorm,cc_a2_renorm,cc_a_sampled,crit_list_fit,SYS0,SIZE_TITLE,SIZE_LABEL,path,lang) #not working yet
    #rv.plot_substract_term(iqa_inter1_fit,iqa_inter1,iqa_inter2_fit,iqa_inter2,iqa_inter,iqa_inter_header1,0,
    #                       cc_a1_renorm,cc_a2_renorm,cc_a_sampled,crit_list_fit,SYS0,SIZE_TITLE,SIZE_LABEL,path,lang) 
    """
    
    if DIFF == True:
        reg_nbo = reg.regnbo(total_energy_wfx, cc_a_sampled, nbo_int, np=POINTS,
                            critical=AUTO_fit, inflex=INFLEX, critical_index=crit_list_fit, verbose=verbose)
        reg_trim = reg.reg_zeroes(reg_nbo, [nbo_header1, nbo_int], verbose=verbose)
        reg_nbo = reg_trim[0]
        if NOFIT == False:
            nbo_header_t = reg_trim[1][0]
        nbo_int = reg_trim[1][1]
    
    ###############################################################################
    #                                                                             #
    #                     WRITE CSV FILES (WITHOUT GROUPS)                        #
    #                                                                             #
    ###############################################################################

    print("===================== WRITING CSV FILES ===================")
        
    #print("I",len(reg_nbo1[0]),len(nbo_header1[0]))
    if WRITE == True:
        """
        if NOFIT == True:
            if len(crit_list1) == 2:
                dataframe_list1=reg.csvnbo(reg_nbo1,nbo_header1,SYS1,path1,True)
                # for i in range(len(dataframe_list1)):
                    # print(dataframe_list1[i])
                dataframe_list2=reg.csvnbo(reg_nbo2,nbo_header2,SYS2,path2,True)
            else:
                dataframe_list1=reg.csvnbo(reg_nbo1,nbo_header1,SYS1,path1,False)
                # for i in range(len(dataframe_list1)):
                    # print(dataframe_list1[i])
                dataframe_list2=reg.csvnbo(reg_nbo2,nbo_header2,SYS2,path2,False)
    
        if FIT == True:
            dataframe_list1_fit=reg.csvnbo(reg_nbo1_fit,nbo_header_t,SYS_f1,path1,False)
            dataframe_list2_fit=reg.csvnbo(reg_nbo2_fit,nbo_header_t,SYS_f2,path2,False)
        """
        if DIFF == True:
            dataframe_list=reg.csvnbo(reg_nbo, nbo_header_t, SYS, path, False, verbose=verbose)
        # bloup
    
        """
        if NOFIT == True:
            reg.excel(iqa_inter1,iqa_inter_header1,iqa_intra1,iqa_intra_header1,iqa_disp1,iqa_disp_header1,SYS1,path1,False)
            reg.excel(iqa_inter2,iqa_inter_header2,iqa_intra2,iqa_intra_header2,iqa_disp2,iqa_disp_header2,SYS2,path2,False)
        if FIT == True:
            reg.excel(iqa_inter1_fit,iqa_inter_header1,iqa_intra1_fit,iqa_intra_header1,iqa_disp1_fit,iqa_disp_header1,SYS_f1,path1,False)
            reg.excel(iqa_inter2_fit,iqa_inter_header2,iqa_intra2_fit,iqa_intra_header2,iqa_disp2_fit,iqa_disp_header2,SYS_f2,path2,False)
        
        if DIFF == True:
            reg.excel(iqa_inter,iqa_inter_header1,iqa_intra,iqa_intra_header1,iqa_disp,iqa_disp_header1,SYS,path,False)
        
        # rv.plot_violin([dataframe_list[i]['R'] for i in range(len(reg_inter[0]))], save=SAVE_FIG, file_name = '../' +'violin.png',lang=lang)
        if NOFIT == True:
            rv.plot_violin([dataframe_list1[i]['R2'] for i in range(
                len(reg_inter1[0]))], save=SAVE_FIG, file_name = path1 +'/'+SYS1+'_violin.png',lang=lang)
            rv.plot_violin([dataframe_list2[i]['R2'] for i in range(
                len(reg_inter2[0]))], save=SAVE_FIG, file_name = path2 +'/'+SYS2+'_violin.png',lang=lang)
        if FIT == True:
            rv.plot_violin([dataframe_list1_fit[i]['R2'] for i in range(len(
                reg_inter1_fit[0]))], save=SAVE_FIG, file_name = path1 +'/'+SYS_f1+'_violin.png',lang=lang)
            rv.plot_violin([dataframe_list2_fit[i]['R2'] for i in range(len(
                reg_inter2_fit[0]))], save=SAVE_FIG, file_name = path2 +'/'+SYS_f2+'_violin.png',lang=lang)
        
        if DIFF == True:
            rv.plot_violin([dataframe_list[i]['R2'] for i in range(
                len(reg_inter[0]))], save=SAVE_FIG, file_name = path1 +'/'+SYS+'_violin.png',lang=lang)
        
        # GENERATE ENERGY AND CONTROL COORDINATE XLSX FILES
        if NOFIT == True:
            tot_en_wfx1=pd.DataFrame(total_energy_wfx1).transpose()
            tot_en_wfx2=pd.DataFrame(total_energy_wfx2).transpose()
            cc_a1d=pd.DataFrame(cc_a1).transpose()
            cc_a2d=pd.DataFrame(cc_a2).transpose()
            #writer object needed to have the lists written in different sheets
            with pd.ExcelWriter(path + '/' + SYS + '_energy.xlsx') as writer:
            # also need the data to be a dataframe to enable writing to XLSX file
                tot_en_wfx1.to_excel(writer, sheet_name='E_WFX_OUT')
                tot_en_wfx2.to_excel(writer, sheet_name='E_WFX_IN')
            #writer object needed to have the lists written in different sheets
            with pd.ExcelWriter(path + '/' + SYS + '_control_coordinates.xlsx') as writer:
            # also need the data to be a dataframe to enable writing to XLSX file
                cc_a1d.to_excel(writer, sheet_name='cc_a1')
                cc_a2d.to_excel(writer, sheet_name='cc_a2')
        if FIT == True:
            tot_en_wfx1_fit=pd.DataFrame(total_energy_wfx1_fit).transpose()
            tot_en_wfx2_fit=pd.DataFrame(total_energy_wfx2_fit).transpose()
            cc_a1rd=pd.DataFrame(cc_a1_renorm).transpose()
            cc_a2rd=pd.DataFrame(cc_a2_renorm).transpose()
            #writer object needed to have the lists written in different sheets
            with pd.ExcelWriter(path + '/' + SYS + '_energy.xlsx') as writer:
            # also need the data to be a dataframe to enable writing to XLSX file
                tot_en_wfx1_fit.to_excel(writer, sheet_name='E_WFX_FIT_OUT')
                tot_en_wfx2_fit.to_excel(writer, sheet_name='E_WFX_FIT_OUT')
            #writer object needed to have the lists written in different sheets
            with pd.ExcelWriter(path + '/' + SYS + '_control_coordinates.xlsx') as writer:
            # also need the data to be a dataframe to enable writing to XLSX file
                cc_a1rd.to_excel(writer, sheet_name='cc_a1_renorm')
                cc_a2rd.to_excel(writer, sheet_name='cc_a2_renorm')
        """
        if DIFF == True:
            tot_en_wfx=pd.DataFrame(total_energy_wfx).transpose()
            cc_ad=pd.DataFrame(cc_a_sampled).transpose()
            #writer object needed to have the lists written in different sheets
            with pd.ExcelWriter(path + '/' + SYS + '_energy.xlsx') as writer:
            # also need the data to be a dataframe to enable writing to XLSX file
                tot_en_wfx.to_excel(writer, sheet_name='E_WFX_DIFF')
            #writer object needed to have the lists written in different sheets
            with pd.ExcelWriter(path + '/' + SYS + '_control_coordinates.xlsx') as writer:
            # also need the data to be a dataframe to enable writing to XLSX file
                cc_ad.to_excel(writer, sheet_name='cc_a_sampled')
    """
    # if AUTO == True:
    #    critical_points = reg.find_critical(total_energy_wfx, cc_a, min_points=POINTS, use_inflex=INFLEX)
    # elif AUTO == False:
    #    critical_points = tp
    # blip
    """
    
    if len(crit_list1) == 2 or len(crit_list2) == 2 or len(crit_list_fit) == 2:
        name='REGNBO_analysis_3seg.png'
    else:
        name='REGNBO_analysis.png'
    """
    if NOFIT == True:
        
        rv.plot_segment(cc_a1, (total_energy_wfx1-total_energy_wfx1[0]), crit_list1, annotate=ANNOTATE, label=LABELS,
                        y_label = r'$\Delta E_{SCF}$ ($kJ$ $mol^{-1}$)', x_label=r'Reaction coordinate (a.m.u)', title = SYS1,
                        size_title=SIZE_TITLE, size_label=SIZE_LABEL, save = SAVE_FIG, file_name=path1+ \
                            '/'+SYS1+name, force = FORCE, diff=False, lang=lang) #_3seg
        rv.plot_segment(cc_a2, (total_energy_wfx2-total_energy_wfx2[0]), crit_list2, annotate=ANNOTATE, label=LABELS,
                        y_label = r'$\Delta E_{SCF}$ ($kJ$ $mol^{-1}$)', x_label=r'Reaction coordinate (a.m.u)', title = SYS2,
                        size_title=SIZE_TITLE, size_label=SIZE_LABEL, save = SAVE_FIG, file_name=path2+ \
                            '/'+SYS2+name, force = FORCE, diff=False, lang=lang) #_3seg
    if FIT == True:
        rv.plot_segment(cc_a_sampled, (total_energy_wfx1_fit-total_energy_wfx1_fit[0]), crit_list_fit, annotate=ANNOTATE, label=LABELS,
                    y_label = r'$\Delta E_{SCF}$ ($kJ$ $mol^{-1}$)', x_label=r'Reaction coordinate (a.m.u)', title = SYS_f1,
                    size_title=SIZE_TITLE, size_label=SIZE_LABEL, save = SAVE_FIG, file_name=path+ \
                        '/'+SYS+name, force = FORCE,diff=False, lang=lang)
        rv.plot_segment(cc_a_sampled, (total_energy_wfx2_fit-total_energy_wfx2_fit[0]), crit_list_fit, annotate=ANNOTATE, label=LABELS,
                    y_label = r'$\Delta E_{SCF}$ ($kJ$ $mol^{-1}$)', x_label=r'Reaction coordinate (a.m.u)', title = SYS_f2,
                    size_title=SIZE_TITLE, size_label=SIZE_LABEL, save = SAVE_FIG, file_name=path+ \
                        '/'+SYS+name, force = FORCE,diff=False, lang=lang)
    """
    if DIFF == True:
        rv.plot_segment(cc_a_sampled, (total_energy_wfx-total_energy_wfx[0]), crit_list_fit, annotate=ANNOTATE, label=LABELS,
                    y_label = r'$\Delta E_{SCF}$ ($kJ$ $mol^{-1}$)', x_label=r'Reaction coordinate (a.m.u)', title = SYS,
                    size_title=SIZE_TITLE, size_label=SIZE_LABEL, save = SAVE_FIG, file_name=path+ \
                        '/'+SYS+name, force = FORCE,diff=True, lang=lang)
    # AC:Changed the energy plot by substracting the energy of the first point
    
    ###############################################################################
    #                                                                             #
    #                           PLOT SEGMENTS AND TERMS                           #
    #                                                                             #
    ###############################################################################

    print("=================== PLOTTING SEGMENTS (2) ==================")
    
    #blip
    if DETAILED_ANALYSIS ==True:
        """
        if NOFIT == True:
             for i in range(len(dataframe_list2)): #AC change: R to R2.
                # print(dataframe_list2[i]['REG'])
                #len(dataframe_list2[i]['REG']) > 1 :
                if dataframe_list2[i]['REG'].dtype != np.dtype(object):
                    if i == 0 and len(crit_list2) == 2 :
                        rv.generate_data_vis(dataframe_list2[i], [dataframe_list2[i]['R2'] for i in range(len(nbo_int2))],
                                                    n_terms, save=SAVE_FIG, file_name=path1+'/'+SYS1+ '-NBO_detailed_seg_1+2.png',
                                                    title = SYS1 +' NBO seg. 1+2', color = COLOR_TERMS, lang=lang) #_3seg
                    else:
                        rv.generate_data_vis(dataframe_list2[i], [dataframe_list2[i]['R2'] for i in range(len(nbo_int2))],
                                                    n_terms, save=SAVE_FIG, file_name=path2+'/'+SYS2+ '-NBO_detailed_seg_'+ str(i+1) + '.png',
                                                    title = SYS2 +' NBO seg. '+str(i+1), color = COLOR_TERMS, lang=lang) #_3seg
            
             for i in range(len(dataframe_list1)): #AC change: R to R2.
                # print(dataframe_list1[i]['REG'])
                if dataframe_list1[i]['REG'].dtype != np.dtype(object) :
                    if i == 0 and len(crit_list1) == 2 :
                        rv.generate_data_vis(dataframe_list1[i], [dataframe_list1[i]['R2'] for i in range(len(nbo_int1))],
                                                    n_terms, save=SAVE_FIG, file_name=path1+'/'+SYS1+ '-NBO_detailed_seg_1+2.png',
                                                    title = SYS1 +' NBO seg. 1+2', color = COLOR_TERMS) #_3seg
                    else:
                        rv.generate_data_vis(dataframe_list1[i], [dataframe_list1[i]['R2'] for i in range(len(nbo_int1))],
                                                    n_terms, save=SAVE_FIG, file_name=path1+'/'+SYS1+ '-NBO_detailed_seg_'+ str(i+1) + '.png',
                                                    title = SYS1 +' NBO seg. '+str(i+1), color = COLOR_TERMS, lang=lang) #_3seg
          
        if FIT == True:
            for i in range(len(dataframe_list1_fit)): #AC change: R to R2.
                if i == 0 and len(crit_list1) == 2 :
                    rv.generate_data_vis(dataframe_list1_fit[i], [dataframe_list1_fit[i]['R2'] for i in range(len(nbo_int1_fit))],
                                                    n_terms, save=SAVE_FIG, file_name=path1+'/'+ SYS_f1+'-NBO_detailed_seg_1+2.png',
                                                    title = SYS_f1 +' NBO seg. 1+2', size_title=SIZE_TITLE, size_label=SIZE_LABEL, 
                                                    color = COLOR_TERMS, lang=lang)
                else:
                    rv.generate_data_vis(dataframe_list1_fit[i], [dataframe_list1_fit[i]['R2'] for i in range(len(nbo_int1_fit))],
                                                    n_terms, save=SAVE_FIG, file_name=path1+'/'+ SYS_f1+'-NBO_detailed_seg_'+ str(i+1) +'.png',
                                                    title = SYS_f1 +' NBO seg. '+str(i+1), size_title=SIZE_TITLE, size_label=SIZE_LABEL, 
                                                    color = COLOR_TERMS, lang=lang)
            for i in range(len(dataframe_list2_fit)): #AC change: R to R2.
                if i == 0 and len(crit_list1) == 2 :    
                    rv.generate_data_vis(dataframe_list2_fit[i], [dataframe_list2_fit[i]['R2'] for i in range(len(nbo_int2_fit))],
                                                            n_terms, save=SAVE_FIG, file_name=path2+'/'+ SYS_f2+'-NBO_detailed_seg_1+2.png',
                                                            title = SYS_f2 +' NBO seg. 1+2', size_title=SIZE_TITLE, size_label=SIZE_LABEL, 
                                                            color = COLOR_TERMS, lang=lang)
                else:
                    rv.generate_data_vis(dataframe_list2_fit[i], [dataframe_list2_fit[i]['R2'] for i in range(len(nbo_int2_fit))],
                                                            n_terms, save=SAVE_FIG, file_name=path2+'/'+ SYS_f2+'-NBO_detailed_seg_'+ str(i+1) +'.png',
                                                            title = SYS_f2 +' NBO seg. '+str(i+1), size_title=SIZE_TITLE, size_label=SIZE_LABEL, 
                                                            color = COLOR_TERMS, lang=lang)
        """
        if DIFF == True:
            for i in range(len(dataframe_list)): #AC change: R to R2.
                if i == 0 and len(crit_list1) == 2 :
                    rv.generate_data_vis(dataframe_list[i], [dataframe_list[i]['R2'] for i in range(len(nbo_int))],
                                                    n_terms, save=SAVE_FIG, file_name=path+'/'+ SYS+'-NBO_detailed_seg_1+2.png',
                                                    title = SYS +' NBO seg. 1+2', size_title=SIZE_TITLE, size_label=SIZE_LABEL, 
                                                    color = COLOR_TERMS, lang=lang, verbose=verbose)
                else:
                    rv.generate_data_vis(dataframe_list[i], [dataframe_list[i]['R2'] for i in range(len(nbo_int))],
                                                    n_terms, save=SAVE_FIG, file_name=path+'/'+ SYS+'-NBO_detailed_seg_'+ str(i+1) +'.png',
                                                    title = SYS +' NBO seg. '+str(i+1), size_title=SIZE_TITLE, size_label=SIZE_LABEL,
                                                    color = COLOR_TERMS, lang=lang, verbose=verbose)
    
    print("=================== REG-NBO V2 COMPLETED ==================")
    #if verbose == True:
        #sys.stdout.close() 

