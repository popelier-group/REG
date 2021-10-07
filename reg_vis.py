#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
reg_vis.py v0.0
L. J. Duarte, A. Cador, P. L. A. Popelier 

Library with function to visualize the results of the IQA REG analysis. 
Check for updates at github.com/ljduarte
For details about the method, please see XXXXXXX

Please, report bugs and issues to leo.j.duarte@hotmail.com.br and ael.cador@univ-rouen.fr
coded by L. J. Duarte
modified and expanded by A. Cador
"""

import numpy as np
import pandas as pd
from scipy.stats import norm
import matplotlib.pyplot as plt
#import matplotlib.image as mpimg
from adjustText import adjust_text
from matplotlib.ticker import MultipleLocator
import re  
#from matplotlib.offsetbox import (TextArea, DrawingArea, OffsetImage, AnnotationBbox)

"""
###########################################################################################################

PROCEDURE LIST

###########################################################################################################
plot_energy
            plot of the energy profiles of two paths (REGdiff)
plot_renorm
            plot of the renormalized energy profiles of two paths (REGdiff)
plot_fit
            plot of the renormalized energy profiles with fitted curves (REGdiff)
plot_substract
            plot of the renormalized, fitted energy profiles and of the difference in energy (REGdiff)
plot_substract_term
            plot of the fitted and unfitted term profiles (REGdiff)
set_axis_style
            set axis style for other graphs            
plot_violin
            violin plot of the REG analysis results (distribution of terms according to their R² coefficients)
generate_data_vis
            plots of the REG analysis for each segment
plot_segment
            plot of the segments (energy and force only)            
"""

"""
###########################################################################################################
(REGdiff)
FUNCTION: plot_energy
            plot of the energy profiles of two paths (REGdiff)

INPUT: total_energy_wfx1, total_energy_wfx2, cc_a1, cc_a2, crit_list1, crit_list2, SYS0, size_title, size_label, path, lang="en"
    total_energy_wfx1   = list of the WFX energy values along the 1st reaction path
    total_energy_wfx2   = list of the WFX energy values along the 2nd reaction path
    cc_a1               = list of the control coordinate of the 1st path
    cc_a2               = list of the control coordinate of the 2nd path
    crit_list1          = list of the critical points of the 1st path
    crit_list2          = list of the critical points of the 2nd path
    SYS0                = system name
    size_title          = int, fontsize of the plot title
    size_label          = int, fontsize of the plot labels
    path                = path to the folder where the plot is saved
    lang                = "en" (default) displays graph text in English
                        = "fr" displays graph text in French

OUTPUT: plot of the energy profiles

ERROR:
###########################################################################################################
"""

def plot_energy(total_energy_wfx1, total_energy_wfx2, cc_a1, cc_a2, crit_list1, crit_list2, SYS0, size_title, size_label,
                path, lang="en"):
    print("Starting to plot the energy profiles...")
    
    plot0 = plt.figure(figsize=(10,7))
    plot1 = plot0.add_subplot(111)
    if lang=="fr":
        plot1.plot(cc_a1,2625.5*(total_energy_wfx1-total_energy_wfx1[0]),'bo--',label='Énergie (out)')
        plot1.plot(cc_a2,2625.5*(total_energy_wfx2-total_energy_wfx2[0]),'ro--',label='Énergie (in)')
    else:
        plot1.plot(cc_a1,2625.5*(total_energy_wfx1-total_energy_wfx1[0]),'bo--',label='Energy (out)')
        plot1.plot(cc_a2,2625.5*(total_energy_wfx2-total_energy_wfx2[0]),'ro--',label='Energy (in)')
    plot1.axvline(cc_a1[0], linestyle = 'dotted', color = 'blue')
    plot1.axvline(cc_a1[-1], linestyle = 'dotted', color = 'blue')
    plot1.axvline(cc_a2[0], linestyle = 'dotted', color = 'red')
    plot1.axvline(cc_a2[-1], linestyle = 'dotted', color = 'red')
    for i in crit_list1: #split the segments
        plot1.axvline(cc_a1[i], linestyle = 'dotted', color = 'blue')
    for i in crit_list2: #split the segments
        plot1.axvline(cc_a2[i], linestyle = 'dotted', color = 'red')
    plot1.legend(fontsize=size_label)
    if lang=="fr":
        plot1.set_title('Normalisation de la coordonnée de réaction ('+SYS0+')',fontsize=size_title)
        plot1.set_xlabel(r'$\xi$',fontsize=size_label) #Coordonnée de réaction
        plot1.set_ylabel('Énergie relative ($kJ.mol^{-1}$)',fontsize=size_label)
    else:
        plot1.set_title('Reaction coordinate renormalization ('+SYS0+')',fontsize=size_title)
        plot1.set_xlabel('Reaction coordinate',fontsize=size_label)
        plot1.set_ylabel('Relative Energy ($kJ.mol^{-1}$)',fontsize=size_label)
    plot1.tick_params(axis='both',labelsize=size_label)
    plot0.savefig(path+'/'+SYS0+'_different_reaction_coordinates_'+lang)
    
    print("Finished plotting the energy profiles!")
    
"""
###########################################################################################################
(REGdiff)
FUNCTION: plot_renorm
            plot of the renormalized energy profiles of two paths (REGdiff)

INPUT: total_energy_wfx1, total_energy_wfx2, cc_a1, cc_a2, cc_a1_renorm, cc_a2_renorm, crit_list1, SYS0, size_title, size_label, path, lang="en"
    total_energy_wfx1   = list of the WFX energy values along the 1st reaction path
    total_energy_wfx2   = list of the WFX energy values along the 2nd reaction path
    cc_a1               = list of the control coordinate of the 1st path
    cc_a2               = list of the control coordinate of the 2nd path
    cc_a1_renorm        = list of the normalized control coordinate of the 1st path
    cc_a2_renorm        = list of the normalized control coordinate of the 2nd path
    crit_list1          = list of the critical points of the 1st path
    SYS0                = system name
    size_title          = int, fontsize of the plot title
    size_label          = int, fontsize of the plot labels
    path                = path to the folder where the plot is saved
    lang                = "en" (default) displays graph text in English
                        = "fr" displays graph text in French

OUTPUT: plot of the renormalized energy profiles

ERROR:
###########################################################################################################
"""

def plot_renorm(total_energy_wfx1, total_energy_wfx2, cc_a1, cc_a2, cc_a1_renorm, cc_a2_renorm, crit_list1, SYS0, size_title,
                size_label, path, lang="en"):
    
    print("Starting to plot the normalized energy profiles...")
    
    plot = plt.figure(figsize=(10,7))
    plot11 = plot.add_subplot(111)
    if lang=="fr":
        plot11.plot(cc_a1,2625.5*(total_energy_wfx1-total_energy_wfx1[0]),'bo--',label='Énergie (out)')
        plot11.plot(cc_a2,2625.5*(total_energy_wfx2-total_energy_wfx2[0]),'ro--',label='Énergie (in)')
    else:
        plot11.plot(cc_a1,2625.5*(total_energy_wfx1-total_energy_wfx1[0]),'bo--',label='Energy (out)')
        plot11.plot(cc_a2,2625.5*(total_energy_wfx2-total_energy_wfx2[0]),'ro--',label='Energy (in)')
    if lang=="fr":
        plot11.set_title('Normalisation de la coordonnée de réaction ('+SYS0+')',fontsize=size_title)
        plot11.set_xlabel(r'$\xi$ (a.m.u.)',fontsize=size_label) #Coordonnée de réaction
        plot11.set_ylabel('Énergie relative ($kJ.mol^{-1}$)',fontsize=size_label)
    else:
        plot11.set_title('Reaction coordinate renormalization ('+SYS0+')',fontsize=size_title)
        plot11.set_xlabel('Reaction coordinate',fontsize=size_label)
        plot11.set_ylabel('Relative Energy ($kJ.mol^{-1}$)',fontsize=size_label)
    plot11.tick_params(axis='both',labelsize=size_label)
    plot11.axvline(cc_a1_renorm[0], linestyle = 'dotted', color = 'black')
    plot11.axvline(cc_a1_renorm[-1], linestyle = 'dotted', color = 'black')
    for i in crit_list1: #split the segments
        plot11.axvline(cc_a1_renorm[i], linestyle = 'dotted', color = 'black')
    if lang=="fr":
        plot11.plot(cc_a1_renorm,2625.5*(total_energy_wfx1-total_energy_wfx1[0]),'o--',label='Énergie (out, renorm.)')
        plot11.plot(cc_a2_renorm,2625.5*(total_energy_wfx2-total_energy_wfx2[0]),'o--',label='Énergie (in, renorm.)')
    else:
        plot11.plot(cc_a1_renorm,2625.5*(total_energy_wfx1-total_energy_wfx1[0]),'o--',label='Energy (out, renorm.)')
        plot11.plot(cc_a2_renorm,2625.5*(total_energy_wfx2-total_energy_wfx2[0]),'o--',label='Energy (in, renorm.)')
    plot11.legend(fontsize=size_label)
    plot.savefig(path+'/'+SYS0+'_reaction_coordinate_renormalization_'+lang)
    
    print("Finished plotting the normalized energy profiles!")
    return

"""
###########################################################################################################
(REGdiff)
FUNCTION: plot_fit
            plot of the renormalized energy profiles with fitted curves (REGdiff)
            
INPUT: total_energy_wfx1, total_energy_wfx2, total_energy_wfx1_fit, total_energy_wfx2_fit, cc_a1_renorm, cc_a2_renorm, 
             cc_a_sampled, crit_list_fit, SYS0, size_title, size_label, path, lang="en"
    total_energy_wfx1       = list of the WFX energy values along the 1st reaction path
    total_energy_wfx2       = list of the WFX energy values along the 2nd reaction path
    total_energy_wfx1_fit   = list of the fitted WFX energy values along the 1st reaction path
    total_energy_wfx2_fit   = list of the fitted WFX energy values along the 2nd reaction path
    cc_a1_renorm            = list of the normalized control coordinate of the 1st path
    cc_a2_renorm            = list of the normalized control coordinate of the 2nd path
    cc_a_sampled            = list of the sampled, normalized control coordinate (common to both paths)
    crit_list_fit           = list of the critical points of the control coordinate after fitting
    SYS0                    = system name
    size_title              = int, fontsize of the plot title
    size_label              = int, fontsize of the plot labels
    path                    = path to the folder where the plot is saved
    lang                    = "en" (default) displays graph text in English
                            = "fr" displays graph text in French

OUTPUT: plot of the renormalized energy profiles with fitted curves

ERROR:
###########################################################################################################
"""

def plot_fit(total_energy_wfx1, total_energy_wfx2, total_energy_wfx1_fit, total_energy_wfx2_fit, cc_a1_renorm, cc_a2_renorm, 
             cc_a_sampled, crit_list_fit, SYS0, size_title, size_label, path, lang="en"):
    print("Starting to plot the fitted and unfitted energy profiles...")
    
    plot_ = plt.figure(figsize=(10,7))
    plot2 = plot_.add_subplot(111)

    if lang=="fr":
        plot2.plot(cc_a1_renorm,2625.5*(total_energy_wfx1-total_energy_wfx1[0]),color='deepskyblue',linestyle='--',marker='o',label='Énergie (out)') #markeredgecolor='darkblue'
        plot2.plot(cc_a2_renorm,2625.5*(total_energy_wfx2-total_energy_wfx2[0]),color='orange',linestyle='--',marker='o',label='Énergie (in)') #,markeredgecolor='red'
    else:
        plot2.plot(cc_a1_renorm,2625.5*(total_energy_wfx1-total_energy_wfx1[0]),color='deepskyblue',linestyle='--',marker='o',label='Energy (out)') #markeredgecolor='darkblue'
        plot2.plot(cc_a2_renorm,2625.5*(total_energy_wfx2-total_energy_wfx2[0]),color='orange',linestyle='--',marker='o',label='Energy (in)') #,markeredgecolor='red'
    plot2.legend(fontsize=size_label)
    if lang=="fr":
        plot2.set_title('Rééchantillonnage de l\'énergie ('+SYS0+')',fontsize=size_title)
        plot2.set_xlabel(r'$\tilde \xi$ (a.m.u.)',fontsize=size_label) #Coordonnée de réaction normalisée
        plot2.set_ylabel('$\Delta E_{mol}$ ($kJ.mol^{-1}$)',fontsize=size_label)
    else:
        plot2.set_title('Energy fitting ('+SYS0+')',fontsize=size_title)
        plot2.set_xlabel('Reaction coordinate',fontsize=size_label)
        plot2.set_ylabel('Relative Energy ($kJ$ $mol^{-1}$)',fontsize=size_label)
    plot2.tick_params(axis='both',labelsize=15)
    plot2.set_xticks([i for i in range(0,len(crit_list_fit)+2)])
    plot2.grid(b=True,axis='x')
    #plot_.savefig(path+'/'+SYS0+'_energy_fitting')
    plot2.plot(cc_a_sampled[0:crit_list_fit[0]+1],2625.5*(total_energy_wfx1_fit[0:crit_list_fit[0]+1]-total_energy_wfx1_fit[0]),
               color='blue',linestyle='--',marker='.')
    plot2.plot(cc_a_sampled[0:crit_list_fit[0]+1],2625.5*(total_energy_wfx2_fit[0:crit_list_fit[0]+1]-total_energy_wfx2_fit[0]),
               color='red',linestyle='--',marker='.')
    #plot_.savefig(path+'/'+SYS0+'_energy_fitting_0')
    for i in range(len(crit_list_fit)-1):
        plot2.plot(cc_a_sampled[crit_list_fit[i]:crit_list_fit[i+1]+1],2625.5*(total_energy_wfx1_fit[crit_list_fit[i]:crit_list_fit[i+1]+1]-total_energy_wfx1_fit[0]),
                   color='blue',linestyle='-')
        plot2.plot(cc_a_sampled[crit_list_fit[i]:crit_list_fit[i+1]+1],2625.5*(total_energy_wfx2_fit[crit_list_fit[i]:crit_list_fit[i+1]+1]-total_energy_wfx2_fit[0]),
                   color='red',linestyle='-')
        plot_.savefig(path+'/'+SYS0+'_energy_fitting_'+str(i+1)+'_'+lang)
    if lang == "fr":
        plot2.plot(cc_a_sampled[crit_list_fit[-1]:],2625.5*(total_energy_wfx1_fit[crit_list_fit[-1]:]-total_energy_wfx1_fit[0]),
                   color='blue',linestyle='-',label='Énergie rééchantillonnée (out)')
        plot2.plot(cc_a_sampled[crit_list_fit[-1]:],2625.5*(total_energy_wfx2_fit[crit_list_fit[-1]:]-total_energy_wfx2_fit[0]),
                   color='red',linestyle='-',label='Énergie rééchantillonnée (in)')
    else:
        plot2.plot(cc_a_sampled[crit_list_fit[-1]:],2625.5*(total_energy_wfx1_fit[crit_list_fit[-1]:]-total_energy_wfx1_fit[0]),
                   color='blue',linestyle='-',label='Energy fit (out)')
        plot2.plot(cc_a_sampled[crit_list_fit[-1]:],2625.5*(total_energy_wfx2_fit[crit_list_fit[-1]:]-total_energy_wfx2_fit[0]),
                   color='red',linestyle='-',label='Energy fit (in)')
    plot_.savefig(path+'/'+SYS0+'_energy_fitting_N_'+lang)
    
    print("Finished plotting the fitted and unfitted energy profiles!")
    return

"""
###########################################################################################################
(REGdiff)
FUNCTION: plot_substract
            plot of the renormalized, fitted energy profiles and of the difference in energy (REGdiff)

INPUT: total_energy_wfx1_fit, total_energy_wfx2_fit, total_energy_wfx, cc_a_sampled, crit_list_fit, SYS0, size_title,
                   size_label, path, lang="en"
    total_energy_wfx1_fit   = list of the fitted WFX energy values along the 1st reaction path
    total_energy_wfx2_fit   = list of the fitted WFX energy values along the 2nd reaction path
    total_energy_wfx        = list of the differences of the fitted WFX energy values along each reaction path    
    cc_a1                   = list of the control coordinate of the 1st path
    cc_a2                   = list of the control coordinate of the 2nd path
    cc_a_sampled            = list of the sampled, normalized control coordinate (common to both paths)
    crit_list_fit           = list of the critical points of the control coordinate after fitting
    SYS0                    = system name
    size_title              = int, fontsize of the plot title
    size_label              = int, fontsize of the plot labels
    path                    = path to the folder where the plot is saved
    lang                    = "en" (default) displays graph text in English
                            = "fr" displays graph text in French

OUTPUT: plot of the renormalized, fitted energy profiles and of the difference in energy

ERROR:
###########################################################################################################
"""

def plot_substract(total_energy_wfx1_fit, total_energy_wfx2_fit, total_energy_wfx, cc_a_sampled, crit_list_fit, SYS0, size_title,
                   size_label, path, lang="en"):
    print("Starting to plot the fitted, unfitted and substracted energy profiles...")
    
    plot_1 = plt.figure(figsize=(10,7))
    plot20 = plot_1.add_subplot(111)
    if lang == "fr":
        plot20.plot(cc_a_sampled,2625.5*(total_energy_wfx1_fit-total_energy_wfx1_fit[0]),color='blue',linestyle='-',marker='.',label='Énergie rééchantillonnée (out)')
        plot20.plot(cc_a_sampled,2625.5*(total_energy_wfx2_fit-total_energy_wfx2_fit[0]),color='red',linestyle='-',marker='.',label='Énergie rééchantillonnée (in)')
        plot20.set_title('Rééchantillonnage de l\'énergie ('+SYS0+')',fontsize=size_title)
        plot20.set_xlabel(r'$\tilde \xi$ (a.m.u.)',fontsize=size_label) #Coordonnée de réaction normalisée
        plot20.set_ylabel('$\Delta E_{mol}$ ($kJ.mol^{-1}$)',fontsize=size_label)
        plot20.set_ylim(-90,120)
        plot20.plot(cc_a_sampled,2625.5*(total_energy_wfx-total_energy_wfx[0]),color='limegreen',linestyle='-',marker='.',label='Énergie (différence in-out)')
    else:    
        plot20.plot(cc_a_sampled,(total_energy_wfx1_fit-total_energy_wfx1_fit[0]),color='blue',linestyle='-',marker='.',label='Energy fit (out)')
        plot20.plot(cc_a_sampled,(total_energy_wfx2_fit-total_energy_wfx2_fit[0]),color='red',linestyle='-',marker='.',label='Energy fit (in)')
        plot20.set_title('Energy fitting ('+SYS0+')',fontsize=size_title)
        plot20.set_xlabel('Normalized reaction coordinate',fontsize=size_label)
        plot20.set_ylabel('Relative Energy ($kJ$ $mol^{-1}$)',fontsize=size_label)
        plot20.plot(cc_a_sampled,(total_energy_wfx-total_energy_wfx[0]),color='limegreen',linestyle='-',marker='.',label='Energy (in-out difference)')
    plot20.tick_params(axis='both',labelsize=size_label)
    plot20.set_xticks([i for i in range(0,len(crit_list_fit)+2)])
    plot20.grid(b=True,axis='x')
    plot20.legend(fontsize=size_label-2)
    plot_1.savefig(path+'/'+SYS0+'_energy_fitting_and_subtraction_'+lang)
    
    print("Finished plotting the fitted, unfitted and substracted energy profiles!")
    return

"""
###########################################################################################################
(REGdiff)
FUNCTION: plot_substract_term
            plot of the fitted and unfitted term profiles (REGdiff)
            

INPUT: term1_fit, term1, term2_fit, term2, term, term_header, index, cc_a1_renorm, cc_a2_renorm,
                        cc_a_sampled, crit_list_fit, SYS0, size_title, size_label, path, lang="en"
    term1_fit               = list of the fitted term values along the 1st reaction path (all terms)
    term1                   = list of the original term values along the 1st reaction path (all terms)
    term2_fit               = list of the fitted term values along the 2nd reaction path (all terms)
    term2                   = list of the original term values along the 2nd reaction path (all terms)
    term                    = list of the differences of the fitted WFX energy values along each reaction path (all terms)    
    term_header             = list of the term names    
    index                   = index of the term to plot    
    cc_a1_renorm            = list of the renormalized control coordinate of the 1st path
    cc_a2_renorm            = list of the renormalized control coordinate of the 2nd path
    cc_a_sampled            = list of the sampled, normalized control coordinate (common to both paths)
    crit_list_fit           = list of the critical points of the control coordinate after fitting
    SYS0                    = system name
    size_title              = int, fontsize of the plot title
    size_label              = int, fontsize of the plot labels
    path                    = path to the folder where the plot is saved
    lang                    = "en" (default) displays graph text in English
                            = "fr" displays graph text in French

OUTPUT: plot of the fitted and unfitted term profiles

ERROR:
    "Index out of range for term number": the index entered doesn't correspond to the list
###########################################################################################################
"""

def plot_substract_term(term1_fit, term1, term2_fit, term2, term, term_header, index, cc_a1_renorm, cc_a2_renorm,
                        cc_a_sampled, crit_list_fit, SYS0, size_title, size_label, path, lang="en"):
    print("Starting to plot the fitted, unfitted and substracted term profiles...")
    if index > len(term1)-1:
        raise ValueError("Index out of range for term number")
    plot__ = plt.figure(figsize=(10,7))
    plot3 = plot__.add_subplot(111)
    if lang == "fr":
        plot3.plot(cc_a_sampled,2625.5*(term1_fit[index]-term1_fit[index][0]),color='blue',linestyle='-',marker='.',label=term_header[0]+' (out) ré-éch.')
        plot3.plot(cc_a_sampled,2625.5*(term2_fit[index]-term2_fit[index][0]),color='red',linestyle='-',marker='.',label=term_header[0]+' (in) ré-éch.')
        plot3.plot(cc_a_sampled,2625.5*(term[index]-term[index][0]),color='limegreen',linestyle='-',marker='.',label=term_header[0]+' (différence in-out)')
    else:
        plot3.plot(cc_a_sampled,2625.5*(term1_fit[index]-term1_fit[index][0]),color='blue',linestyle='-',marker='.',label=term_header[0]+' (out) fit')
        plot3.plot(cc_a_sampled,2625.5*(term2_fit[index]-term2_fit[index][0]),color='red',linestyle='-',marker='.',label=term_header[0]+' (in) fit')
        plot3.plot(cc_a_sampled,2625.5*(term[index]-term[index][0]),color='limegreen',linestyle='-',marker='.',label=term_header[0]+' (in-out difference)')
    plot3.plot(cc_a1_renorm,2625.5*(term1[index]-term1[index][0]),color='deepskyblue',linestyle='-',marker='o',label=term_header[0]+' (out)')
    plot3.plot(cc_a2_renorm,2625.5*(term2[index]-term2[index][0]),color='orange',linestyle='-',marker='o',label=term_header[0]+' (in)')
    
    plot3.legend(fontsize=size_label)
    if lang == "fr":
        plot3.set_title('Rééchantillonnage de '+term_header[index]+' ('+SYS0+')',fontsize=size_title)
        plot3.set_xlabel(r'$\tilde \xi$ (a.m.u.)',fontsize=size_label) #Coordonnée de réaction normalisée
        plot3.set_ylabel('$\Delta E_{mol}$ ($kJ.mol^{-1}$)',fontsize=size_label)
    else:    
        plot3.set_title('Energy fitting of '+term_header[index]+' ('+SYS0+')',fontsize=size_title)
        plot3.set_xlabel('Normalized reaction coordinate',fontsize=size_label)
        plot3.set_ylabel('Relative Energy ($kJ$ $mol^{-1}$)',fontsize=size_label)
    plot3.tick_params(axis='both',labelsize=size_label)
    plot3.set_xticks([i for i in range(0,len(crit_list_fit)+2)])
    plot3.grid(b=True,axis='x')
    """
    plot3.set_xlim(0,1)
    ylim0=2625.5*(min(term1_fit[index][:crit_list_fit[0]])-term1_fit[index][0])
    ylim1=2625.5*(min(term2_fit[index][:crit_list_fit[0]])-term2_fit[index][0])
    ylim2=2625.5*(min(term[index][:crit_list_fit[0]])-term[index][0])
    ylim3=2625.5*(max(term1_fit[index][:crit_list_fit[0]])-term1_fit[index][0])
    ylim4=2625.5*(max(term2_fit[index][:crit_list_fit[0]])-term2_fit[index][0])
    ylim5=2625.5*(max(term[index][:crit_list_fit[0]])-term[index][0])
    #print(ylim0,ylim1,ylim2,ylim3,ylim4,ylim5)
    ylim00=min(ylim0,ylim1,ylim2)
    ylim01=max(ylim3,ylim4,ylim5)
    plot3.set_ylim(ylim00,ylim01)
    """
    plot__.savefig(path+'/'+SYS0+'_energy_fitting_and_subtraction_'+term_header[index]+'_'+lang)
    
    print("Finished plotting the fitted, unfitted and substracted term profiles!")
    return
"""
###########################################################################################################
FUNCTION: set_axis_style
            set axis style for other graphs

INPUT: ax, labels
    ax      = subplot of a figure
    labels  = labels to be formatted

OUTPUT: ax (reformatted)

ERROR:
###########################################################################################################
"""

def set_axis_style(ax, labels):
    ax.get_xaxis().set_tick_params(direction='out',labelsize=13)
    ax.get_yaxis().set_tick_params(labelsize=13)
    ax.xaxis.set_ticks_position('bottom')
    ax.set_xticks(np.arange(1, len(labels) + 1))
    ax.set_xticklabels(labels)
    ax.set_xlim(0.25, len(labels) + 0.75)

"""
###########################################################################################################
FUNCTION: plot_violin
            violin plot of the REG analysis results (distribution of terms according to their R² coefficients)

INPUT: file_list = [], save = False, file_name = 'reg_violin.png', size_label = 16, lang = "en"
    file_list = []                  = list of the R² coefficient values for all terms over each segment
    save = False,                   = bool,
                                    True = save the plot as an image
                                    False = don't save the plot
    file_name = 'reg_violin.png'    = name of the saved file if save = True)
    size_label = 16                 = int, fontsize of the plot labels
    lang                            = "en" (default) displays graph text in English
                                    = "fr" displays graph text in French

OUTPUT: violin plot

ERROR:
###########################################################################################################
"""

def plot_violin(file_list = [], save = False, file_name = 'reg_violin.png', size_label = 16, lang = "en"):
    print("Starting to plot the global violin plot...")
    labels = [i+1 for i in range(len(file_list))]
    fig=plt.figure()
    if lang == "fr":
        fig.suptitle(r'Distribution de $R^2$ dans chaque segment')
    else:
        fig.suptitle(r'Distribution of $R^2$ in each segment')
    ax = fig.add_subplot(111)
    ax.violinplot(file_list, showmeans='True')
    set_axis_style(ax, labels)
    ax.set_xlabel('Segments',fontsize=size_label)
 
    if save == True:
        fig.savefig(file_name, dpi=300)
    print("Finished plotting the global violin plot!")
    return ax

"""
###########################################################################################################
FUNCTION: generate_data_vis
            plots of the REG analysis for each segment
INPUT: file, file_list, n_term, save = False, file_name = "data_vis.png", title = "REG data visualization", 
                      size_title = 18, size_label = 16, color = False, color_sign = False, sign_term = None, lang="en", verbose=False                     
    file                = full dataframe of R², REG coefficients and signs of terms
    file_list           = dataframe slice containing all R² coefficients of terms
    n_term              = number of significant terms to be considered
    save = False,       : bool,
                            True = save the plot as an image
                            False = don't save the plot
    file_name           = "data_vis.png" (name of the file to be saved)
    title               = "REG data visualization" (plot title)
    size_title          = int, fontsize of the plot title
    size_label          = int, fontsize of the plot labels
    color = False       : bool,
                            True = use colors depending on the fragment the terms belong to (for intermolecular processes) - hardcoded at the moment
                            False = use default colors for the stem plots
    color_sign = False  : bool,
                            True = use colors depending on the sign of the term to plot the stem plots
                            False = use default colors for the stem plots
    sign_term = None    = if the term signs have to be used to color the plot instead of the sign of their REG coefficient
    lang                = "en" (default) displays graph text in English
                        = "fr" displays graph text in French
    verbose             = False (default) : Boolean to ask the procedure to say everything (True) or be quiet (False)
    
OUTPUT: plots of the REG analysis for each segment
ERROR:
###########################################################################################################
"""

def generate_data_vis(file, file_list, n_term, save = False, file_name = "data_vis.png", title = "REG data visualization", 
                      size_title = 18, size_label = 16, color = False, color_sign = False, sign_term = None, lang="en", verbose=False):

    print("Starting to plot REG analysis graphs...")
    if verbose == True:
        print("Setting up the graphs...")
    (mu, sigma) = norm.fit(file['R2']) #AC change: R to R2.
    fig = plt.figure(constrained_layout=False, figsize=[10,10]) #default is 10,8 / change to 10,12
    fig.suptitle(title,fontsize=size_title,weight='bold')       
    gs1 = fig.add_gridspec(nrows=6, ncols=6, top=0.9, bottom=0.65, hspace=0, wspace = 0.40) #top=0.92, bottom=0.77, hspace=0, wspace = 0.40) #def top=0.9, bottom=0.65
    ax1 = fig.add_subplot(gs1[0:6,0:2]) # violin
    ax2 = fig.add_subplot(gs1[0:3,2:6]) # hist
    ax3 = fig.add_subplot(gs1[3:6,2:6]) #steam
    gs2 = fig.add_gridspec(nrows=6, ncols=6, top=0.55, bottom=0.05, hspace=0.05) #top=0.73, bottom=0.34, hspace=0.05) #def top=0.55, bottom=0.05
    ax4 = fig.add_subplot(gs2[:,:]) # data
    #if color == True:
        #gs3 = fig.add_gridspec(nrows=6, ncols=6, top=0.25, bottom=0.0, hspace=0.05)
        #ax5 = fig.add_subplot(gs3[:,:]) # data
    
##Violin_plot
    if verbose == True:
        print("Plotting the violin plot...")
    labels = [i+1 for i in range(len(file_list))]
    parts = ax1.violinplot(file_list, showmeans='True')
    set_axis_style(ax1, labels)
    for i in range(len(parts['bodies'])):
        if file.R2.equals(file_list[i]) == True: #AC change: R to R2.
            parts['bodies'][i].set_color('green')
            parts['bodies'][i].set_facecolor('green')
            parts['bodies'][i].set_edgecolor('green')
    ax1.set_xlabel('Segments',fontsize=size_label)
    ax1.grid(axis='y', alpha=0.75)
    if lang == "fr":
        ax1.set_title(r'Distribution des $R^2$',fontsize=size_title)
    else:
        ax1.set_title(r'Distribution of $R^2$',fontsize=size_title)

##histogram plot
    if verbose == True:
        print("Plotting the histogram plot...")
    ax2.hist(file['R2'], bins=[0.1*i for i in range(0,11,1)], 
             density=0 , color = 'green', edgecolor='black', linewidth=1.2, alpha=0.4) #AC change: R to R2.
    ax2.grid(axis='y', alpha=0.75)
    ax2.grid(axis='x', alpha=0.75)
    ax2.set_xlim(0,1)
    ax2.xaxis.set_ticks(np.arange(0, 1.1, 0.1))
    ax2.xaxis.set_ticks_position('top')
    ax2.xaxis.set_label_position('top')
    if lang == "fr":
        ax2.set_ylabel('Nombre de termes')
    else:
        ax2.set_ylabel('Count Number')
    ax2.set_xlabel(r'$R^2$')
    ax2.yaxis.set_ticks_position('right')
    ax2.axvline(mu, ls='--', color='r')

##stem plot1
    if verbose == True:
        print("Plotting the global stem plot...")
    ax3.set_xlim(0, 1)
    ax3.stem(file['R2'], file['REG'],'k', markerfmt=' ', use_line_collection ='True') #AC change: R to R2.
    ax3.grid(axis='x', alpha=0.75)
    if lang == "fr":
        ax3.set_ylabel('Valeur REG')
    else:
        ax3.set_ylabel('REG value')
    ax3.xaxis.set_ticks(np.arange(0, 1.1, 0.1))
    ax3.yaxis.set_ticks([])
    ax3.set_xlabel(r'$R^2$',fontsize=size_label)
    ax3.axvline(mu, ls='--', color='r')

    
##Stem plot2
    if verbose == True:
        print("Plotting the zoomed stem plot...")
        print("Adding positive REG values to the zoomed stem plot...")
    LIM=0.9 #display lower bound for the 2nd stem plot - could be an option of the function in the future
    pos = pd.DataFrame(columns=file.columns)
    neg = pd.DataFrame(columns=file.columns)    
    cond = file.REG < 0
    rows = file.loc[cond,:]
    neg = neg.append(rows, ignore_index=True)    
    cond = file.REG > 0
    rows = file.loc[cond,:]
    pos = pos.append(rows, ignore_index=True)
    if color_sign == False:        
        markerline, stemlines, baseline = ax4.stem(pos['R2'], pos['REG'], 'b', use_line_collection ='True') #AC change: R to R2.
        markerline.set_markerfacecolor('b')
        markerline.set_markeredgecolor('b')    
        markerline2, stemlines2, baseline2 = ax4.stem(neg['R2'], neg['REG'], 'r', use_line_collection ='True') #AC change: R to R2.
        markerline2.set_markerfacecolor('r')
        markerline2.set_markeredgecolor('r')    
    ax4.set_xlim(LIM, 1.0) #AC change: R to R2, mu to 0.8 then 0.9, max(file['R2']) to 1.
    
    bbox_props = dict(boxstyle="round,pad=0.1", fc="w", ec="k", lw=0, alpha=0.8)
    cond = pos.R2 > LIM #mu #AC change: R to R2. mu to 0.8
    temp = pos.loc[cond,:]
    n1 = pd.DataFrame()
    text_p = []
    if temp['REG'].dtype != np.dtype(object) :
        n = temp.nlargest(n_term, 'REG').reset_index(drop=True)
        n1 = temp.nlargest(1, 'REG').reset_index(drop=True) # Added by Aël Cador (4 April 2021)
        if n1.empty == False:
            r1 = n1['REG'][0] # Added by Aël Cador (4 April 2021)
    #Desactivate to get back to color mode. Modified by Aël Cador (24 June 2020)    
        
        global nmaxp,nminp
        if color == False and color_sign == False:
            text_p = [ax4.text(n['R2'][i], n['REG'][i], n['TERM'][i], bbox=bbox_props, fontweight='bold',fontsize=12) for i in range(len(n))]    #AC change: R to R2.
    
    # Added by Aël Cador (6 March 2020)
        if color == True:
            if verbose == True:
                print("Adding fragment-dependent color to the zoomed stem plot...")
            if n.empty == False:
                nmaxp=max(n['R2'])
                nminp=min(n['R2'])    
                #AC:color for the labels of the most important terms. Warning: this is hardcoded for a 29-atom system with 3 fragments of 7, 12 and 10 atoms respectively
                for i in range(len(n)):
                    atom_list=re.findall('\d+',n['TERM'][i])
                    if len(atom_list)==1:
                        if int(atom_list[0])<8: #Intra 1
                            text_p.append(ax4.text(n['R2'][i], n['REG'][i], n['TERM'][i], bbox=bbox_props, color="red", fontweight='bold',fontsize=12))
                        elif int(atom_list[0])>19: #Intra 3
                            text_p.append(ax4.text(n['R2'][i], n['REG'][i], n['TERM'][i], bbox=bbox_props, color="blue", fontweight='bold',fontsize=12))
                        else: #Intra 2
                            text_p.append(ax4.text(n['R2'][i], n['REG'][i], n['TERM'][i], bbox=bbox_props, color="orangered", fontweight='bold',fontsize=12))               
                    if len(atom_list)>1:
                        if int(atom_list[0])<8 and int(atom_list[1])<8: #Inter 1-1
                            text_p.append(ax4.text(n['R2'][i], n['REG'][i], n['TERM'][i], bbox=bbox_props, color="red", fontweight='bold',fontsize=12))
                        elif int(atom_list[0])<8 and int(atom_list[1])>19: # Inter 1-3
                            text_p.append(ax4.text(n['R2'][i], n['REG'][i], n['TERM'][i], bbox=bbox_props, color="purple", fontweight='bold',fontsize=12))                
                        elif int(atom_list[0])<8 and int(atom_list[1])>7 and int(atom_list[1])<20: #Inter 1-2
                            text_p.append(ax4.text(n['R2'][i], n['REG'][i], n['TERM'][i], bbox=bbox_props, color="orangered", fontweight='bold',fontsize=12))
                        elif int(atom_list[0])>7 and int(atom_list[0])<20 and int(atom_list[1])>19: #Inter 2-3
                            text_p.append(ax4.text(n['R2'][i], n['REG'][i], n['TERM'][i], bbox=bbox_props, color="purple", fontweight='bold',fontsize=12))                 
                        elif int(atom_list[0])>19 and int(atom_list[1])>19: #Inter 3-3
                            text_p.append(ax4.text(n['R2'][i], n['REG'][i], n['TERM'][i], bbox=bbox_props, color="blue", fontweight='bold',fontsize=12))
                        else: #Inter 2-2
                            text_p.append(ax4.text(n['R2'][i], n['REG'][i], n['TERM'][i], bbox=bbox_props, color="green", fontweight='bold',fontsize=12))
    
    # Added by Aël Cador (5 April 2021)    
        if color_sign == True: #for color according to sign - redo if time!
            if verbose == True:
                print("Adding sign-dependent color to the zoomed stem plot...")
            if n.empty == False:
                nmaxp=max(n['R2'])
                nminp=min(n['R2'])
                
                if max(file['SIGN']) == 2 or min(file['SIGN']) == -2:
                    colors=["gold","orange","red","blue","green"]
                    if lang == "fr":
                        color_labels=["Ambigu","Répulsif / Ambigu","Toujours répulsif","Toujours attractif","Attractif / Ambigu"]
                    else:
                        color_labels=["Ambiguous","Repulsive / Ambiguous","Always repulsive","Always attractive","Attractive / Ambiguous"]
                    #colorsl=["y","o","r","b","g"]
                    #AC:color for the lines and markers of all terms.
                    for i in range(-2,3):
                        cols = pd.DataFrame(columns=file.columns)    
                        cond = file.SIGN == float(i)
                        rows = file.loc[cond,:]
                        cols = cols.append(rows, ignore_index=True)
                        if cols.empty == False:
                            markerline, stemlines, baseline = ax4.stem(cols['R2'], cols['REG'], colors[i], use_line_collection ='True', label=color_labels[i], basefmt='k-')
                            markerline.set_markerfacecolor(colors[i])
                            markerline.set_markeredgecolor(colors[i])
                    baseline.set_markerfacecolor('k')
                    baseline.set_markeredgecolor('k')
                    
                    #AC:color for the labels of the most important REG-positive terms.
                    for i in range(len(n)):
                        #print(i,n['TERM'][i])
                        #print(n)
                        if 'Intra' in n['TERM'][i] :
                            text_p.append(ax4.text(n['R2'][i], n['REG'][i], n['TERM'][i], bbox=bbox_props, color=colors[int(n['SIGN'][i])], fontweight='bold',fontsize=12)) #sign_seg_intra[i]
                        elif n['TERM'][i] == 'Dispersion':
                            text_p.append(ax4.text(n['R2'][i], n['REG'][i], n['TERM'][i], bbox=bbox_props, color=colors[int(n['SIGN'][i])], fontweight='bold',fontsize=12))
                        else:
                            text_p.append(ax4.text(n['R2'][i], n['REG'][i], n['TERM'][i], bbox=bbox_props, color=colors[int(n['SIGN'][i])], fontweight='bold',fontsize=12)) #['index']
                
                else: #crappy criterion to see if you have two lists or one
                    colors=["green","red","blue"]
                    if lang == "fr":
                        color_labels=["Ambigu","Déstabilisant","Stabilisant"]
                        #color_labels=["Ambigu","Répulsif","Attractif"]
                    else:
                        color_labels=["Ambiguous","Destabilising","Stabilising"]
                    #AC:color for the lines and markers of all terms.
                    for i in range(-1,2):
                        cols = pd.DataFrame(columns=file.columns)    
                        cond = file.SIGN == float(i)
                        rows = file.loc[cond,:]
                        cols = cols.append(rows, ignore_index=True)
                        if cols.empty == False:
                            markerline, stemlines, baseline = ax4.stem(cols['R2'], cols['REG'], colors[i], use_line_collection ='True', label=color_labels[i], basefmt='k-')
                            markerline.set_markerfacecolor(colors[i])
                            markerline.set_markeredgecolor(colors[i])
                    baseline.set_markerfacecolor('k')
                    baseline.set_markeredgecolor('k')
                    
                    #AC:color for the labels of the most important REG-positive terms.
                    for i in range(len(n)):
                        #print(i,n['TERM'][i])
                        #print(n)
                        if 'Intra' in n['TERM'][i] :
                            text_p.append(ax4.text(n['R2'][i], n['REG'][i], n['TERM'][i], bbox=bbox_props, color=colors[int(n['SIGN'][i])], fontweight='bold',fontsize=12))
                        elif n['TERM'][i] == 'Dispersion':
                            text_p.append(ax4.text(n['R2'][i], n['REG'][i], n['TERM'][i], bbox=bbox_props, color=colors[int(n['SIGN'][i])], fontweight='bold',fontsize=12))
                        else:
                            text_p.append(ax4.text(n['R2'][i], n['REG'][i], n['TERM'][i], bbox=bbox_props, color=colors[int(n['SIGN'][i])], fontweight='bold',fontsize=12)) #['index']

        #print(text_p)
# End of modifications by Aël Cador (6 March 2020)
               
    #adjust_text(text_p, arrowprops=dict(arrowstyle='->', color='black'))#, lw=0.1)) 
    if verbose == True:
        print("Adding negative REG values to the zoomed stem plot...")
    cond = neg.R2 > LIM #mu  #AC change: R to R2. mu to 0.8
    temp = neg.loc[cond,:]
    if temp['REG'].dtype != np.dtype(object) :
        n = temp.nsmallest(n_term, 'REG').reset_index(drop=True)
        n2 = temp.nsmallest(1, 'REG').reset_index(drop=True) # Added by Aël Cador (4 April 2021)
        if n1.empty == False:
            r1 = n1['REG'][0]
            if n2.empty == False:
                r2 = n2['REG'][0] # Added by Aël Cador (4 April 2021)
                pad = 0.1*max(abs(r1),abs(r2)) # Added by Aël Cador (4 April 2021)
                ax4.set_ylim(r2-pad,r1+pad) # Added by Aël Cador (4 April 2021) this is just to have better visibilty of terms regardless of their value
    #Desactivate to get back to color mode. Modified by Aël Cador (24 June 2020)   
        text_n = []
        global nmaxn,nminn
        if color == False and color_sign == False:
            text_n = [ax4.text(n['R2'][i], n['REG'][i], n['TERM'][i], bbox=bbox_props,fontweight='bold',fontsize=12) for i in range(len(n))] #AC change: R to R2.
            adjust_text(text_p+text_n, arrowprops=dict(arrowstyle='->', color='black'))#, lw=0.1))
            
    # Added by Aël Cador (6 March 2020 / 24 June 2020)
        if color == True:
            if verbose == True:
                print("Adding fragment-dependent color to the zoomed stem plot...")
            nmaxn=max(n['R2'])
            nminn=min(n['R2'])
            #AC:color for the labels of the most important terms
            for i in range(len(n)):
                atom_list=re.findall('\d+',n['TERM'][i])
                #print(n['R2'][i], n['REG'][i])
                if len(atom_list)==1:#intra
                    if int(atom_list[0])<8: #Intra 1
                        text_n.append(ax4.text(n['R2'][i], n['REG'][i], n['TERM'][i], bbox=bbox_props, color="red", fontweight='bold',fontsize=12))
                    elif int(atom_list[0])>19: #Intra 3
                        text_n.append(ax4.text(n['R2'][i], n['REG'][i], n['TERM'][i], bbox=bbox_props, color="blue", fontweight='bold',fontsize=12))
                    else: #Intra 2
                        text_n.append(ax4.text(n['R2'][i], n['REG'][i], n['TERM'][i], bbox=bbox_props, color="orange", fontweight='bold',fontsize=12))               
                if len(atom_list)>1: #inter
                    if int(atom_list[0])<8 and int(atom_list[1])<8: #inter 1-1
                        text_n.append(ax4.text(n['R2'][i], n['REG'][i], n['TERM'][i], bbox=bbox_props, color="red", fontweight='bold',fontsize=12))
                    elif int(atom_list[0])<8 and int(atom_list[1])>19: #inter 1-3
                        text_n.append(ax4.text(n['R2'][i], n['REG'][i], n['TERM'][i], bbox=bbox_props, color="purple", fontweight='bold',fontsize=12))                
                    elif int(atom_list[0])<8 and int(atom_list[1])>7 and int(atom_list[1])<20: #inter 1-2
                        text_n.append(ax4.text(n['R2'][i], n['REG'][i], n['TERM'][i], bbox=bbox_props, color="orangered", fontweight='bold',fontsize=12))
                    elif int(atom_list[0])>7 and int(atom_list[0])<20 and int(atom_list[1])>19: #inter 2-3
                        text_n.append(ax4.text(n['R2'][i], n['REG'][i], n['TERM'][i], bbox=bbox_props, color="purple", fontweight='bold',fontsize=12))                 
                    elif int(atom_list[0])>19 and int(atom_list[1])>19: #inter  3-3
                        text_n.append(ax4.text(n['R2'][i], n['REG'][i], n['TERM'][i], bbox=bbox_props, color="blue", fontweight='bold',fontsize=12))
                    else: #inter 2-2
                        text_n.append(ax4.text(n['R2'][i], n['REG'][i], n['TERM'][i], bbox=bbox_props, color="green", fontweight='bold',fontsize=12))
    # End of modifications by Aël Cador (6 March 2020)

    # Added by Aël Cador (5 April 2021)
        if color_sign == True: #for color according to sign - redo if time
            if verbose == True:
                print("Adding sign-dependent color to the zoomed stem plot...")
            if n.empty == False:
                
                nmaxn=max(n['R2'])
                nminn=min(n['R2'])
                
                #AC:color for the labels of the most important REG-negative terms.
                if max(file['SIGN']) == 2 or min(file['SIGN']) == -2:
                    colors=["gold","orange","red","blue","green"]
                    for i in range(len(n)):
                        #print(i,n['TERM'][i])
                        #print(n)
                        if 'Intra' in n['TERM'][i] :
                            text_n.append(ax4.text(n['R2'][i], n['REG'][i], n['TERM'][i], bbox=bbox_props, color=colors[int(n['SIGN'][i])], fontweight='bold',fontsize=12))
                        elif n['TERM'][i] == 'Dispersion':
                            text_n.append(ax4.text(n['R2'][i], n['REG'][i], n['TERM'][i], bbox=bbox_props, color=colors[int(n['SIGN'][i])], fontweight='bold',fontsize=12))
                        else:
                            text_n.append(ax4.text(n['R2'][i], n['REG'][i], n['TERM'][i], bbox=bbox_props, color=colors[int(n['SIGN'][i])], fontweight='bold',fontsize=12)) #['index']
                else: #crappy criterion to see if you have two lists or one
                    colors=["green","red","blue"]
                    for i in range(len(n)):
                        #print(i,n['TERM'][i])
                        #print(n)
                        if 'Intra' in n['TERM'][i] :
                            text_p.append(ax4.text(n['R2'][i], n['REG'][i], n['TERM'][i], bbox=bbox_props, color=colors[int(n['SIGN'][i])], fontweight='bold',fontsize=12))
                        elif n['TERM'][i] == 'Dispersion':
                            text_p.append(ax4.text(n['R2'][i], n['REG'][i], n['TERM'][i], bbox=bbox_props, color=colors[int(n['SIGN'][i])], fontweight='bold',fontsize=12))
                        else:
                            text_p.append(ax4.text(n['R2'][i], n['REG'][i], n['TERM'][i], bbox=bbox_props, color=colors[int(n['SIGN'][i])], fontweight='bold',fontsize=12)) #['index']

    if temp['REG'].dtype != np.dtype(object) :
        adjust_text(text_p+text_n, arrowprops=dict(arrowstyle='->', color='black'))#, lw=0.1))
        #adjust_text(text_n, arrowprops=dict(arrowstyle='->', color='black'))#, lw=1.5))

    if lang=="fr":
        ax4.set_title('Contributions significatives',fontsize=size_title)
        ax4.set_ylabel('Valeur REG',fontsize=size_label)
    else:
        ax4.set_title('Relevant contributions',fontsize=size_title)
        ax4.set_ylabel('REG value',fontsize=size_label)
    ax4.set_xlabel(r'$R^2$',fontsize=size_label)
    ax4.tick_params(axis='both',labelsize=13)
    if color_sign == True:
        ax4.legend(loc='upper left',fontsize=size_label-2) #'best' #'upper left'
# Added by Aël Cador (6 March 2020 / 24 June 2020)    
# Added a legend to see which color is what + automatic positioning depending on where the REG terms are.   
    if color == True: # or color_sign == True:
        if nminn > 0.875: #mu+(max(file['R2'])-mu)/3:
            legx=0.025
            legy=0.05
        elif nminp > 0.875: #> mu+(max(file['R2'])-mu)/3:
            legx=0.025
            legy=0.75        
        elif nmaxn < 0.925: #< mu+2*(max(file['R2'])-mu)/3:
            legx=0.65
            legy=0.05         
        elif nmaxp < 0.925: #< mu+2*(max(file['R2'])-mu)/3:
            legx=0.65
            legy=0.75
        else:
            legx=0.025
            legy=0.05         
        if color == True:
            if lang == "fr":
                ax4.text(legx,legy,' Red: Fragment 1 only\n Orange: Interactions 1-2\n Green: Fragment 2 only\n Blue: Fragment 3 only\n Purple: Interactions 1-3 and 2-3',transform=ax4.transAxes)
            else:
                ax4.text(legx,legy,' Rouge : Fragment 1 uniquement\n Orange : Interactions 1-2\n Vert : Fragment 2 uniquement\n Bleu : Fragment 3 uniquement\n Violet : Interactions 1-3 et 2-3',transform=ax4.transAxes)
        
        #elif color_sign == True:
        #    ax4.text(legx,legy,' Red: Always destabilising\n Orange: Destabilising / ambiguous\n Golden: Ambiguous\n Green: Stabilising / ambiguous\n Blue: Always stabilising',transform=ax4.transAxes)
        
        """
        arr_6ts = mpimg.imread('../6-B3LYP-SI/6-B3LYP-SE_TS_borderless.png')
    
        imagebox = OffsetImage(arr_6ts, zoom=0.32)
    
        ab = AnnotationBbox(imagebox, (0.5,0.5))
    
        #ax5.add_artist(ab)
        """ 
# End of modifications by Aël Cador (6 March 2020) 
   
    if save == True:
        fig.savefig(file_name, dpi=300)
        
    plt.close(fig)
    
    print("Finished plotting all plots!")
    return

"""
###########################################################################################################
FUNCTION: plot_segment
            plot of the segments (energy and force only)
INPUT: coordinate, wfn_energy, critical_points, label=False, color=True, annotate=True, title='REG segments', 
                 y_label='Energy', x_label='Coordinate', size_title = 18, size_label = 16, save=False, 
                 file_name='segments.png', force=True, diff=False, lang="en", ylim=(-70,85), verbose=False
     coordinate                 = list of the control coordinate to use for the plots
     wfn_energy                 = list of the energy values on the path to use for the plots
     critical_points            = list of the critical points on the pathway
     label=False                = bool:
                                 True = add labels on the points on the segments
                                 False = don't add labels
     color=True                 = bool:
                                 True = color the segments
                                 False = leave segments white
     annotate=True              = bool:
                                 True = number the segments
                                 False = don't add the segment number
     title='REG segments'       = title of the plot
     y_label='Energy'           = Y axis title
     x_label='Coordinate'       = X axis title
     size_title                 = int, fontsize of the plot title
     size_label                 = int, fontsize of the plot labels
     save=False                 : bool,
                                True = save the plot as an image
                                False = don't save the plot 
     file_name='segments.png'   = file name
     force=True                 : bool,
                                True = plot the force
                                False = don't plot the force
     diff=False                 : bool,
                                True = this is an energy difference analysis (REGdiff) -> constrain the scales
                                False = normal REG analysis                 
    lang                        = "en" (default) displays graph text in English
                                = "fr" displays graph text in French
    ylim=(-70,85)               = Y axis limits (usually in kJ mol-1)
    verbose                     = False (default) : Boolean to ask the procedure to say everything (True) or be quiet (False)

OUTPUT: plot of the segments (energy and force only)
ERROR:
###########################################################################################################
"""

def plot_segment(coordinate, wfn_energy, critical_points, label=False, color=True, annotate=True, title='REG segments', 
                 y_label='Energy', x_label='Coordinate', size_title = 18, size_label = 16, save=False, 
                 file_name='segments.png', force=True, diff=False, lang="en", ylim=(-70,85), verbose=False):
    
    print("Starting to plot segmented energy profile...")
    color_list=['blue', 'darkgreen', 'darkred', 'yellow', 'cyan', 'magenta', 'grey', 'salmon',
                'seagreen', 'aquamarine', 'lightgreen', 'silver', 'lime', 'indigo', 'indianred']

# Modified by Aël Cador (6 March 2020)
# Slightly bigger figure to allow for better visibility    
    fig = plt.figure(figsize=(10,7)) #figure size in inches. 12,10 def is 8,5 #Change sometime to accomodate for X points automatically
# End of modifications by Aël Cador (6 March 2020)
    graph = plt.subplot(111)
    if verbose==True:
        print("Plotting the energy...")
    graph.plot(coordinate, wfn_energy, '--', color='black')
       
# Modified by Aël Cador (15 February 2020)
    if verbose==True:
        print("Setting up axes and labels...")     
    if coordinate[0]==0.0:
        graph.xaxis.set_major_locator(MultipleLocator(0.5*round(10*abs(coordinate[0]-coordinate[1])+1))) #For CH4D:0.2/0.1 #AC: changed the divisions of the x axis so that they don't overlap
        graph.xaxis.set_minor_locator(MultipleLocator(0.2*round(10*abs(coordinate[0]-coordinate[1])+1)))
    else:
        graph.xaxis.set_major_locator(MultipleLocator(2*round(abs(coordinate[0]-coordinate[1])+1))) #For 6 #AC: changed the divisions of the x axis so that they don't overlap
        graph.xaxis.set_minor_locator(MultipleLocator(round(abs(coordinate[0]-coordinate[1])+1)))
# End of modifications by Aël Cador (15 February 2020)      
    graph.set_title(title, fontsize=size_title, weight='bold')  
    graph.set_xlim([graph.get_xlim()[0], graph.get_xlim()[1]]) #wtf
    #if diff==True:
    graph.set_ylim(ylim[0],ylim[1]) # (-10,30)#(-70,85) #  #Added by Aël Cador (4 April 2021)
    if lang == "fr":
        graph.set_ylabel('$\Delta E_{mol}$ ($kJ.mol^{-1}$)', fontsize=size_label+2, color='green') #$\Delta E_{SCF}$ / Énergie relative
        if diff == True:
            graph.set_xlabel(r'$\tilde \xi$ (a.m.u.)',fontsize=size_label+2) # 'Coordonnée de réaction (a.m.u.)' #'$\Xi$' #Coordonnée de réaction '$\tilde{\Xi}$' #Coordonnée de réaction normalisée
        else:
            graph.set_xlabel(r'$\xi$ (a.m.u.)',fontsize=size_label+2)
    else:
        graph.set_ylabel(y_label, fontsize=size_label+2)
        graph.set_xlabel(x_label, fontsize=size_label+2)
    graph.tick_params(axis='both',labelsize=15)
    fig.savefig(file_name.replace(".png","_vanilla.png"), dpi=300)
    graph.plot(coordinate, wfn_energy, 'o--', color='black')
    if label == True:
        for i in range(len(wfn_energy)):
# Modified by Aël Cador (6 March 2020)            
            graph.annotate(str(i), [coordinate[i], wfn_energy[i]],xytext=(-5,8),textcoords=('offset pixels')) #AC: Added xytext and textcoords to set the label directly above the point
# End of modifications by Aël Cador (6 March 2020)       
        fig.savefig(file_name.replace(".png","_points.png"), dpi=300)

# Added by Aël Cador (6 March 2020)
# Added the derivative of the energy (the force)
    force_wfx=[]
    if force == True:
        if verbose==True:
            print("Computing the reaction force...")
        
        force_wfx=[-1*((wfn_energy[1]-wfn_energy[0])/(coordinate[1]-coordinate[0]))]
        for i in range(1,len(wfn_energy)-1): #don't forget the minus sign for the force!
            force_wfx.append(-0.5*(((wfn_energy[i+1]-wfn_energy[i])/(coordinate[i+1]-coordinate[i]))+((wfn_energy[i]-wfn_energy[i-1])/(coordinate[i]-coordinate[i-1]))))
        force_wfx.append(-1*((wfn_energy[-1]-wfn_energy[-2])/(coordinate[-1]-coordinate[-2])))
        force_wfx = np.array(force_wfx)
        if verbose==True:
            print("Reaction force values...",force_wfx)
        """ #force, old version
        force_wfx=[]
        for i in range(0,len(wfn_energy)-1):
            if coordinate[i]==0.0:
                force_wfx.append(0.5*(((wfn_energy[i+1]-wfn_energy[i])/(coordinate[i+1]-coordinate[i]))+((wfn_energy[i]-wfn_energy[i-1])/(coordinate[i]-coordinate[i-1]))))
            else:
                force_wfx.append((wfn_energy[i+1]-wfn_energy[i])/(coordinate[i+1]-coordinate[i]))
            
        force_wfx.append((wfn_energy[-1]-wfn_energy[-2])/(coordinate[-1]-coordinate[-2]))
        force_wfx = np.array(force_wfx)
        """
        if verbose==True:
            print("Adding the reaction force to the graph...")    
        axis2 = graph.twinx() #AC: Force axis
        axis2.set_ylim(graph.get_ylim()[0]/1.5, graph.get_ylim()[1]/1.5)
        axis2.plot(coordinate, force_wfx, '--', color='blue') #AC: Force plot
        if lang == "fr":
            axis2.set_ylabel(r'$F(\xi)$ (a.m.u.)', color = 'blue', fontsize=size_label+2) #AC: Force axis label / 'Force de réaction' / r'$F(\Xi)$'
        else:
            axis2.set_ylabel('Reaction force', color = 'blue', fontsize=size_label+2) #AC: Force axis label
        axis2.axhline(linestyle = 'dotted', color = 'blue')
        axis2.tick_params(axis='both',labelsize=15)
# End of modifications by Aël Cador (6 March 2020) 

    graph.axvline(coordinate[0], linestyle = 'dotted', color = 'black') #the first segment starts at the global minimum.
    graph.axvline(coordinate[-1], linestyle = 'dotted', color = 'black') #the first segment starts at the global minimum.
    
    for i in critical_points: #split the segments
        graph.axvline(coordinate[i], linestyle = 'dotted', color = 'black')
          
    if color == True: #Color_segments
        if verbose==True:
            print("Adding segment colors to the graph...")
        #graph.axvspan(graph.get_xlim()[0], graph.get_xlim()[1], facecolor = color_list[len(critical_points)], alpha = 0.2)
        start = 0
        for i in range(len(critical_points)):
            if coordinate[start] > coordinate[critical_points[i]]:
                graph.axvspan(coordinate[critical_points[i]], coordinate[start], facecolor='white', alpha =1)
                graph.axvspan(coordinate[critical_points[i]], coordinate[start], facecolor=color_list[i], alpha =0.2)
                fig.savefig(file_name.replace(".png","_"+str(i)+".png"), dpi=300)
            start=critical_points[i]            
        
        start = 0   
        for i in range(len(critical_points)):
            if coordinate[start] < coordinate[critical_points[i]]:
                graph.axvspan(coordinate[start], coordinate[critical_points[i]], facecolor='white', alpha =1)
                graph.axvspan(coordinate[start], coordinate[critical_points[i]], facecolor=color_list[i], alpha =0.2)
                fig.savefig(file_name.replace(".png","_"+str(i)+".png"), dpi=300)
            start=critical_points[i]
        graph.axvspan(coordinate[critical_points[-1]], coordinate[-1], facecolor=color_list[len(critical_points)], alpha =0.2)       

    if annotate== True: #Write segment number:
        if verbose==True:
            print("Adding the segment numbers to the graph...")
        text = []
        start=0
        y_pos = graph.get_ylim()[1] 
        x_pos = coordinate[start]
        text.append(graph.text(x_pos,  y_pos, '1',fontsize=14))
        
        for i in range(len(critical_points)):
            start=critical_points[i]
            x_pos = coordinate[start]
            text.append(graph.text(x_pos, y_pos,  str(i+2), fontsize=14))
            start=critical_points[i]
        adjust_text(text)
    
    if save == True:
        if force == True:
            fig.savefig(file_name.replace(".png","_force.png"), dpi=300)
        else:
            fig.savefig(file_name, dpi=300)
    print("Finished plotting segmented energy profile!")
    
    return
    

