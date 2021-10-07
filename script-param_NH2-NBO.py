# -*- coding: utf-8 -*-
"""
Created on Sun May  2 23:34:09 2021

@author: Ezan
"""

# IMPORT MODULES
import script_REGdiff_common as rd
import script_REGdiff_RMSD_common as rdrm
import script_REGNBO_common as nb
import script_REGNBO_common_v2 as nb2
import script_REGNBO_common_v3 as nb3

##############################    OPTIONS    ##################################
#SYSTEM NAME
SYS0 = 'NH2'

#DEFINE PATHS AND FILES FOR PATH 1: 
path1 = SYS0+"-out" #Path to the folder with .wfn 
F1=[i for i in range(0,67,1)] #the Forward SP/AIM file numbers ; don't forget to add 5 to the last one (e.g. F0-F70 -> 0,75)
R1=[i for i in range(1,133,1)] #the Reverse SP/AIM file numbers ; don't forget to add 5 to the last one

#DEFINE PATHS AND FILES FOR PATH 2:
path2 = SYS0+"-in" #Path to the folder with .wfn 
F2=[i for i in range(0,74,1)] #the Forward SP/AIM file numbers ; don't forget to add 5 to the last one
R2=[i for i in range(1,125,1)] #the Reverse SP/AIM file numbers ; don't forget to add 5 to the last one

#CALCULATION OPTIONS
inv_1=True #True if the file numbering is inverted for path 1 (Forward->Reactant and Reverse->Product)
inv_2=True #True if the file numbering is inverted for path 2 (Forward->Reactant and Reverse->Product)
groups=[['n10','h11','h12']]
POINTS = 4 #minimal number of points on each segment
INFLEX = True #if the inflexion points have to be considered as critical points or not

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
LABELS = False  # if you want to have point labels on the complete energy profile plot or not
# number of terms to be displayed on each segment graph (=number of significant terms)
n_terms = 5
lang="en"

"""
nb.REGNBO(SYS0, path1, F1, R1, path2, F2, R2, inv_1, inv_2, groups, POINTS, INFLEX, WRITE, AUTO, AUTO_fit, FORCE, COLOR_TERMS, COLOR_SIGN, SIZE_TITLE,
            SIZE_LABEL, SAVE_FIG, ANNOTATE, DETAILED_ANALYSIS, LABELS, n_terms)
nb2.REGNBO(SYS0, path1, F1, R1, path2, F2, R2, inv_1, inv_2, groups, POINTS, INFLEX, WRITE, AUTO, AUTO_fit, FORCE, COLOR_TERMS, COLOR_SIGN, SIZE_TITLE,
            SIZE_LABEL, SAVE_FIG, ANNOTATE, DETAILED_ANALYSIS, LABELS, n_terms)
"""
nb3.REGNBO(SYS0, path1, F1, R1, path2, F2, R2, inv_1, inv_2, groups, POINTS, INFLEX, WRITE, AUTO, AUTO_fit, FORCE, COLOR_TERMS, COLOR_SIGN, SIZE_TITLE,
            SIZE_LABEL, SAVE_FIG, ANNOTATE, DETAILED_ANALYSIS, LABELS, n_terms, lang)