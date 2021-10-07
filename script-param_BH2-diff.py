# -*- coding: utf-8 -*-
"""
Created on Sun Apr  4 14:37:52 2021

@author: Ezan
"""

# IMPORT MODULES
"""
import reg 
import aimall_utils
import numpy as np
import numpy.polynomial.polynomial as poly
import matplotlib.pyplot as plt
import pandas as pd
import reg_vis as rv
"""
import script_REGdiff_common as rd
import script_REGdiff_RMSD_common as rdrm #something I've tried with the RMSD of the structures 
#(compared to the first one) as control coordinate instead of the reaction coordinate
#import script_REGNBO_common as nb

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
FORCE = True #If you want the force, or not
COLOR_TERMS = False # If you want your term names to be colored in the segment analysis, or not
COLOR_SIGN = True # If you want your term names to be colored according to the sign of the terms in the segment analysis, or not
SIZE_TITLE = 18 #fontsize of plot titles
SIZE_LABEL = 18 #fontsize of plot axis and tick labels
SAVE_FIG = True #if you want to save plot figures to the disk
ANNOTATE = True #if you want to write the segment number on the complete energy profile plot or not
DETAILED_ANALYSIS = True #if you want to have the plots of significant terms for each segment or not
LABELS = True #if you want to have point labels on the complete energy profile plot or not
n_terms = 5 #number of terms to be displayed on each segment graph (=number of significant terms)
LANG = "fr" #"en"

rd.REGdiff(SYS0,path1,F1,R1,path2,F2,R2,inv_1,inv_2,groups,POINTS,INFLEX,WRITE,AUTO,AUTO_fit,FORCE,COLOR_TERMS,COLOR_SIGN,SIZE_TITLE,
            SIZE_LABEL,SAVE_FIG,ANNOTATE,DETAILED_ANALYSIS,LABELS,n_terms,LANG)
#rdrm.REGdiff
#nb.REGNBO(SYS0,path1,F1,R1,path2,F2,R2,inv_1,inv_2,groups,POINTS,INFLEX,WRITE,AUTO,AUTO_fit,FORCE,COLOR_TERMS,COLOR_SIGN,SIZE_TITLE,
#            SIZE_LABEL,SAVE_FIG,ANNOTATE,DETAILED_ANALYSIS,LABELS,n_terms)