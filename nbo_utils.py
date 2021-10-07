#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Apr 21 11:51:34 2021

@author: Aël
"""
"""
nbo_utils.py v0.0
A. Cador, P. L. A. Popelier 

Library with functions to get NBO properties values from both G09 and G16, NBO3 (Gaussian default) and NBO7.
Check for updates at github.com/ljduarte
For details about the method, please see XXXXXXX

Please, report bugs and issues to ael.cador@univ-rouen.fr
coded by A. Cador
"""

import copy

"""
###########################################################################################################

PROCEDURE LIST

###########################################################################################################
nbo_property_from_g09_file
          get NBO interactions from g09 files (NBO3)
nbo7_property_from_g16_file
          get NBO interactions from g16 files (NBO7) 
nbo_trim
          removes inconsistent NBO interactions (NBO3)           
nbo_trim_nbo7
          removes inconsistent NBO interactions (NBO7)
nbo_trimdiff
          removes NBO interactions which are not common to two pathways (NBO3)
nbo_trimdiff_nbo7
          removes NBO interactions which are not common to two pathways (NBO7)
sign_segm_nbo
        Returns the sign of each term on each segment:
            +1 if the term is always positive
             0 if the term is sometimes negative, sometimes positive
            -1 if the term is always negative
"""


"""
###########################################################################################################
Added by Aël Cador for REGNBO (April 2021)
FUNCTION: nbo_property_from_g09_file
          get NBO interactions from g09 files (NBO3)
          
          (It probably works with G16 outputs as long as it uses NBO3.)

INPUT: folders
    folders = path to Gaussian output files     
    verbose = False (default) : Boolean to ask the procedure to say everything (True) or be quiet (False)
    
OUTPUT: [contributions_list, nbo_properties_E2, nbo_properties_dE, nbo_properties_F]
    contributions_list  = list of E(2) term names
    nbo_properties_E2   = list of E(2) values (acceptor/donor orbital interaction energies)
    nbo_properties_dE   = list of dE values (acceptor/donor orbital energy difference)
    nbo_properties_F    = list of Fock matrix elements for each significant acceptor/donor interaction
    
ERROR:
    "File is empty or does not exist: " no file found or no NBO analysis found in file    
###########################################################################################################
""" 

def nbo_property_from_g09_file(folders, verbose=False): #, prop, atom_list):
    
    print("Starting reading NBO3 data from G09 files...")
    
    #INTERNAL VARIABLES:
    temp1 = [] #Temporary array (header)
    temp2 = [] #Temporary array (E(2))
    temp3 = [] #Temporary array (E(i)-E(j))
    temp4 = [] #Temporary array (F(i,j))
    contributions_list = [] #Output
    nbo_properties_E2 = [] #Output
    nbo_properties_dE = [] #Output
    nbo_properties_F = [] #Output
    
    for path in folders:
        file = open(path, "r")
        lines = file.readlines()
        file.close()
        for i in lines:
            #if 'NBO Version' not in lines:
            #    raise ValueError( "No NBO analysis found in: " + path )    
            if 'Second Order Perturbation Theory Analysis of Fock Matrix in NBO Basis' in i:
                start = lines.index(i)
            elif 'Natural Bond Orbitals (Summary):' in i:
                end = lines.index(i)
        if end >= len(lines) : #Checks the .int file.
            raise ValueError( "File is empty or does not exist: " + path ) 
        lines = [lines[i] for i in range(start+7, end)]
        for i in lines:
            if i != '\n' and 'within unit' not in i and 'from unit' not in i and 'None' not in i: #avoid lines with None above threshold
                temp1.append(''.join(i.split()[0:-3]))
                temp2.append(-4.184*float(i.split()[-3])) #kcal/mol to kJ/mol conversion
                temp3.append(float(i.split()[-2]))
                temp4.append(float(i.split()[-1]))
                if verbose==True:
                    print("Term name",path,''.join(i.split()[0:-3]))
                    print("E(2) interaction energy",-4.184*float(i.split()[-3]))
                    print("Orbital energy difference",float(i.split()[-2]))
                    print("Fock matrix element",float(i.split()[-1]))
        contributions_list.append(temp1)
        nbo_properties_E2.append(temp2)
        nbo_properties_dE.append(temp3)
        nbo_properties_F.append(temp4)
        #print(len(temp1))
        temp1 = [] #Temporary array (header)
        temp2 = [] #Temporary array (E(2))
        temp3 = [] #Temporary array (E(i)-E(j))
        temp4 = [] #Temporary array (F(i,j))
    
    print("Finished reading NBO3 data from G09 files!")
        
    return [contributions_list, nbo_properties_E2, nbo_properties_dE, nbo_properties_F]

"""
###########################################################################################################
Added by Aël Cador for REGNBO (April 2021)
FUNCTION: nbo7_property_from_g16_file
          get NBO interactions from g16 files (NBO7)

INPUT: folders
    folders = path to Gaussian output files     
    verbose = False (default) : Boolean to ask the procedure to say everything (True) or be quiet (False)
    
OUTPUT: [contributions_list, nbo_properties_E2, nbo_properties_dE, nbo_properties_F, 
         contributions_listE4,nbo_properties_SE4,nbo_properties_E4]
    contributions_list      = list of E(2) term names
    nbo_properties_E2       = list of E(2) values (acceptor/donor orbital interaction energies)
    nbo_properties_dE       = list of dE values (acceptor/donor orbital energy difference)
    nbo_properties_F        = list of Fock matrix elements for each significant acceptor/donor interaction
    contributions_listE4    = list of E(4) term names
    nbo_properties_SE4      = list of S values (occupied orbital overlap)
    nbo_properties_E4       = list of E(4) values (occupied orbital interaction energies)
    
ERROR:
    "File is empty or does not exist: " no file found or no NBO analysis found in file    
###########################################################################################################
"""

def nbo7_property_from_g16_file(folders, verbose=False): #, prop, atom_list):
    
    print("Starting reading NBO7 data from G16 files...")
    
    #INTERNAL VARIABLES:
    temp1 = [] #Temporary array (header)        (E(2))
    temp2 = [] #Temporary array (E(2))          (E(4))
    temp3 = [] #Temporary array (E(i)-E(j))     (E(4))
    temp4 = [] #Temporary array (F(i,j))        (E(4))
    temp5 = [] #Temporary array (header)        (E(4))
    temp6 = [] #Temporary array (PNLMO S(i,j))  (E(4))
    temp7 = [] #Temporary array (dE(i,j))       (E(4))
    contributions_listE2 = [] #Output
    nbo_properties_E2 = [] #Output
    nbo_properties_dE = [] #Output
    nbo_properties_F = [] #Output
    contributions_listE4 = [] #Output
    nbo_properties_SE4 = [] #Output
    nbo_properties_E4 = [] #Output
    
    for path in folders:
        if verbose==True:
            print("Reading file:",path)
        file = open(path, "r")
        lines = file.readlines()
        file.close()
        for i in lines:
            #if 'NBO Version' not in lines:
            #    raise ValueError( "No NBO analysis found in: " + path )    
            if 'SECOND ORDER PERTURBATION THEORY ANALYSIS OF FOCK MATRIX IN NBO BASIS' in i:
                start2 = lines.index(i)
            elif 'NATURAL BOND ORBITALS (Summary):' in i:
                end2 = lines.index(i)
            elif 'Pairwise steric exchange energies dE(i,j)' in i:
                start4 = lines.index(i)
            elif 'Total disjoint NLMO steric exchange energy from pairwise sum:' in i:
                end4 = lines.index(i)    
        if end2 >= len(lines) or end4 >= len(lines): #Checks the .log file.
            raise ValueError( "File is empty or does not exist: " + path )        
        
        lines2 = [lines[i] for i in range(start2+7, end2-2)]
        for i in lines2:
            if i != '\n' and 'within unit' not in i and 'between units' not in i and 'from unit' not in i and 'None' not in i: #avoid lines with None above threshold
                cont = i.split()[0:-3]
                for j in range(1,len(cont)):
                    if '.' in cont[j]:
                        cont.insert(j,'/')
                        break
                temp1.append(''.join(cont))
                temp2.append(-4.184*float(i.split()[-3])) #kcal/mol to kJ/mol conversion
                temp3.append(float(i.split()[-2]))
                temp4.append(float(i.split()[-1]))
                if verbose==True:
                    print("Term name",cont,temp1[-1])
                    print("Read from",path,''.join(i.split()[0:-3]))
                    print("E(2) interaction energy",-4.184*float(i.split()[-3]))
                    print("Orbital energy difference",float(i.split()[-2]))
                    print("Fock matrix element",float(i.split()[-1]))
        lines4 = [lines[i] for i in range(start4+7, end4-4)]
        for i in lines4:
            if i != '\n' and len(i.split())>1 and 'within unit' not in i and 'from unit' not in i and 'between units' not in i and 'None' not in i: #avoid lines with None above threshold
                cont = i.split()[0:-2]
                for j in range(1,len(cont)):
                    if '.' in cont[j]:
                        cont.insert(j,'/')
                        break
                temp5.append(''.join(cont))    
                temp6.append(float(i.split()[-2]))
                temp7.append(4.184*float(i.split()[-1])) #kcal/mol to kJ/mol conversion
                if verbose==True:
                    print("Term name",path,''.join(i.split()[0:-2]))
                    print("Orbital overlap",path,float(i.split()[-2]))
                    print("E(4) interaction energy",4.184*float(i.split()[-1]))
        contributions_listE2.append(temp1)
        nbo_properties_E2.append(temp2)
        nbo_properties_dE.append(temp3)
        nbo_properties_F.append(temp4)
        contributions_listE4.append(temp5)
        nbo_properties_SE4.append(temp6)
        nbo_properties_E4.append(temp7)
        #print(len(temp1))
        temp1 = [] #Temporary array (header)        (E(2))
        temp2 = [] #Temporary array (E(2))          (E(4))
        temp3 = [] #Temporary array (E(i)-E(j))     (E(4))
        temp4 = [] #Temporary array (F(i,j))        (E(4))
        temp5 = [] #Temporary array (header)        (E(4))
        temp6 = [] #Temporary array (PNLMO S(i,j))  (E(4))
        temp7 = [] #Temporary array (dE(i,j))       (E(4))
        
    print("Finished reading NBO7 data from G16 files!")
        
    return [contributions_listE2, nbo_properties_E2, nbo_properties_dE, nbo_properties_F,
            contributions_listE4, nbo_properties_SE4, nbo_properties_E4]

"""
###########################################################################################################
Added by Aël Cador for REGNBO (April 2021)
FUNCTION: nbo_trim
          removes inconsistent NBO interactions (NBO3)
          
          Basically it removes all NBO interactions that aren't found on the entire segment.
          If it finds the interaction for all files inside a given segment, it's kept, else it's removed.

INPUT: nbo_prop, crit_list
    nbo_prop    = list containing NBO3 interaction data lists, formatted as such and in this order:
        [contributions_list, nbo_properties_E2, nbo_properties_dE, nbo_properties_F]
        where each element is:
            contributions_list  = list of E(2) term names
            nbo_properties_E2   = list of E(2) values (acceptor/donor orbital interaction energies)
            nbo_properties_dE   = list of dE values (acceptor/donor orbital energy difference)
            nbo_properties_F    = list of Fock matrix elements for each significant acceptor/donor interaction
    crit_list   = critical point list of the total energy profile    
    verbose = False (default) : Boolean to ask the procedure to say everything (True) or be quiet (False)
    
OUTPUT: [contributions_list_final,nbo_properties_E2_final,nbo_properties_dE_final,nbo_properties_F_final]
    contributions_list_final    = list of E(2) term names
    nbo_properties_E2_final     = list of E(2) values (acceptor/donor orbital interaction energies)
    nbo_properties_dE_final     = list of dE values (acceptor/donor orbital energy difference)
    nbo_properties_F_final      = list of Fock matrix elements for each significant acceptor/donor interaction
    
ERROR:
    "The critical point indices do not match the length of the contributions list" : critical point out of range
    "The critical point list length does not match the length of the contributions list" : 
        too many critical points, or not enough data in contribution lists
###########################################################################################################
"""

def nbo_trim(nbo_prop, crit_list, verbose=False):
    
    print("Starting search for consistent terms on each segment...")
    
    contributions_list = nbo_prop[0]
    nbo_properties_E2 = nbo_prop[1]
    nbo_properties_dE = nbo_prop[2]
    nbo_properties_F = nbo_prop[3]
    contributions_list_final = [] #Output
    nbo_properties_E2_final = [] #Output
    nbo_properties_dE_final = [] #Output
    nbo_properties_F_final = [] #Output
    
    if crit_list[-1]>len(contributions_list):
        raise ValueError("The critical point indices do not match the length of the contributions list")        
    if len(crit_list)>len(contributions_list):
        raise ValueError("The critical point list length does not match the length of the contributions list")

    print("crit_list trim",crit_list)
    for h in range(len(crit_list)-1):
        indice_removal_0=[]
        if verbose==True:
            print("SEGMENT",h)
        start=crit_list[h]
        end=crit_list[h+1]+1
        if verbose==True:
            print("SEGMENT FROM POINTS",start,"TO",end)
        contributions_list_h = contributions_list[start:end].copy()
        nbo_properties_E2_h = nbo_properties_E2[start:end].copy()
        nbo_properties_dE_h = nbo_properties_dE[start:end].copy()
        nbo_properties_F_h = nbo_properties_F[start:end].copy()
        for i in range(len(contributions_list_h[0])):
            a=0
            for j in range(1,end-start):#len(contributions_list)):
                #for k in range(len(contributions_list[j])):
                    #if contributions_list[0][i] == contributions_list[j][k]:   
                
                if (contributions_list_h[0][i] not in contributions_list_h[j]):
                    if verbose==True:
                        print("Inconsistent contribution:",contributions_list[0][i])
                        print("at index",i)
                    indice_removal_0.append(i)    
                    a+=1
        indice_removal_0 = sorted(list(set(indice_removal_0)),reverse=True)
        if verbose==True:
            print("Contributions to remove",indice_removal_0,len(indice_removal_0),len(contributions_list[0]))
        if len(indice_removal_0) != 0: 
            for i in range(len(indice_removal_0)):
                if verbose==True:
                    print("Removing step:",i,indice_removal_0[i],len(contributions_list_h[0]))
                contributions_list_h[0].pop(indice_removal_0[i])
                nbo_properties_E2_h[0].pop(indice_removal_0[i])
                nbo_properties_dE_h[0].pop(indice_removal_0[i])
                nbo_properties_F_h[0].pop(indice_removal_0[i])
        for j in range(1,end-start):
            indice_removal_k = []
            for k in range(len(contributions_list_h[j])):
                if (contributions_list_h[j][k] not in contributions_list_h[0]):
                    indice_removal_k.append(k)
            if len(indice_removal_k) != 0:
                indice_removal_k = sorted(list(set(indice_removal_k)),reverse=True)
                if verbose==True:
                    print("LENGTH OF I.R.K.",len(indice_removal_k)) #,len(contributions_list[j]))
                    print(indice_removal_k,len(indice_removal_k),len(contributions_list[j]))
                for k in range(len(indice_removal_k)):
                    contributions_list_h[j].pop(indice_removal_k[k])
                    nbo_properties_E2_h[j].pop(indice_removal_k[k])
                    nbo_properties_dE_h[j].pop(indice_removal_k[k])
                    nbo_properties_F_h[j].pop(indice_removal_k[k])
                if verbose==True:
                    print("H",h,"LENGTH OF C.L.H.",len(contributions_list_h[j]))
        contributions_list_final.append(list(map(list, zip(*contributions_list_h))))
        if verbose==True:
            print(nbo_properties_E2_h[0])
        nbo_properties_E2_final.append(list(map(list, zip(*nbo_properties_E2_h)))) #list(map(list,zip(*l))) used to transpose without resorting to Numpy arrays
        nbo_properties_dE_final.append(list(map(list, zip(*nbo_properties_dE_h))))
        nbo_properties_F_final.append(list(map(list, zip(*nbo_properties_F_h))))
        if verbose==True:
            print(len(nbo_properties_E2_final[-1]),len(nbo_properties_E2_final[-1][0]))
        contributions_list_h = [] #Output
        nbo_properties_E2_h = [] #Output
        nbo_properties_dE_h = [] #Output
        nbo_properties_F_h = [] #Output
        indice_removal_k = []
    
    for i in range(len(contributions_list_final)):
        for j in range(len(contributions_list_final[i])):
            #print(contributions_list_final[i][j],contributions_list_final[i][j][0])
            contributions_list_final[i][j] = contributions_list_final[i][j][0]
    
    print("Search for consistent terms on each segment completed!")
    
    return [contributions_list_final, nbo_properties_E2_final, nbo_properties_dE_final, nbo_properties_F_final]

"""
###########################################################################################################
Added by Aël Cador for REGNBO (April 2021)
FUNCTION: nbo_trim_nbo7
          removes inconsistent NBO interactions (NBO7)
          
          Basically it removes all NBO interactions that aren't found on the entire segment.
          If it finds the interaction for all files inside a given segment, it's kept, else it's removed.

INPUT: nbo_prop, crit_list
    nbo_prop    = list containing NBO7 interaction data lists, formatted as such and in this order:
        [contributions_list, nbo_properties_E2, nbo_properties_dE, nbo_properties_F, 
         contributions_listE4, nbo_properties_SE4, nbo_properties_E4]
        where each element is:
            contributions_list      = list of E(2) term names
            nbo_properties_E2       = list of E(2) values (acceptor/donor orbital interaction energies)
            nbo_properties_dE       = list of dE values (acceptor/donor orbital energy difference)
            nbo_properties_F        = list of Fock matrix elements for each significant acceptor/donor interaction
            contributions_listE4    = list of E(4) term names
            nbo_properties_SE4      = list of S values (occupied orbital overlap)
            nbo_properties_E4       = list of E(4) values (occupied orbital interaction energies)
    crit_list   = critical point list of the total energy profile    
    verbose                     = False (default) : Boolean to ask the procedure to say everything (True) or be quiet (False)
    
OUTPUT: [contributions_listE2_final, nbo_properties_E2_final, nbo_properties_dE_final, nbo_properties_F_final,
            contributions_listE4_final, nbo_properties_SE4_final, nbo_properties_E4_final]
    contributions_list_final    = list of E(2) term names
    nbo_properties_E2_final     = list of E(2) values (acceptor/donor orbital interaction energies)
    nbo_properties_dE_final     = list of dE values (acceptor/donor orbital energy difference)
    nbo_properties_F_final      = list of Fock matrix elements for each significant acceptor/donor interaction
    contributions_listE4_final  = list of E(4) term names
    nbo_properties_SE4_final    = list of S values (occupied orbital overlap)
    nbo_properties_E4_final     = list of E(4) values (occupied orbital interaction energies)

ERROR: (duplicates for E(2) and E(4))
    "The critical point indices do not match the length of the E(X) contributions list" : 
        critical point out of range of the E(X) list
    "The critical point list length does not match the length of the E(X) contributions list" : 
        too many critical points, or not enough data in E(X) contribution lists
###########################################################################################################
"""

def nbo_trim_nbo7(nbo_prop, crit_list, verbose=False):
    
    print("Starting search for consistent terms on each segment...")
    
    contributions_listE2 = nbo_prop[0]
    nbo_properties_E2 = nbo_prop[1]
    nbo_properties_dE = nbo_prop[2]
    nbo_properties_F = nbo_prop[3]
    contributions_listE4 = nbo_prop[4]
    nbo_properties_SE4 = nbo_prop[5]
    nbo_properties_E4 = nbo_prop[6]
    
    contributions_listE2_final = [] #Output
    nbo_properties_E2_final = [] #Output
    nbo_properties_dE_final = [] #Output
    nbo_properties_F_final = [] #Output
    contributions_listE4_final = [] #Output
    nbo_properties_SE4_final = [] #Output
    nbo_properties_E4_final = [] #Output
    
    if crit_list[-1]>len(contributions_listE2):
        raise ValueError("The critical point indices do not match the length of the E(2) contributions list")        
    if len(crit_list)>len(contributions_listE2):
        raise ValueError("The critical point list length does not match the length of the E(2) contributions list")
    if crit_list[-1]>len(contributions_listE4):
        raise ValueError("The critical point indices do not match the length of the E(4) contributions list")        
    if len(crit_list)>len(contributions_listE4):
        raise ValueError("The critical point list length does not match the length of the E(4) contributions list")

    print("crit_list trim",crit_list)
    for h in range(len(crit_list)-1):
        indice_removal_0=[]
        if verbose==True:
            print("SEGMENT",h)
        start=crit_list[h]
        end=crit_list[h+1]+1
        if verbose==True:
            print("SEGMENT FROM POINTS",start,"TO",end)
        contributions_listE2_h = contributions_listE2[start:end].copy()
        nbo_properties_E2_h = nbo_properties_E2[start:end].copy()
        nbo_properties_dE_h = nbo_properties_dE[start:end].copy()
        nbo_properties_F_h = nbo_properties_F[start:end].copy()
        contributions_listE4_h = contributions_listE4[start:end].copy()
        nbo_properties_SE4_h = nbo_properties_SE4[start:end].copy()
        nbo_properties_E4_h = nbo_properties_E4[start:end].copy()
    
        for i in range(len(contributions_listE2_h[0])):
            a=0
            for j in range(1,end-start):#len(contributions_list)):
                #for k in range(len(contributions_list[j])):
                    #if contributions_list[0][i] == contributions_list[j][k]:   
                
                if (contributions_listE2_h[0][i] not in contributions_listE2_h[j]):
                    if verbose==True:
                        print("Inconsistent contribution:",contributions_listE2[0][i])
                        print("at index",i)
                    indice_removal_0.append(i)    
                    a+=1
        indice_removal_0 = sorted(list(set(indice_removal_0)),reverse=True)
        if verbose==True:
            print("Contributions to remove",indice_removal_0,len(indice_removal_0),len(contributions_listE2[0]))
        if len(indice_removal_0) != 0: 
            for i in range(len(indice_removal_0)):
                if verbose==True:
                    print("Removing step:",i,indice_removal_0[i],len(contributions_listE2_h[0]))
                contributions_listE2_h[0].pop(indice_removal_0[i])
                nbo_properties_E2_h[0].pop(indice_removal_0[i])
                nbo_properties_dE_h[0].pop(indice_removal_0[i])
                nbo_properties_F_h[0].pop(indice_removal_0[i])
        for j in range(1,end-start):
            indice_removal_k = []
            for k in range(len(contributions_listE2_h[j])):
                if (contributions_listE2_h[j][k] not in contributions_listE2_h[0]):
                    indice_removal_k.append(k)
            if len(indice_removal_k) != 0:
                indice_removal_k = sorted(list(set(indice_removal_k)),reverse=True)
                if verbose==True:
                    print("LENGTH OF I.R.K.",len(indice_removal_k)) #,len(contributions_list[j]))
                    print(indice_removal_k,len(indice_removal_k),len(contributions_listE2[j]))
                for k in range(len(indice_removal_k)):
                    contributions_listE2_h[j].pop(indice_removal_k[k])
                    nbo_properties_E2_h[j].pop(indice_removal_k[k])
                    nbo_properties_dE_h[j].pop(indice_removal_k[k])
                    nbo_properties_F_h[j].pop(indice_removal_k[k])
                if verbose==True:
                    print("H",h,"LENGTH OF C.L.H.",len(contributions_listE2_h[j]))
        contributions_listE2_final.append(list(map(list, zip(*contributions_listE2_h))))
        nbo_properties_E2_final.append(list(map(list, zip(*nbo_properties_E2_h)))) #list(map(list,zip(*l))) used to transpose without resorting to Numpy arrays
        if verbose==True:
            print("Final E2 first term",nbo_properties_E2_final[0])
        nbo_properties_dE_final.append(list(map(list, zip(*nbo_properties_dE_h))))
        nbo_properties_F_final.append(list(map(list, zip(*nbo_properties_F_h))))
        if verbose==True:
            print("Final last E2 term",len(nbo_properties_E2_final[-1]),len(nbo_properties_E2_final[-1][0]))
        contributions_listE2_h = [] #Output
        nbo_properties_E2_h = [] #Output
        nbo_properties_dE_h = [] #Output
        nbo_properties_F_h = [] #Output
        indice_removal_k = []
        
        indice_removal_0=[]
        for i in range(len(contributions_listE4_h[0])):
            a=0
            for j in range(1,end-start):
                if (contributions_listE4_h[0][i] not in contributions_listE4_h[j]):
                    if verbose==True:
                        print("Inconsistent contribution:",contributions_listE4[0][i])
                        print("at index",i)
                    indice_removal_0.append(i)    
                    a+=1
        indice_removal_0 = sorted(list(set(indice_removal_0)),reverse=True)
        if verbose==True:
            print("Contributions to remove",indice_removal_0,len(indice_removal_0),len(contributions_listE4[0]))
        if len(indice_removal_0) != 0: 
            for i in range(len(indice_removal_0)):
                if verbose==True:
                    print("Removing step:",i,indice_removal_0[i],len(contributions_listE4_h[0]))
                contributions_listE4_h[0].pop(indice_removal_0[i])
                nbo_properties_SE4_h[0].pop(indice_removal_0[i])
                nbo_properties_E4_h[0].pop(indice_removal_0[i])
        for j in range(1,end-start):
            indice_removal_k = []
            for k in range(len(contributions_listE4_h[j])):
                if (contributions_listE4_h[j][k] not in contributions_listE4_h[0]):
                    indice_removal_k.append(k)
            if len(indice_removal_k) != 0:
                indice_removal_k = sorted(list(set(indice_removal_k)),reverse=True)
                if verbose==True:
                    print("LENGTH OF I.R.K.",len(indice_removal_k)) #,len(contributions_list[j]))
                    print(indice_removal_k,len(indice_removal_k),len(contributions_listE4[j]))
                for k in range(len(indice_removal_k)):
                    contributions_listE4_h[j].pop(indice_removal_k[k])
                    nbo_properties_SE4_h[j].pop(indice_removal_k[k])
                    nbo_properties_E4_h[j].pop(indice_removal_k[k])
                if verbose==True:
                    print("H",h,"LENGTH OF C.L.H.",len(contributions_listE4_h[j]))
        contributions_listE4_final.append(list(map(list, zip(*contributions_listE4_h))))
        if verbose==True:
            print("First element of E4 terms:",contributions_listE4_final[0])
        nbo_properties_SE4_final.append(list(map(list, zip(*nbo_properties_SE4_h)))) #list(map(list,zip(*l))) used to transpose without resorting to Numpy arrays
        nbo_properties_E4_final.append(list(map(list, zip(*nbo_properties_E4_h))))
        if verbose==True and len(contributions_listE4_final) > 0:
            print("Length of first and last element of E4 terms",len(contributions_listE4_final[0]),len(contributions_listE4_final[-1]))
        elif verbose==True:
            print("Length of last element of E4 terms is zero")
        contributions_listE4_h = [] #Output
        nbo_properties_SE4_h = [] #Output
        nbo_properties_E4_h = [] #Output
        indice_removal_k = []
    
    for i in range(len(contributions_listE2_final)):
        if verbose==True:
            print(len(contributions_listE2_final[i]))
        for j in range(len(contributions_listE2_final[i])):
            #[j]),contributions_list_final[i][j][0])
            contributions_listE2_final[i][j] = contributions_listE2_final[i][j][0]
    for i in range(len(contributions_listE4_final)):
        if verbose==True:
            print(len(contributions_listE4_final[i]))
        for j in range(len(contributions_listE4_final[i])):
            if verbose==True:
                print(contributions_listE4_final[i][j], contributions_listE4_final[i][j][0])
            contributions_listE4_final[i][j] = contributions_listE4_final[i][j][0]
    
    print("Search for consistent terms on each segment completed!")
    
    return [contributions_listE2_final, nbo_properties_E2_final, nbo_properties_dE_final, nbo_properties_F_final,
            contributions_listE4_final, nbo_properties_SE4_final, nbo_properties_E4_final]

"""
###########################################################################################################
Added by Aël Cador for REGNBO (April 2021)
FUNCTION: nbo_trimdiff
          removes NBO interactions which are not common to two pathways (NBO3)
          
          Removes all NBO interactions that aren't found on two pathways.
          If it finds the interaction for all files on both pathways, it's kept, else it's removed.

INPUT: contributions_list1, nbo_properties_E21, nbo_properties_dE1, nbo_properties_F1,
        contributions_list2, nbo_properties_E22, nbo_properties_dE2, nbo_properties_F2, verbose=False
    contributions_list1  = list of E(2) term names (pathway 1)
    nbo_properties_E21   = list of E(2) values (acceptor/donor orbital interaction energies) (pathway 1)
    nbo_properties_dE1   = list of dE values (acceptor/donor orbital energy difference) (pathway 1)
    nbo_properties_F1    = list of Fock matrix elements for each significant acceptor/donor interaction (pathway 1)
    contributions_list2  = list of E(2) term names (pathway 2)
    nbo_properties_E22   = list of E(2) values (acceptor/donor orbital interaction energies) (pathway 2)
    nbo_properties_dE2   = list of dE values (acceptor/donor orbital energy difference) (pathway 2)
    nbo_properties_F2    = list of Fock matrix elements for each significant acceptor/donor interaction (pathway 2)    
    verbose                     = False (default) : Boolean to ask the procedure to say everything (True) or be quiet (False)
    
OUTPUT: [contributions_list1_t,nbo_properties_E21_t,nbo_properties_dE1_t,nbo_properties_F1_t,
                 contributions_list2_t,nbo_properties_E22_t,nbo_properties_dE2_t,nbo_properties_F2_t]
    contributions_list1_t  = list of E(2) term names (pathway 1)
    nbo_properties_E21_t   = list of E(2) values (acceptor/donor orbital interaction energies) (pathway 1)
    nbo_properties_dE1_t   = list of dE values (acceptor/donor orbital energy difference) (pathway 1)
    nbo_properties_F1_t    = list of Fock matrix elements for each significant acceptor/donor interaction (pathway 1)
    contributions_list2_t  = list of E(2) term names (pathway 2)
    nbo_properties_E22_t   = list of E(2) values (acceptor/donor orbital interaction energies) (pathway 2)
    nbo_properties_dE2_t   = list of dE values (acceptor/donor orbital energy difference) (pathway 2)
    nbo_properties_F2_t    = list of Fock matrix elements for each significant acceptor/donor interaction (pathway 2)

ERROR: (duplicates for E(2) and E(4))
    "The length of the contribution lists don't match" : the lengths of the two term name lists (paths 1/2) don't match
    "The length of the E(2) value lists don't match" : the lengths of the two E(2) value lists (paths 1/2) don't match
###########################################################################################################
"""

def nbo_trimdiff(contributions_list1, nbo_properties_E21, nbo_properties_dE1, nbo_properties_F1,
                 contributions_list2, nbo_properties_E22, nbo_properties_dE2, nbo_properties_F2, verbose=False):
    
    print("Starting search for terms common to both pathways...")
    
    contributions_list1_t = copy.deepcopy(contributions_list1) #without deepcopy, since the items inside the arrays are treated as objects,
    nbo_properties_E21_t = copy.deepcopy(nbo_properties_E21) #the standard copy doesn't work well (the objects stay linked)
    nbo_properties_dE1_t = copy.deepcopy(nbo_properties_dE1)
    nbo_properties_F1_t = copy.deepcopy(nbo_properties_F1)
    contributions_list2_t = copy.deepcopy(contributions_list2)
    nbo_properties_E22_t = copy.deepcopy(nbo_properties_E22)
    nbo_properties_dE2_t = copy.deepcopy(nbo_properties_dE2)
    nbo_properties_F2_t = copy.deepcopy(nbo_properties_F2)
    
    if len(contributions_list1) > len(contributions_list2):
        raise ValueError("The length of the contribution lists don't match")        
    if len(nbo_properties_E21) > len(nbo_properties_E22):
        raise ValueError("The length of the E(2) value lists don't match")
    if verbose==True:
        print("C1",len(contributions_list1))
        for i in range(len(contributions_list1)):
            print("C1i",i,len(contributions_list1[i]))
        print("E21",len(nbo_properties_E21))
        for i in range(len(nbo_properties_E21)):
            print("E21i",i,len(nbo_properties_E21[i]))
            for j in range(len(nbo_properties_E21[i])):
                print("E21ij",i,j,len(nbo_properties_E21[i][j]))
    
    for i in range(len(contributions_list1_t)):
        if verbose==True:
            print("C1.1",i,len(contributions_list1[i]))
            print("E21.1",i,len(nbo_properties_E21[i]))
            print("C2.1",i,len(contributions_list2[i]))
        indice_removal_0 = []
        for j in range(len(contributions_list1_t[i])):
            if contributions_list1_t[i][j] not in contributions_list2_t[i]:
                indice_removal_0.append(j)
                if verbose==True:
                    print(contributions_list1[i][j],j)
            indice_removal_0 = sorted(list(set(indice_removal_0)),reverse=True)    
        if len(indice_removal_0) != 0: 
            for k in range(len(indice_removal_0)):
                if verbose==True:
                    print(k,indice_removal_0[k],len(contributions_list1_t[i]))
                contributions_list1_t[i].pop(indice_removal_0[k])
                nbo_properties_E21_t[i].pop(indice_removal_0[k])
                nbo_properties_dE1_t[i].pop(indice_removal_0[k])
                nbo_properties_F1_t[i].pop(indice_removal_0[k])
        if verbose==True:
            print("C1.2",i,len(contributions_list1[i]))
            print("E21.2",i,len(nbo_properties_E21[i]))
            print("C2.2",i,len(contributions_list2[i]))
    for i in range(len(contributions_list2_t)):
        indice_removal_0 = []
        for j in range(len(contributions_list2_t[i])):
            if contributions_list2_t[i][j] not in contributions_list1_t[i]:
                indice_removal_0.append(j)
                if verbose==True:
                    print(contributions_list2[i][j],j)
            indice_removal_0 = sorted(list(set(indice_removal_0)),reverse=True)    
        if len(indice_removal_0) != 0: 
            for k in range(len(indice_removal_0)):
                if verbose==True:
                    print(k,indice_removal_0[k],len(contributions_list2_t[i]))
                contributions_list2_t[i].pop(indice_removal_0[k])
                nbo_properties_E22_t[i].pop(indice_removal_0[k])
                nbo_properties_dE2_t[i].pop(indice_removal_0[k])
                nbo_properties_F2_t[i].pop(indice_removal_0[k])
            if verbose==True:
                print("C1.3",len(contributions_list1[i]))
                print("E2.3",len(nbo_properties_E22[i]))
                print("C2.3",len(contributions_list2[i]))
    if verbose==True:
        print("C1",len(contributions_list1))
        for i in range(len(contributions_list1)):
            print("C1i",i,len(contributions_list1[i]))
        print("E21",len(nbo_properties_E21))
        for i in range(len(nbo_properties_E21)):
            print("E21i",i,len(nbo_properties_E21[i]))
            for j in range(len(nbo_properties_E21[i])):
                print("E21ij",i,j,len(nbo_properties_E21[i][j]))
        print("C2",len(contributions_list2))
        for i in range(len(contributions_list2)):
            print("C2i",i,len(contributions_list2[i]))
        print("E22",len(nbo_properties_E22))
        for i in range(len(nbo_properties_E22)):
            print("E22i",i,len(nbo_properties_E22[i]))
            for j in range(len(nbo_properties_E22[i])):
                print("E22ij",i,j,len(nbo_properties_E22[i][j]))
    
    print("Search for terms common to both pathways completed!")
    return [contributions_list1_t,nbo_properties_E21_t,nbo_properties_dE1_t,nbo_properties_F1_t,
                 contributions_list2_t,nbo_properties_E22_t,nbo_properties_dE2_t,nbo_properties_F2_t]

"""
###########################################################################################################
Added by Aël Cador for REGNBO (April 2021)
FUNCTION: nbo_trimdiff_nbo7
          removes NBO interactions which are not common to two pathways (NBO7)
          
          Removes all NBO interactions that aren't found on two pathways.
          If it finds the interaction for all files on both pathways, it's kept, else it's removed.

INPUT: contributions_listE21, nbo_properties_E21, nbo_properties_dE1, nbo_properties_F1,
        contributions_listE41, nbo_properties_SE41, nbo_properties_E41,
        contributions_listE22, nbo_properties_E22, nbo_properties_dE2, nbo_properties_F2,
        contributions_listE42, nbo_properties_SE42, nbo_properties_E42, verbose=False
    contributions_list1     = list of E(2) term names (pathway 1)
    nbo_properties_E21      = list of E(2) values (acceptor/donor orbital interaction energies) (pathway 1)
    nbo_properties_dE1      = list of dE values (acceptor/donor orbital energy difference) (pathway 1)
    nbo_properties_F1       = list of Fock matrix elements for each significant acceptor/donor interaction (pathway 1)
    contributions_listE41   = list of E(4) term names (pathway 1)
    nbo_properties_SE41     = list of S values (occupied orbital overlap) (pathway 1)
    nbo_properties_E41      = list of E(4) values (occupied orbital interaction energies) (pathway 1)
    contributions_list2     = list of E(2) term names (pathway 2)
    nbo_properties_E22      = list of E(2) values (acceptor/donor orbital interaction energies) (pathway 2)
    nbo_properties_dE2      = list of dE values (acceptor/donor orbital energy difference) (pathway 2)
    nbo_properties_F2       = list of Fock matrix elements for each significant acceptor/donor interaction (pathway 2)
    contributions_listE42   = list of E(4) term names (pathway 2)
    nbo_properties_SE42     = list of S values (occupied orbital overlap) (pathway 2)
    nbo_properties_E42      = list of E(4) values (occupied orbital interaction energies) (pathway 2)
    verbose                     = False (default) : Boolean to ask the procedure to say everything (True) or be quiet (False)

OUTPUT: [contributions_listE21_t, nbo_properties_E21_t, nbo_properties_dE1_t, nbo_properties_F1_t,
        contributions_listE41_t, nbo_properties_SE41_t, nbo_properties_E41_t,
        contributions_listE22_t, nbo_properties_E22_t, nbo_properties_dE2_t, nbo_properties_F2_t,
        contributions_listE42_t, nbo_properties_SE42_t, nbo_properties_E42_t]
    contributions_list1_t       = list of E(2) term names (pathway 1)
    nbo_properties_E21_t        = list of E(2) values (acceptor/donor orbital interaction energies) (pathway 1)
    nbo_properties_dE1_t        = list of dE values (acceptor/donor orbital energy difference) (pathway 1)
    nbo_properties_F1_t         = list of Fock matrix elements for each significant acceptor/donor interaction (pathway 1)
    contributions_listE41_t     = list of E(4) term names (pathway 1)
    nbo_properties_SE41_t       = list of S values (occupied orbital overlap) (pathway 1)
    nbo_properties_E41_t        = list of E(4) values (occupied orbital interaction energies) (pathway 1)
    contributions_list2_t       = list of E(2) term names (pathway 2)
    nbo_properties_E22_t        = list of E(2) values (acceptor/donor orbital interaction energies) (pathway 2)
    nbo_properties_dE2_t        = list of dE values (acceptor/donor orbital energy difference) (pathway 2)
    nbo_properties_F2_t         = list of Fock matrix elements for each significant acceptor/donor interaction (pathway 2)
    contributions_listE42_t     = list of E(4) term names (pathway 2)
    nbo_properties_SE42_t       = list of S values (occupied orbital overlap) (pathway 2)
    nbo_properties_E42_t        = list of E(4) values (occupied orbital interaction energies) (pathway 2)
    

ERROR: (duplicates for E(2) and E(4))
    "The length of the E(X) contribution lists don't match" : the lengths of the two term name lists (paths 1/2) don't match
    "The length of the E(X) value lists don't match" : the lengths of the two value lists (paths 1/2) don't match
###########################################################################################################
"""

def nbo_trimdiff_nbo7(contributions_listE21, nbo_properties_E21, nbo_properties_dE1, nbo_properties_F1,
                 contributions_listE41, nbo_properties_SE41, nbo_properties_E41,
                 contributions_listE22, nbo_properties_E22, nbo_properties_dE2, nbo_properties_F2,
                 contributions_listE42, nbo_properties_SE42, nbo_properties_E42, verbose=False):
    
    print("Starting search for terms common to both pathways...")
    contributions_listE21_t = copy.deepcopy(contributions_listE21) #without deepcopy, since the items inside the arrays are treated as objects,
    nbo_properties_E21_t = copy.deepcopy(nbo_properties_E21) #the standard copy doesn't work well (the objects stay linked)
    nbo_properties_dE1_t = copy.deepcopy(nbo_properties_dE1)
    nbo_properties_F1_t = copy.deepcopy(nbo_properties_F1)
    contributions_listE41_t = copy.deepcopy(contributions_listE41)
    nbo_properties_SE41_t = copy.deepcopy(nbo_properties_SE41)
    nbo_properties_E41_t = copy.deepcopy(nbo_properties_E41)
    
    contributions_listE22_t = copy.deepcopy(contributions_listE22)
    nbo_properties_E22_t = copy.deepcopy(nbo_properties_E22)
    nbo_properties_dE2_t = copy.deepcopy(nbo_properties_dE2)
    nbo_properties_F2_t = copy.deepcopy(nbo_properties_F2)
    contributions_listE42_t = copy.deepcopy(contributions_listE42)
    nbo_properties_SE42_t = copy.deepcopy(nbo_properties_SE42)
    nbo_properties_E42_t = copy.deepcopy(nbo_properties_E42)
    
    if len(contributions_listE21) > len(contributions_listE22):
        raise ValueError("The length of the E(2) contribution lists don't match")        
    if len(nbo_properties_E21) > len(nbo_properties_E22):
        raise ValueError("The length of the E(2) value lists don't match")
    if len(contributions_listE41) > len(contributions_listE42):
        raise ValueError("The length of the E(4) contribution lists don't match")        
    if len(nbo_properties_E41) > len(nbo_properties_E42):
        raise ValueError("The length of the E(4) value lists don't match")
    if verbose==True:
        print("C1",len(contributions_listE21_t))
        for i in range(len(contributions_listE21_t)):
            print("C1i",i,len(contributions_listE21_t[i]))
            
        print("E21",len(nbo_properties_E21_t))
        for i in range(len(nbo_properties_E21_t)):
            print("E21i",i,len(nbo_properties_E21_t[i]))
            for j in range(len(nbo_properties_E21_t[i])):
                print("E21ij",i,j,len(nbo_properties_E21_t[i][j]))
    
    #remove E(2) terms found only in pathway 1
    for i in range(len(contributions_listE21_t)):
        if verbose==True:
            print("C21.1",i,len(contributions_listE21_t[i]))
            print("E21.1",i,len(nbo_properties_E21_t[i]))
            print("C22.1",i,len(contributions_listE22_t[i]))
        indice_removal_0 = []
        for j in range(len(contributions_listE21_t[i])):
            if contributions_listE21_t[i][j] not in contributions_listE22_t[i]:
                indice_removal_0.append(j)
                if verbose==True:
                    print(contributions_listE21_t[i][j],j)
            indice_removal_0 = sorted(list(set(indice_removal_0)),reverse=True)    
        if len(indice_removal_0) != 0: 
            for k in range(len(indice_removal_0)):
                if verbose==True:
                    print(k,indice_removal_0[k],len(contributions_listE21_t[i]))
                contributions_listE21_t[i].pop(indice_removal_0[k])
                nbo_properties_E21_t[i].pop(indice_removal_0[k])
                nbo_properties_dE1_t[i].pop(indice_removal_0[k])
                nbo_properties_F1_t[i].pop(indice_removal_0[k])
        if verbose==True:
            print("C21.2",i,len(contributions_listE21_t[i]))
            print("E21.2",i,len(nbo_properties_E21_t[i]))
            print("C22.2",i,len(contributions_listE22_t[i]))
    
    #remove E(2) terms found only in pathway 2
    for i in range(len(contributions_listE22_t)):
        indice_removal_0 = []
        for j in range(len(contributions_listE22_t[i])):
            if contributions_listE22_t[i][j] not in contributions_listE21_t[i]:
                indice_removal_0.append(j)
                if verbose==True:
                    print(contributions_listE22_t[i][j],j)
            indice_removal_0 = sorted(list(set(indice_removal_0)),reverse=True)    
        if len(indice_removal_0) != 0: 
            for k in range(len(indice_removal_0)):
                #for j in range(len(contributions_list2[i])):
                if verbose==True:
                    print(k,indice_removal_0[k],len(contributions_listE22_t[i]))
                contributions_listE22_t[i].pop(indice_removal_0[k])
                nbo_properties_E22_t[i].pop(indice_removal_0[k])
                nbo_properties_dE2_t[i].pop(indice_removal_0[k])
                nbo_properties_F2_t[i].pop(indice_removal_0[k])
            if verbose==True:
                print("C21.3",len(contributions_listE21_t[i]))
                print("E22.3",len(nbo_properties_E22_t[i]))
                print("C22.3",len(contributions_listE22_t[i]))
    if verbose==True:
        print("C1",len(contributions_listE21_t))
        for i in range(len(contributions_listE21_t)):
            print("C1i",i,len(contributions_listE21_t[i]))
        print("E21",len(nbo_properties_E21_t))
        for i in range(len(nbo_properties_E21_t)):
            print("E21i",i,len(nbo_properties_E21_t[i]))
            for j in range(len(nbo_properties_E21_t[i])):
                print("E21ij",i,j,len(nbo_properties_E21_t[i][j]))
        print("C2",len(contributions_listE22_t))
        for i in range(len(contributions_listE22_t)):
            print("C2i",i,len(contributions_listE22_t[i]))
        print("E22",len(nbo_properties_E22_t))
        for i in range(len(nbo_properties_E22_t)):
            print("E22i",i,len(nbo_properties_E22_t[i]))
            for j in range(len(nbo_properties_E22_t[i])):
                print("E22ij",i,j,len(nbo_properties_E22_t[i][j]))
    
    #remove E(4) terms found only in pathway 1
    for i in range(len(contributions_listE41_t)):
        if verbose==True:
            print("C41.1",i,len(contributions_listE41_t[i]))
            print("E41.1",i,len(nbo_properties_E41_t[i]))
            print("C42.1",i,len(contributions_listE42_t[i]))
        indice_removal_0 = []
        for j in range(len(contributions_listE41_t[i])):
            if contributions_listE41_t[i][j] not in contributions_listE42_t[i]:
                indice_removal_0.append(j)
                if verbose==True:
                    print(contributions_listE41_t[i][j],j)
            indice_removal_0 = sorted(list(set(indice_removal_0)),reverse=True)    
        if len(indice_removal_0) != 0: 
            for k in range(len(indice_removal_0)):
                if verbose==True:
                    print(k,indice_removal_0[k],len(contributions_listE41_t[i]))
                contributions_listE41_t[i].pop(indice_removal_0[k])
                nbo_properties_SE41_t[i].pop(indice_removal_0[k])
                nbo_properties_E41_t[i].pop(indice_removal_0[k])
        if verbose==True:
            print("C41.2",i,len(contributions_listE41_t[i]))
            print("E41.2",i,len(nbo_properties_E41_t[i]))
            print("C42.2",i,len(contributions_listE42_t[i]))
            
    #remove E(4) terms found only in pathway 2
    for i in range(len(contributions_listE42_t)):
        if verbose==True:
            print("C41.1",i,len(contributions_listE41_t[i]))
            print("E42.1",i,len(nbo_properties_E42_t[i]))
            print("C42.1",i,len(contributions_listE42_t[i]))
        indice_removal_0 = []
        for j in range(len(contributions_listE42_t[i])):
            if contributions_listE42_t[i][j] not in contributions_listE41_t[i]:
                indice_removal_0.append(j)
                if verbose==True:
                    print(contributions_listE42_t[i][j],j)
            indice_removal_0 = sorted(list(set(indice_removal_0)),reverse=True)    
        if len(indice_removal_0) != 0: 
            for k in range(len(indice_removal_0)):
                if verbose==True:
                    print(k,indice_removal_0[k],len(contributions_listE41_t[i]))
                contributions_listE42_t[i].pop(indice_removal_0[k])
                nbo_properties_SE42_t[i].pop(indice_removal_0[k])
                nbo_properties_E42_t[i].pop(indice_removal_0[k])
        if verbose==True:
            print("C41.2",i,len(contributions_listE41_t[i]))
            print("E42.2",i,len(nbo_properties_E42_t[i]))
            print("C42.2",i,len(contributions_listE42_t[i]))
    print("Search for terms common to both pathways completed!.")
    
    return [contributions_listE21_t,nbo_properties_E21_t,nbo_properties_dE1_t,nbo_properties_F1_t,
                 contributions_listE41_t,nbo_properties_SE41_t,nbo_properties_E41_t,
                 contributions_listE22_t,nbo_properties_E22_t,nbo_properties_dE2_t,nbo_properties_F2_t,
                 contributions_listE42_t,nbo_properties_SE42_t,nbo_properties_E42_t]

"""
###########################################################################################################
Added by Aël Cador (03-06-21)
FUNCTION: sign_segm_nbo
        Returns the sign of each term on each segment:
            +1 if the term is always positive
             0 if the term is sometimes negative, sometimes positive
            -1 if the term is always negative
            
INPUT: term, crit_list, verbose=False
    term        = list of the term values over the complete path, organised as [[term1],[term2],...]
    crit_list   = list of the index of the critical points 
    verbose     = False (default) : Boolean to ask the procedure to say everything (True) or be quiet (False)

OUTPUT: sign_seg_term
    sign_seg_term = list of the signs of the terms over each segment, organised as:
    [[sign of term 1 over segment 1, sign of term 2 over segment 1,...] 
     [sign of term 1 over segment 2, sign of term 2 over segment 2,...]]
    
ERROR:
    "Invalid critical point index" : index of critical point out of range in the function array.
###########################################################################################################
"""

def sign_segm_nbo(term, crit_list, verbose=False):
    print("Starting NBO sign check...")
    len_term = 0
    NullTerm = 0
    for i in range(len(term)):
        if len(term[i])>0: 
            len_term += len(term[i][0])
            if verbose==True:
                print(len(term[i][0]),len_term)
        else:
            NullTerm+=1
    if crit_list[-1] >= len_term - len(crit_list) and NullTerm == 0: ## Checks if index of critical point is in range of function
        if verbose==True:
            print(crit_list[-1], len_term)
        raise ValueError("Invalid critical point index")
    split_term=term
    #split_intra1=reg.split_segm(iqa_intra1, crit_list1)
    sign_seg_term=[]
    for i in range(len(split_term)):
        sign_seg_term_i=[]
        for j in range(len(split_term[i])):
            sign_seg_term_j=[]
            for k in range(len(split_term[i][j])):
                if split_term[i][j][k]>0:
                    sign_seg_term_j.append(1)
                else:
                    sign_seg_term_j.append(-1)
            if sum(sign_seg_term_j) == len(split_term[i][j]):
                sign_seg_term_i.append(1)
            elif sum(sign_seg_term_j) == -len(split_term[i][j]):
                sign_seg_term_i.append(-1)
            else:
                sign_seg_term_i.append(0)
        sign_seg_term.append(sign_seg_term_i)
    
    if verbose==True:
        print(sign_seg_term)
    print("NBO sign check completed!")
    return sign_seg_term


                    