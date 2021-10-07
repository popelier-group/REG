#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
reg.py v0.0
L. J. Duarte, A. Cador, P. L. A. Popelier 

Library with function to perform the IQA and NBO REG analysis. 
Check for updates at github.com/ljduarte
For details about the method, please see XXXXXXX

Please, report bugs and issues to leo.j.duarte@hotmail.com.br and ael.cador@univ-rouen.fr
coded by L. J. Duarte
modified and expanded by A. Cador
"""

import numpy as np
import numpy.polynomial.polynomial as poly
#from sklearn.linear_model import LinearRegression
import pandas as pd
#import matplotlib.pyplot as plt
import warnings #used as a patch to hide numpy.polyfit error messages in reg.fit

"""
###########################################################################################################

PROCEDURE LIST

###########################################################################################################

regression
          Perform a linear regression between A and B (B = slope*A + intercept)
find_critical
          Takes the X and Y values of a function and return the critical points
          that fit the desirable criteria.
crit_trim
          Removes critical points from critical point lists if the lengths of the two lists differ
          (for REGdiff, this is useful to keep the same number of segments).
          The choice of the critical points is left to the user.
normalize
          Takes the control coordinate and normalizes it according to the critical points:
          for each segment, the control coordinate is set to N and N+1 at its limits and 
          is normalized (e.g. the segment [-6,-4,-3,-2] is normalized to [0.0, 0.5, 0.75, 1.0])
norm_diff
          Outputs the critical points of the original profile corresponding to the critical points of the difference profile.         
          In the case where the difference profile is non-monotonous, useful to change the critical points.
fit
          Fits a 4th-order polynomial to the term values with a resampling of the points 
          (the number of points equals the length of the longest segment), with corrections 
          of higher orders up to half the number of points.          
split_segm
          Takes the A array and divide it into N arrays according with the number of 
          critical points. Each array corresponds to a segment of the REG analysis.
sign_segm
          Returns the sign of each term on each segment:
            +1 if the term is always positive
             0 if the term is sometimes negative, sometimes positive
            -1 if the term is always negative
reg
          Perform the REG analysis over all contributions inside "terms" array
regnbo
          Perform the REG-NBO analysis over all contributions inside "terms" array
reg_zeroes
          removes terms with null REG coefficients (constant terms)
integration_error
          calculates integration error for each PES point
group_intra
          group IQA intra-atomic terms into the user-defined groups
group_inter
          group IQA interatomic terms into the user-defined groups
group_vcvx
          group Vc and Vx IQA interatomic terms
csv
          generates .csv and .xlsx files containing the results of the REG analysis
csvnbo
          generates .csv and .xlsx files containing the results of the REG-NBO3 analysis (Gaussian default)
csvnbo7
          generates .csv and .xlsx files containing the results of the REG-NBO7 analysis (most recent NBO version)
excel
          generates .xlsx files containing the rawa data of the REG analysis
excelnbo
          generates .xlsx files containing the raw data of the REG-NBO3 analysis
excelnbo7
          generates .xlsx files containing the raw data of the REG-NBO7 analysis
           
###########################################################################################################
"""

"""
###########################################################################################################
FUNCTION: regression
          Perform a linear regression between A and B (B = slope*A + intercept)

INPUT: A, B and mode
    A    = X values
    B    = Y values
    mode = None (default) : for regular linear regression
           "norm"         : use normalized values of A and B
           "std"          : use standardized values for A and B
           
OUTPUT: [slope, intercept, person]
    slope       : angular coefficient of the linear regression
    intercept   : linear coefficient of the linear regression
    pearson     : Pearson correlation coefficient
    
ERROR:
    "Arrays must have the same size" : A and B have different sizes.
    "Mode not recognized. Use None, 'norm' or, 'std'" : invalide value for 'mode'
###########################################################################################################
"""

def regression(A, B, mode=None):
    #ERRORS:  
    if len(A) != len(B): #Checks if A and B have the same size
        return ValueError("Arrays must have the same size") 
    if mode != None and mode != "norm" and mode != "std": #Checks if mode in valid
        return ValueError("Mode not recognized. Use None, 'norm' or, 'std'")
    
    #MODIFY A AND B IF REQUESTED   
    if mode == "norm":
        avgA = sum(A)/len(A)
        avgB = sum(B)/len(B)
        A = [a-avgA for a in A ]
        B = [b-avgB for b in B ]
    if mode == "std":
        avgA = sum(A)/len(A)
        sqA = [(a-avgA)**2 for a in A]
        stdevA = (sum(sqA)/(len(A)-1))**(0.5)
        A = [(a-avgA)/stdevA for a in A]
        avgB = sum(B)/len(B)
        sqB = [(b-avgB)**2 for b in B]
        stdevB = (sum(sqB)/(len(B)-1))**(0.5)
        B = [(b-avgB)/stdevB for b in B]
        
    #PERFORM THE LINEAR REGRESSION
    sum_AB = sum([A[i]*B[i] for i in range(len(A))])
    sum_A_sum_B = sum(A)*sum(B)
    sqsum_A = sum(A)**2
    sum_sqA = sum([a**2 for a in A])
    sum_sqB = sum([b**2 for b in B])
    slope = (len(A)*sum_AB - sum_A_sum_B)/(len(A)*sum_sqA - sqsum_A) 
    intercept = (sum(B)*sum_sqA - sum(A)*sum_AB)/(len(A)*sum_sqA-sqsum_A)
    pearson = (len(A)*sum_AB - sum(A)*sum(B))/((len(A)*sum_sqA-sum(A)**2)*(len(A)*sum_sqB-sum(B)**2))**(0.5)
    
    return slope, intercept, pearson

#RT_inv = RT_inv.reshape((-1, 1))
#model = LinearRegression().fit(RT_inv,k6_log)
#E6 = -model.coef_
#A6 = numpy.exp(model.intercept_)


"""
###########################################################################################################
FUNCTION: find_critical
          Takes the X and Y values of a function and return the critical points
          that fit the desirable criteria.
          
INPUT: Y, X, min_points=3, use_inflex=False, verbose=False
    Y           = ordinate values of a function
    X           = abcissa values of a function
    min_point   = 3 (default) : minimum amount of points between two critical points
    use_inflex  = False (default) : Do not search for inflexion points (second derivative = 0)
                  True            : Search for inflexion points (second derivative = 0)
    verbose     = False (default) : Boolean to ask the procedure to say everything (True) or be quiet (False)
           
OUTPUT: critical_point_list
    critical_point_list = Array containing the index of critical points (zero indexed)    

ERROR:
    "Arrays must have the same size" : X and Y have different sizes.
    "Invalid value. Use True or False" : invalid value for use_inflex
    "Too many points between two critical points" : min_points must be lower than the X array length

A. Cador modifications: derivative formula (from the slope to the next point, to the average of the slopes on each side)
###########################################################################################################
"""

def find_critical(Y, X, min_points=3, use_inflex=False, verbose=False):
    print("Starting critical point search...")
    #ERRORS:  
    if len(Y) != len(X): ## Checks if X and Y have the same size
        raise ValueError("Arrays must have the same size") 
    if use_inflex != True and use_inflex != False: ## Checks if use_inflex is valid
        raise ValueError("Invalid value. Use True or False")
    if min_points >= len(X):
        raise ValueError("Too many points between two critical points") ## Check if min_points is valid.
   
    #INTERNAL VARIABLES:
    critical_point = [] ## temporary array
    critical_point_list = [] ## final array
    
    #Y_prime = [(Y[i+1]-Y[i])/(X[i+1]-X[i]) for i in range(len(Y)-1)] ## First order derivative
    Y_prime = [(Y[1]-Y[0])/(X[1]-X[0])]
    for i in range (1,len(Y)-2):
        Y_prime.append(0.5*(((Y[i+1]-Y[i])/(X[i+1]-X[i]))+((Y[i]-Y[i-1])/(X[i]-X[i-1])))) ## First order derivative 
    Y_prime.append(((Y[-1]-Y[-2])/(X[-1]-X[-2])))

    #FIND CRITICAL POINTS (FIRST DERIVATIVE = 0)    
    if use_inflex == False:
        if verbose==True:
            print("Looking for maxima and minima...")
        for i in range (1, len(Y)-1):
            if Y[i] < Y[i-1] and Y[i] < Y[i+1]:               
                critical_point.append(i)
            if Y[i] > Y[i-1] and Y[i] > Y[i+1]:               
                critical_point.append(i)

        for i in range(1, len(Y_prime)-1):
           if Y_prime[i] == 0:
               critical_point.append(i+1)
           if Y_prime[i] > 0 and Y_prime[i+1] < 0:
               critical_point.append(i+1)
           if Y_prime[i] < 0 and Y_prime[i+1] > 0:
               critical_point.append(i+1)                        
    
    #FIND CRITICAL POINTS (FIRST DERIVATIVE = 0, SECOND DERIVATIVE = 0)   
    if use_inflex == True:
        if verbose==True:
            print("Looking for maxima, inflexion points and minima...")
        Y_2prime = [(Y_prime[i+1]-Y_prime[i])/(X[i+2]-X[i+1]) for i in range(len(Y_prime)-1)] ## First order derivative  

        for i in range (1, len(Y)-1):
            if Y[i] < Y[i-1] and Y[i] < Y[i+1]:               
                critical_point.append(i)
            if Y[i] > Y[i-1] and Y[i] > Y[i+1]:               
                critical_point.append(i)

        for i in range(0, len(Y_prime)-1):
            if Y_prime[i] == 0:
                critical_point.append(i+1)
            if Y_prime[i] > 0 and Y_prime[i+1] < 0:
                critical_point.append(i+1)
            if Y_prime[i] < 0 and Y_prime[i+1] > 0:
                critical_point.append(i+1)
                
        for i in range(0, len(Y_2prime)-1):
            if Y_2prime[i] == 0:
                critical_point.append(i+1)
            if Y_2prime[i] < 0 and Y_2prime[i+1] > 0:
                critical_point.append(i+1)               
            if Y_2prime[i] > 0 and Y_2prime[i+1] < 0:
                critical_point.append(i+1)
                             
    critical_point = list(set(critical_point)) ## Remove duplicates and sort in crescent order
    critical_point.sort()
    
    if len(critical_point)==0: ## checks if critical points were found. 
        raise ValueError("No critical point found") 

    #CHECKS MIN_POINTS CRITERIA
    while critical_point[0] < min_points-1 : critical_point.remove(critical_point[0])
    critical_point_list.append(critical_point[0])
    for i in range(1, len(critical_point)):
        if critical_point[i] - critical_point_list[-1] + 1 >= min_points:
           critical_point_list.append(critical_point[i])
    
    if verbose==True:
        print("Critical points found:",critical_point_list)
    print("Critical point search succesfully completed.")
    return critical_point_list

"""
###########################################################################################################
FUNCTION: crit_trim
        Removes critical points from critical point lists if the lengths of the two lists differ
        (for REGdiff, this is useful to keep the same number of segments).
        The choice of the critical points is left to the user.
          
INPUT: crit_list1, crit_list2, SYS1, SYS2
    crit_list1  = list of the critical points of the 1st energy profile
    crit_list2  = list of the critical points of the 2nd energy profile
    SYS1        = name of the first system
    SYS2        = name of the second system
    
OUTPUT: crit_list1,crit_list2
    crit_list1 = modified (if necessary) list of the critical points of the 1st energy profile
    crit_list2 = modified (if necessary) list of the critical points of the 2nd energy profile

ERROR:
    Add error check if the point entered by the user is not in the critical point list
    
Added by Aël Cador for REGdiff (04-04-21)
###########################################################################################################
"""

def crit_trim(crit_list1, crit_list2, SYS1, SYS2):
    if len(crit_list1) != len(crit_list2):
        crit_diff=abs(len(crit_list1)-len(crit_list2))
        print("The number of critical points differs on each pathway:")
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
        return [crit_list1,crit_list2]
    else:
        print("The two pathways have the same number of critical points.")
        return [crit_list1,crit_list2]

"""
###########################################################################################################
FUNCTION: normalize
          Takes the control coordinate and normalizes it according to the critical points:
          for each segment, the control coordinate is set to N and N+1 at its limits and 
          is normalized (e.g. the segment [-6,-4,-3,-2] is normalized to [0.0, 0.5, 0.75, 1.0])
          
INPUT: coord, critical_point_list, verbose=False
    coord               = list of the control coordinate
    critical_point_list = list of the critical points of the energy profile
    verbose             = False (default) : Boolean to ask the procedure to say everything (True) or be quiet (False)
    
OUTPUT: normalized_coord
    normalized_coord = the normalized control coordinate

ERROR:
    "The critical point indices do not match the length of the control coordinate list": 
        the indices of the critical points are out of range of the control coordinate list
    "The critical point list length does not match the length of the control coordinate list": 
        there are more critical points than points in the control coordinate list

Added by Aël Cador for REGdiff (03-03-21)
###########################################################################################################
"""

def normalize(coord, critical_point_list, verbose=False):
    print("Starting normalization...")
    if critical_point_list[-1]>len(coord):
        raise ValueError("The critical point indices do not match the length of the control coordinate list")        
    if len(critical_point_list)>len(coord):
        raise ValueError("The critical point list length does not match the length of the control coordinate list")
        
    critical_list = critical_point_list.copy() #copy to be able to edit it
    normalized_coord = []
    start=0 #define the limits of the 1st segment
    crit=critical_list[0]
    for i in range(len(coord)):
        #print(start,crit)
        if i <= crit:
            renorm=((coord[i]-coord[start])/(coord[crit]-coord[start])) + len(critical_point_list) - len(critical_list)
            #print(renorm)
            #the first part ensures the coordinate is in [0,1] depending on its position on the segment
            #the second one adds 1 for every segment because the elements in critical_list are gradually removed
            normalized_coord.append(renorm)
        if i==crit:
            start=crit #the end of segment I is the start of segment I+1
            if len(critical_list)>0:
                critical_list.pop(0) #remove the element to not have to worry about changing index ; also used to increment the normalized control coordinate by 1
            if len(critical_list)>=1:
                crit=critical_list[0] #
            else:
                crit=len(coord)-1 #on the last segment, the last point isn't in the critical point list
    if verbose==True:
        print("Original control coordinate:",coord)
        print("Normalized control coordinate:",normalized_coord)
    print("Normalization successfully completed.")
    return normalized_coord

"""
###########################################################################################################
FUNCTION: norm_diff
        Outputs the critical points of the original profile corresponding to the critical points of the difference profile.         
        In the case where the difference profile is non-monotonous, useful to change the critical points.
          
INPUT: coord1, coord2, coord_sampled, crit_list_diff2, crit_list1_0N, crit_list2_0N, verbose=False
    coord1               = list of the control coordinate of the first profile
    coord2               = list of the control coordinate of the second profile
    coord_sampled        = list of the control coordinate of the fitted profile
    crit_list_diff2      = list of the critical points of the difference of the energy profiles
    crit_list1_0N        = list of the critical points of the first energy profile (including first and last points)
    crit_list2_0N        = list of the critical points of the second energy profile (including first and last points)
    verbose              = False (default) : Boolean to ask the procedure to say everything (True) or be quiet (False)
    
OUTPUT: [crit_list_redone1,crit_list_redone2]
    crit_list_redone1 = list of the points of the first profile matching the critical points of the difference of the energy profiles
    crit_list_redone2 = list of the points of the first profile matching the critical points of the difference of the energy profiles

ERROR:
    "The critical point indices do not match the length of the control coordinate list": 
        the indices of the critical points are out of range of the control coordinate list
    "The critical point list length does not match the length of the control coordinate list": 
        there are more critical points than points in the control coordinate list

Added by Aël Cador for REGdiff (29-05-21)
###########################################################################################################
"""

def norm_diff(coord1, coord2, coord_sampled, crit_list_diff2, crit_list1_0N, crit_list2_0N, verbose=False):
    print("Starting the search for the points of the original profiles corresponding to the difference critical points...")
    
    crit_list_redone1=[]
    crit_list_redone2=[]
    for i in range(len(crit_list_diff2)):
        coord = coord_sampled[crit_list_diff2[i]] #where is the new critical point
        seg_number = int(coord) #since the nth segment (in 0-n-1 notation) is from 
        rem = coord - seg_number
        
        start1 = coord1[crit_list1_0N[seg_number]]
        end1 = coord1[crit_list1_0N[seg_number+1]]
        coord_corr1 = start1 + rem*(end1-start1)
        cred1 = np.array(coord1) - coord_corr1
        crit_corr1 = np.argmin(np.fabs(cred1))
        crit_list_redone1.append(crit_corr1)
        crit_list_redone1 = sorted(list(set(crit_list_redone1)))
        
        start2 = coord2[crit_list2_0N[seg_number]]
        end2 = coord2[crit_list2_0N[seg_number+1]]
        coord_corr2 = start2 + rem*(end2-start2)
        cred2 = np.array(coord2) - coord_corr2
        crit_corr2 = np.argmin(np.fabs(cred2))
        crit_list_redone2.append(crit_corr2)
        crit_list_redone2 = sorted(list(set(crit_list_redone2)))
    
    if verbose==True:
        print("Difference critical points:",crit_list_diff2)
        print("Corresponding points in profile 1:",crit_list_redone1)
        print("Corresponding points in profile 2:",crit_list_redone2)
    print("Critical point search successfully completed.")
    
    return([crit_list_redone1,crit_list_redone2])

"""
###########################################################################################################
FUNCTION: fit
          Fits a 4th-order polynomial to the term values with a resampling of the points 
          (the number of points equals the length of the longest segment), with corrections 
          of higher orders up to half the number of points.
          
INPUT: term, crit_list_0N, seg_long, coord, verbose=False
    term            = array of the term over the profile
    crit_list_0N    = list of the critical points of the energy profile preceded by 0 and followed by the total number of points - 1
    maxlong         = the length of the longest segment
    coord           = list of the control coordinate
    verbose         = boolean (True to ask for all details (orders, rmse, correction steps))
    
OUTPUT: [term_fit, xn_tot]
    term_fit    = the fitted polynomial values
    xn_tot      = the complete sampling abscissa

ERROR:
    Add an error check or at least a warning when the fit isn't good (instead of the fabulous error messages of numpy)
    (Desactivated at the moment)

Added by Aël Cador for REGdiff (03-03-21)
###########################################################################################################
"""

def fit(term, crit_list_0N, maxlong, coord, verbose=False):
    print("Starting fitting...")
    warnings.simplefilter('ignore') #just to ignore errors with "bad" fits
    
    order=4
    #print(type(term))
    if type(term) == np.ndarray:
        dim=term.ndim
    else:
        dim=0
    xn_tot=np.linspace(coord[crit_list_0N[0]],coord[crit_list_0N[1]],num=maxlong)
    if dim==1 and term.dtype != object:#and term.dtype!=np.float64 :
        term_fit_temp=np.zeros(((len(crit_list_0N)-1),maxlong))
        for i in range(len(crit_list_0N)-1):
            npoints=crit_list_0N[i+1]+1 - crit_list_0N[i] # number of points on the segment
            xn=np.linspace(coord[crit_list_0N[i]],coord[crit_list_0N[i+1]],num=maxlong) #points at which the polynomial will be evaluated
            rmse=0
            #print(len(coord[crit_list_0N[i]:crit_list_0N[i+1]+1]),len(term[crit_list_0N[i]:crit_list_0N[i+1]+1]))
            fit = poly.polyfit(coord[crit_list_0N[i]:crit_list_0N[i+1]+1],term[crit_list_0N[i]:crit_list_0N[i+1]+1],order) #Fit of the 4th order polynomial (change the order at the very end)
            if verbose==True:
                print("No correction",i,order)
            term_fit_tempk = poly.polyval(coord[crit_list_0N[i]:crit_list_0N[i+1]+1],fit)#[0:10]
            for j in range(npoints):
                rmse+=abs(term_fit_tempk[j]-term[crit_list_0N[i]+j])
            rmse=2625.5*rmse/npoints #RMSE in kJ/mol
            if verbose==True:
                print("Energy fit uncorrected",order,rmse)
            
            if rmse>0.01:
                mael=[]
                for order1 in range(1,npoints-3):
                    rmse=0
                    fit1 = poly.polyfit(coord[crit_list_0N[i]:crit_list_0N[i+1]+1],term[crit_list_0N[i]:crit_list_0N[i+1]+1],order1) #Fit of the 4th order polynomial (change the order at the very end)
                    if verbose==True:
                        print(i,order1)
                    term_fit_tempk = poly.polyval(coord[crit_list_0N[i]:crit_list_0N[i+1]+1],fit1)#[0:10]
                    for j in range(npoints):
                        rmse+=pow(term_fit_tempk[j]-term[crit_list_0N[i]+j],2)
                    rmse=2625.5*rmse/(2*npoints) #RMSE in kJ/mol
                    mael.append(rmse)
                    if verbose==True:
                        print(order1,rmse)
                order1=1+mael.index(min(mael))
                if verbose==True:
                    print("After correction",i,order1)
                    print("Energy fit corrected (1) - Chosen order:",order1," (up to order: ",npoints-1,") ",min(mael))
                """ #some corrections that caused problems and that I would like to rework someday.
                if min(mael)>0.01:
                    maell1,maelo1=[],[]
                    maell2,maelo2=[],[]
                    for rmid in range(2,npoints-2):
                        mid=crit_list_0N[i]+rmid
                        term1=term[crit_list_0N[i]:mid+1].copy()
                        term2=term[mid:crit_list_0N[i+1]+1].copy()
                        mael1,mael2=[],[]
                        for order1 in range(1,rmid):
                            mae1=0
                            fit1 = poly.polyfit(coord[crit_list_0N[i]:mid+1],term1,order1)
                            term_fit_tempk1 = poly.polyval(coord[crit_list_0N[i]:mid+1],fit1)
                            for j in range(rmid):
                                mae1+=abs(term_fit_tempk1[j]-term1[j])
                            mae1=2625.5*mae1/(rmid)
                            #print(order1,mae1)
                            mael1.append(mae1)
                        #print(order1,min(mael1))
                        maell1.append(min(mael1))
                        maelo1.append(1+mael1.index(min(mael1)))
                        for order2 in range(1,npoints-rmid):
                            mae2=0
                            fit2 = poly.polyfit(coord[mid:crit_list_0N[i+1]+1],term2,order2)
                            term_fit_tempk2 = poly.polyval(coord[mid:crit_list_0N[i+1]+1],fit2)
                            for j in range(crit_list_0N[i+1]+1-mid):
                                mae2+=abs(term_fit_tempk2[j]-term2[j])
                            mae2=2625.5*mae2/(npoints-rmid)
                            #print(order2,mae2)
                            mael2.append(mae2)
                        #print(order2,min(mael2))
                        maell2.append(min(mael2))
                        maelo2.append(1+mael2.index(min(mael2)))
                    mid=crit_list_0N[i]+2+maell1.index(min(maell1))
                    order1=maelo1[mid-2-crit_list_0N[i]]
                    order2=maelo2[mid-2-crit_list_0N[i]]
                    rmid=2+maell1.index(min(maell1))
                    print("Et donc (2)",rmid,mid,order1,min(maell1),order2,min(maell2))
                    term1=term[crit_list_0N[i]:mid+1].copy()
                    term2=term[mid:crit_list_0N[i+1]+1].copy()
                    fit1 = poly.polyfit(coord[crit_list_0N[i]:mid+1],term1,order1)
                    fit2 = poly.polyfit(coord[mid:crit_list_0N[i+1]+1],term2,order2)
                    xn1=np.linspace(coord[crit_list_0N[i]],coord[mid],num=int(maxlong*(rmid/npoints)))
                    xn2=np.linspace(coord[mid],coord[crit_list_0N[i+1]],num=maxlong-int(maxlong*(rmid/npoints)))
                    term_fit_temp1 = poly.polyval(xn1,fit1)
                    term_fit_temp2 = poly.polyval(xn2,fit2)
                    term_fit_temp[i] = np.concatenate((term_fit_temp1,term_fit_temp2))
                    """
                #else:
                fit2 = poly.polyfit(coord[crit_list_0N[i]:crit_list_0N[i+1]+1],term[crit_list_0N[i]:crit_list_0N[i+1]+1],order1) #Fit of the 4th order polynomial (change the order at the very end)
                term_fit_temp[i] = poly.polyval(xn,fit2)#[0:10]
            
            else:
                term_fit_temp[i] = poly.polyval(xn,fit)
            if i>0:
                xn_tot=np.concatenate((xn_tot,xn[1:]))
        if len(crit_list_0N)>2:
            for i in range(len(term_fit_temp)-1):
                term_fit_temp[i][-1] = 0.5*(term_fit_temp[i][-1]+term_fit_temp[i+1][0])
            term_fit=np.concatenate((term_fit_temp[0],term_fit_temp[1][1:]))
            for i in range(2,len(crit_list_0N)-1):
                term_fit=np.concatenate((term_fit,term_fit_temp[i][1:]))
    elif dim == 2 and term[0].ndim == 1:
        order=4
        term_fit_temp=np.zeros(((len(crit_list_0N)-1),len(term),maxlong)) #uses 5 times the number of points of the longest segment for all segment sampling
        for i in range(len(crit_list_0N)-1):
            npoints=crit_list_0N[i+1]+1 - crit_list_0N[i]
            xn=np.linspace(coord[crit_list_0N[i]],coord[crit_list_0N[i+1]],num=maxlong) #points at which the polynomial will be evaluated
            for j in range(len(term)):
                rmse=0
                fit = poly.polyfit(coord[crit_list_0N[i]:crit_list_0N[i+1]+1],term[j][crit_list_0N[i]:crit_list_0N[i+1]+1],order) #Fit of the 4th order polynomial (change the order at the very end)
                if verbose==True:
                    print("No correction",i,j,order)
                term_fit_tempk = poly.polyval(coord[crit_list_0N[i]:crit_list_0N[i+1]+1],fit)#[0:10]
                for k in range(npoints):
                    rmse+=pow(term_fit_tempk[k]-term[j][crit_list_0N[i]+k],2)
                rmse=2625.5*rmse/(npoints*abs(max(term[j][crit_list_0N[i]:crit_list_0N[i+1]+1])-min(term[j][crit_list_0N[i]:crit_list_0N[i+1]+1]))) #RMSE in kJ/mol
                #if i==0:
                    #print(i,j,rmse)
                    
                if rmse > 0.0001:
                    if verbose==True:
                        print(i,j,rmse,"Refitting...")
                    
                    mael=[]
                    if npoints<8:
                        maxo=5
                    else:
                        maxo=int(npoints/2)
                    for order1 in range(1,maxo):
                        rmse=0
                        fit1 = poly.polyfit(coord[crit_list_0N[i]:crit_list_0N[i+1]+1],term[j][crit_list_0N[i]:crit_list_0N[i+1]+1],order1) #Fit of the 4th order polynomial (change the order at the very end)
                        term_fit_templ = poly.polyval(coord[crit_list_0N[i]:crit_list_0N[i+1]+1],fit1)#[0:10]
                        for k in range(npoints):
                            rmse+=pow(term_fit_templ[k]-term[j][crit_list_0N[i]+k],2)
                        rmse=2625.5*rmse/npoints #MAE in kJ/mol
                        if verbose==True:
                            print(order1,rmse)
                        mael.append(rmse)
                    order1=1+mael.index(min(mael))
                    if verbose==True:
                        print("Chosen order:",order1," (up to order: ",maxo-1,") ",min(mael))
                    fit2 = poly.polyfit(coord[crit_list_0N[i]:crit_list_0N[i+1]+1],term[j][crit_list_0N[i]:crit_list_0N[i+1]+1],order1) #Fit of the 4th order polynomial (change the order at the very end)
                    if verbose==True:
                        print("With correction",i,j,order1)
                    term_fit_temp[i][j] = poly.polyval(xn,fit2)#[0:10]
                  
                else:
                    term_fit_temp[i][j] = poly.polyval(xn,fit)
            if i>0:
                xn_tot=np.concatenate((xn_tot,xn[1:]))
        if len(crit_list_0N)>2:               
            term_fit=[]
            for j in range(len(term_fit_temp[0])):
                for i in range(len(term_fit_temp)-1):
                    #print(i,j)
                    if i==0:
                        temp=(term_fit_temp[i][j].tolist()+term_fit_temp[i+1][j][1:].tolist())                            
                    else:
                        temp=(temp+(term_fit_temp[i+1][j][1:].tolist()))
                    #print(temp)
                term_fit.append(temp)
            term_fit=np.array(term_fit)
    else: #nbo
        order=4
        term_fit=[] #np.zeros(((len(crit_list_0N)-1),len(term),maxlong)) #uses 5 times the number of points of the longest segment for all segment sampling
        for i in range(len(term)):
            npoints = crit_list_0N[i+1]+1 - crit_list_0N[i]
            term_fit_temp = []
            for j in range(len(term[i])): #len(crit_list_0N)-1):
                xn = np.linspace(coord[crit_list_0N[i]],coord[crit_list_0N[i+1]],num=maxlong) #points at which the polynomial will be evaluated
                #for j in range(len(term[h][i])):
                rmse=0
                #print(len(coord[crit_list_0N[i]:crit_list_0N[i+1]+1]),len(term[i][j]))
                fit = poly.polyfit(coord[crit_list_0N[i]:crit_list_0N[i+1]+1],term[i][j],order) #[j][crit_list_0N[i]:crit_list_0N[i+1]+1],order) #Fit of the 4th order polynomial (change the order at the very end)
                if verbose==True:
                    print("No correction",i,j,order)
                term_fit_tempk = poly.polyval(coord[crit_list_0N[i]:crit_list_0N[i+1]+1],fit)#[0:10]
                for k in range(npoints):
                    rmse+=pow(term_fit_tempk[k]-term[i][j][k],2)
                rmse=2625.5*rmse/(npoints*abs(max(term[i][j])-min(term[i][j]))) #RMSE in kJ/mol
                #if i==0:
                    #print(i,j,rmse)
                    
                if rmse > 0.0001:
                    if verbose==True:
                        print(i,j,rmse,"Refitting")
                    
                    mael=[]
                    if npoints<8:
                        maxo=5
                    else:
                        maxo=int(npoints/2)
                    for order1 in range(1,maxo):
                        rmse=0
                        fit1 = poly.polyfit(coord[crit_list_0N[i]:crit_list_0N[i+1]+1],term[i][j],order1) #Fit of the 4th order polynomial (change the order at the very end)
                        term_fit_templ = poly.polyval(coord[crit_list_0N[i]:crit_list_0N[i+1]+1],fit1)#[0:10]
                        for k in range(npoints):
                            rmse+=pow(term_fit_templ[k]-term[i][j][k],2)
                        rmse=2625.5*rmse/npoints #MAE in kJ/mol
                        if verbose==True:
                            print(order1,rmse)
                        mael.append(rmse)
                    order1=1+mael.index(min(mael))
                    if verbose==True:
                        print("Chosen order:",order1," (up to order: ",maxo-1,") ",min(mael))
                    fit2 = poly.polyfit(coord[crit_list_0N[i]:crit_list_0N[i+1]+1],term[i][j],order1) #Fit of the 4th order polynomial (change the order at the very end)
                    term_fit_temp.append(poly.polyval(xn,fit2))
                  
                else:
                    term_fit_temp.append(poly.polyval(xn,fit))
            term_fit_temp = np.array(term_fit_temp)
            term_fit.append(term_fit_temp)
        term_fit = np.array(term_fit)
    
    print("Fitting successfully completed.")
    return [term_fit, xn_tot]

"""
###########################################################################################################
FUNCTION: split_segm
          Takes the A array and divide it into N arrays according with the number of 
          critical points. Each array corresponds to a segment of the REG analysis.
          
INPUT: A, critical_point_list, verbose=False
    A                       = Array containing the abcissa values of a function 
    critical_points_list    = list of the index of the critical points (zero indexed)
    verbose                 = False (default) : Boolean to ask the procedure to say everything (True) or be quiet (False)
 
OUTPUT: segm
    segm = Array of arrays containing the values of A for each segment    

ERROR:
    "Invalid critical point index" : index of critical point out of range in the function array.
###########################################################################################################
"""

def split_segm(A, critical_point_list, verbose=False):
    if verbose==True:
        print("Splitting segment...")
    
    critical_point_list = list(set(critical_point_list)) ## Remove duplicates and sort in crescent order
    #ERRORS:
    if critical_point_list[-1] >= len(A): ## Checks if index of critical point is in range of function
        raise ValueError("Invalid critical point index") 

    #INTERNAL VARIABLES:
    segm = [] ## array of arrays
    start = 0
    end = len(A)
    
    #DIVIDE FUNCTION INTO SEGMENTS:
    while start < sorted(critical_point_list)[-1]:
        for i in sorted(critical_point_list):
            segm.append([A[j] for j in range(start, i+1)])
            start = i
    segm.append([A[j] for j in range(start, end)])
    
    if verbose==True:
        print("Original segment:",A)
        print("Split segments:",segm)
        print("Segment successfully split.")
    return segm

"""
###########################################################################################################
FUNCTION: sign_segm
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

Added by Aël Cador (09-04-21)
###########################################################################################################
"""

def sign_segm(term,crit_list, verbose=False):
    print("Starting sign check...")
    
    if crit_list[-1] >= len(term[0]): ## Checks if index of critical point is in range of function
        raise ValueError("Invalid critical point index")
    split_term=[]
    for A in term:
        split_term.append(split_segm(A, crit_list))
    #split_intra1=reg.split_segm(iqa_intra1, crit_list1)
    sign_seg_term=[]
    for j in range(len(split_term[0])):
        sign_seg_term_j=[]
        for i in range(len(split_term)):
            sign_seg_term_i=[]
            for k in range(len(split_term[i][j])):
                if split_term[i][j][k]>0:
                    sign_seg_term_i.append(1)
                else:
                    sign_seg_term_i.append(-1)
            if sum(sign_seg_term_i) == len(split_term[i][j]):
                sign_seg_term_j.append(1)
            elif sum(sign_seg_term_i) == -len(split_term[i][j]):
                sign_seg_term_j.append(-1)
            else:
                sign_seg_term_j.append(0)
        sign_seg_term.append(sign_seg_term_j)
    
    if verbose==True:
        print(sign_seg_term)
    print("Sign check completed.")
    return sign_seg_term
    
"""
###########################################################################################################
FUNCTION: reg
          Perform the REG analysis over all contributions inside "terms" array
          
INPUT:  wfn_energy, control_coord, terms, critical, np, inflex, critical_index, mode, verbose=False
    wfn_energy      = Array containing the energy values for each point (PES - Energy)
    control_coord   = Array containing the control coodinate values for each point (PES - Coordinate)
    terms           = Array of arrays, each one corresponding to one IQA contribution. 
    critical        = True (default):       Search for critical points
                    = False:                user must provide the critical point index array
    np              = 5 (default):          minimum amount of points between two critical points (used if critical == True)
    inflex      = False (default):      Do not search for inflexion points (second derivative = 0) (used if critical == True)
                    = True:                 Search for inflexion points (second derivative = 0 (used if critical == True)
    critical_index  = [] (default):         list of critical points provided by the user. (used if critical == False)
    mode            = None:                 for regular linear regression
                    = "norm" (default):     use normalized values of A and B
                    = "std":                use standardized values for A and B
    verbose         = False (default):      standard output (quiet)
                    = True:                 output everything
            
OUTPUT: [REG_values/pearson_values][segment]
    REG_values = array of REG coefficients for each contribution. Each array correspond to a segment
    pearson_values = array of Pearson coefficients for each contribution. Each array correspond to a segment

ERROR:
     "PES abcissa and ordenate must have same number of points" : wfn_energy and control_coord have differnt sizes
     "Contributions arrays must have same size": arrays inside terms have diferent sizes
###########################################################################################################
"""
    
def reg(wfnx_energy, control_coord, terms, critical=True, np=5, inflex=True, critical_index=[], mode="norm", verbose=False):
    print("STARTING REG ANALYSIS...")
    
    #ERRORS:
    if len(wfnx_energy) != len(control_coord): # Checks if abcissa and ordenate of PES have same size
        raise ValueError("PES abcissa and ordenate must have same number of points")  
    for contribution in terms:
        if len(wfnx_energy) != len(contribution):  # Checks if all contributions arrays have same size
            raise ValueError("Contributions arrays must have same size")
    
    #INTERNAL VARIABLES:      
    REG_values = []
    pearson_values = []
    split_terms = []  #split_terms[contribution][segment]
    split_wfnx = []
    temp_REG = []
    temp_pearson = []

    #FIND CRITICAL POINTS:
    if critical == True: # Search for critical points automatically
        critical_points_list = find_critical(wfnx_energy, control_coord, min_points=np, use_inflex=inflex, verbose=verbose)
    
    else:
        critical_points_list = critical_index #User provides the critical points index

    #DIVIDE WFN ENERGY INTO SEGMENTS:
    split_wfnx = split_segm(wfnx_energy, critical_points_list, verbose=verbose)

    #DIVIDE TERMS INTO SEGMENTS:
    if verbose==True:
        print("First term:",terms[0])
    print("Starting segment splitting...")
    for term in terms:
        split_terms.append(split_segm(term, critical_points_list, verbose=verbose))
    print("Segment splitting successfully completed.")
    if verbose==True:
        print("First split term:",split_terms[0])  
    #LINEAR REGRESSION FOR CONTRIBUTION INSIDE EACH SEGMENT:
    for i in range(len(critical_points_list)+1):
        for split_term in split_terms:
            temp_REG.append(regression(split_wfnx[i], split_term[i], mode=mode)[0])
            temp_pearson.append(regression(split_wfnx[i], split_term[i], mode=mode)[2])
        REG_values.append(temp_REG)
        pearson_values.append(temp_pearson)
        #print(temp_REG,temp_pearson)
        temp_REG = []
        temp_pearson = []
    
    if verbose==True:
        print("REG values of first term:",REG_values[0])
        print("R² values of first term:",pearson_values[0])
    print("REG ANALYSIS COMPLETED.")
    return REG_values, pearson_values

"""
###########################################################################################################
FUNCTION: regnbo
          Perform the REG-NBO analysis over all contributions inside "terms" array
          
INPUT:  wfn_energy, control_coord, terms, critical, np, inflex, critical_index, mode, verbose=False
    wfn_energy      = Array containing the energy values for each point (PES - Energy)
    control_coord   = Array containing the control coodinate values for each point (PES - Coordinate)
    terms           = Array of arrays, each one corresponding to one IQA contribution. 
    critical        = True (default):       Search for critical points
                    = False:                user must provide the critical point index array
    np              = 5 (default):          minimum amount of points between two critical points (used if critical == True)
    inflex      = False (default):      Do not search for inflexion points (second derivative = 0) (used if critical == True)
                    = True:                 Search for inflexion points (second derivative = 0 (used if critical == True)
    critical_index  = [] (default):         list of critical points provided by the user. (used if critical == False)
    mode            = None:                 for regular linear regression
                    = "norm" (default):     use normalized values of A and B
                    = "std":                use standardized values for A and B
    verbose         = False (default):      standard output (quiet)
                    = True:                 output everything
           
OUTPUT: [REG_values/pearson_values][segment]
    REG_values = array of REG coefficients for each contribution. Each array correspond to a segment
    pearson_values = array of Pearson coefficients for each contribution. Each array correspond to a segment

ERROR:
     "PES abcissa and ordenate must have same number of points" : wfn_energy and control_coord have differnt sizes
     "Contributions arrays must have same size": arrays inside terms have diferent sizes

Added by Aël Cador for REGNBO (28-04-21)
###########################################################################################################
"""
    
def regnbo(wfnx_energy, control_coord, terms, critical=True, np=5, inflex=True, critical_index=[], mode="norm", verbose=False):
    print("STARTING REG-NBO ANALYSIS...")
    
    #ERRORS:
    if len(wfnx_energy) != len(control_coord): # Checks if abcissa and ordenate of PES have same size
        raise ValueError("PES abcissa and ordenate must have same number of points")  
    #for contribution in terms:
    
    #INTERNAL VARIABLES:      
    REG_values = []
    pearson_values = []
    split_wfnx = []
    temp_REG = []
    temp_pearson = []

    #FIND CRITICAL POINTS:
    if critical == True: # Search for critical points automatically
        critical_points_list = find_critical(wfnx_energy, control_coord, min_points=np, use_inflex=inflex, verbose=False)    
    else:
        critical_points_list = critical_index #User provides the critical points index
    """
    term_len=0
    for i in range(len(terms)):
        print(i)
        if terms[i].size != 0:
            if terms[i][0].size == 0:
                print("Empty element at index ",i,0)
            else:
                term_len += len(terms[i][0])
        else:
            print("Empty element at index ",i)
        print(len(wfn_energy),len(critical_index),term_len)
    if len(wfn_energy)+len(critical_index) != term_len:  # Checks if all contributions arrays have same size ; +len(critical_index) because the 
        raise ValueError("Contributions arrays must have same size")
    """
    #DIVIDE WFN ENERGY INTO SEGMENTS:
    split_wfnx = split_segm(wfnx_energy, critical_points_list, verbose=False)
  
    #LINEAR REGRESSION FOR CONTRIBUTION INSIDE EACH SEGMENT:
    for i in range(len(critical_points_list)+1):
        for j in range(len(terms[i])):
            #print(i,j)
            temp_REG.append(regression(split_wfnx[i], terms[i][j], mode=mode)[0])
            temp_pearson.append(regression(split_wfnx[i], terms[i][j], mode=mode)[2])
        REG_values.append(temp_REG)
        #print(temp_REG,temp_pearson)
        pearson_values.append(temp_pearson)
        temp_REG = []
        temp_pearson = []
    
    if verbose==True:
        print("REG values of first term:",REG_values[0])
        print("R² values of first term:",pearson_values[0])
    print("REG ANALYSIS COMPLETED.")
    return REG_values, pearson_values

"""
###########################################################################################################
FUNCTION: reg_zeroes
          removes terms with null REG coefficients (constant terms)

INPUT: REG_table, term_list, verbose=False
    REG_table       = table of both REG coefficients and R², resulting from reg.reg or reg.regnbo
    term_list       = list containing all terms which have to be analyzed (e.g.:
                           - Intra OR Inter OR Disp terms and header for REG-IQA, 
                           - E(2), E(2) header, dE, F OR E(4), E(4) header, and S for REG-NBO). 
                    Probably easier to read with header first.
    verbose         = False (default):      standard output (quiet)
                    = True:                 output everything       

OUTPUT: [error, RMSE]
    error   = list of wfn_energies - IQA_energies
    RMSE    = Root mean square error
    
ERROR:
    "Arrays must have the same size" : wfn_energies and IQA_energies have different sizes.
    "Segments must have the same size"
Added by Aël Cador for REGNBO (21-09-21)
###########################################################################################################
"""

def reg_zeroes(REG_table, term_list, verbose=False):
    print("Starting searching for zero-REG terms...")
    indice_removal=[]
    for i in range(len(REG_table[0])):
        for j in range(len(REG_table[0][i])):
            if REG_table[0][i][j] == 0:
                if verbose==True:
                    print("Zero-REG term at indices",i,j,[term[i][j] for term in term_list])
                    print("REG coefficient:",REG_table[0][i][j],"R²:",REG_table[1][i][j])
                indice_removal.append([i,j])
    indice_removal.reverse()
    for i in range(len(indice_removal)):
        REG_table[0][indice_removal[i][0]].pop(indice_removal[i][1])   #REG coefficients
        REG_table[1][indice_removal[i][0]].pop(indice_removal[i][1])   #R²
        for term in term_list:
            term[indice_removal[i][0]].pop(indice_removal[i][1])
    
    if verbose==True:
        print("Removed",len(indice_removal),"terms.")
        if len(indice_removal)>0:
            print("Indices removed:",indice_removal.reverse())     
    print("Search for zero-REG terms completed.")
    
    return REG_table, term_list

"""
###########################################################################################################
FUNCTION: integration_error
          calculates integration error for each PES point

INPUT: wfn_energies, IQA_energies
    wfn_energies = list of total energies from wfn files
    IQA_energies = list of total energy obtained from IQA terms       

OUTPUT: [error, RMSE]
    error   = list of wfn_energies - IQA_energies
    RMSE    = Root mean square error
    
ERROR:
    "Arrays must have the same size" : wfn_energies and IQA_energies have different sizes.
    "Segments must have the same size"
###########################################################################################################
""" 

def integration_error(wfn_energies, IQA_energies):
    #ERRORS:  
    if len(wfn_energies) != len(IQA_energies): ## Checks if X and Y have the same size
        raise ValueError("Arrays must have the same size")
    if wfn_energies.ndim != IQA_energies.ndim: ## Checks if X and Y have the same size
        raise ValueError("Segments must have the same size")
    if wfn_energies.ndim == 2 and len(wfn_energies[0]) != len(IQA_energies[0]): ## Checks if X and Y have the same size
        raise ValueError("Segments must have the same size")
   
    if wfn_energies.ndim == 1 and IQA_energies.ndim == 1:
        #CREATE LIST FOR EACH POINT ERROR    
        error = [wfn_energies[i] - IQA_energies[i] for i in range(len(wfn_energies))] # error values list
        #CALCULATE RMSE
        temp = [a**2 for a in error]
        RMSE = 2625.5*(sum(temp)/len(temp))**0.5
    if wfn_energies.ndim == 2 and IQA_energies.ndim == 2:
        error=[]
        temp=[]
        for i in range(len(wfn_energies)):
            error1 = [wfn_energies[i][j] - IQA_energies[i][j] for j in range(len(wfn_energies[i]))]
            error.append(error1)
            temp1 = [a**2 for a in error1]
            temp.append(sum(temp1))
        RMSE = 2625.5*(sum(temp)/(len(wfn_energies)*len(wfn_energies[0])))**0.5
        
    return error, RMSE
    

"""
###########################################################################################################
FUNCTION: group_intra
          group IQA intra-atomic terms into the user-defined groups
          
INPUT: intra_energies, intra_energies_header, groups, verbose=False
    intra_energies          = list of IQA intra-atomic terms
    intra_energies_header   = header of IQA intra atomic terms.
    groups                  = list of list of atomic labels, each internal list corresponds to a different group
    verbose                 = False (default) : Boolean to ask the procedure to say everything (True) or be quiet (False)
    
OUTPUT: [new_terms, new_header]
    new_terms   = list of the grouped IQA terms values
    new_header  = list of the grouped terms header
    
ERROR:
    
###########################################################################################################
"""

def group_intra(intra_energies, intra_energies_header, groups, verbose=False):
    print("Grouping Intra terms...")
    if verbose == True:
        print("Intra terms table length:",len(intra_energies),len(intra_energies[0]))
    
    #INTERNAL VARIABLES
    new_terms = []
    new_header = []
    temp_header = []
    
    for i in range(len(intra_energies_header)): # copy IQA energy values list
        temp_header.append(intra_energies_header[i])
    
    #CREATE GROUPS LABELS
    for group in groups:
        group_label = 'G' + str(groups.index(group)+1)
        for atom in group:
            size = -1*len(atom)
            for i in range(len(temp_header)):
                if atom in temp_header[i][size:]:
                    temp_header[i] = temp_header[i].replace(atom, group_label)
    
    #SUM TERMS WITH SAME GROUP LABEL        
    for term in temp_header:
        temp = 0
        if term not in new_header:
            new_header.append(term)
            for i in range(len(temp_header)):
                if term == temp_header[i]:
                    temp = temp + intra_energies[i]
            new_terms.append(temp)
    
    if verbose == True:
        print("Intra terms table length:",len(new_terms),len(new_terms[0]))
        
    print("Intra terms successfully grouped.")                
    return new_terms, new_header

"""
###########################################################################################################
FUNCTION: group_inter
          group IQA interatomic terms into the user-defined groups
          
INPUT: inter_energies, inter_energies_header, groups, verbose=False
    inter_energies          = list of IQA interatomic terms
    inter_energies_header   = header of IQA interatomic terms.
    groups                  = list of list of atomic labels, each internal list corresponds to a different group
    verbose                 = False (default) : Boolean to ask the procedure to say everything (True) or be quiet (False)
    
OUTPUT: [new_terms, new_header]
    new_terms  = list of the grouped IQA terms values
    new_header = list of the grouped terms header
    
ERROR:

AC modification: modified the criteria to search for redundant term headers since the term header format was
                changed from VC_IQA(A,B)-a1_a2 to VC-a1_a2, so it only checks if the term starts with VC or VX
                and a given pair of atoms. The old version is simply commented.
###########################################################################################################
"""

def group_inter(inter_energies, inter_energies_header, groups, verbose=False):
    print("Grouping Inter terms...")
    
    #INTERNAL VARIABLES
    temp_terms = []
    new_header = []
    temp_header = []
    temp_header2 = []
    new_terms = []
    redundancies = []

    for i in range(len(inter_energies_header)): # copy IQA energy values list
        temp_header.append(inter_energies_header[i])
  
    #CREATE GROUPS LABELS
    for group in groups:
        group_label = 'G' + str(groups.index(group)+1) 
        for atom in group:
            atom_label_1 = "-" + atom + "_"
            size = -1*len(atom)

            for i in range(len(temp_header)):
                if atom_label_1 in temp_header[i]:
                    temp_header[i] = temp_header[i].replace(atom_label_1, '-' + group_label + '_')
                    
            for i in range(len(temp_header)):
                if atom in temp_header[i][size:]:
                    temp_header[i] = temp_header[i].replace(atom, group_label)
    
    if verbose == True:
        print("temp_header",temp_header)
    
    #SUM IQA CONTRIBUTIONS OF SAME LABEL
    for term in temp_header:
        temp =0
        if term not in temp_header2: #"-"+term+"_"
            temp_header2.append(term)
            for i in range(len(temp_header)):
                if term == temp_header[i]:
                    temp = temp + inter_energies[i]
            temp_terms.append(temp)
    if verbose == True:
        print("temp_header2",temp_header2)
        
    # ORGANIZE AND REMOVE REDUNDANCIES            
    for i in range(len(temp_header2)):
        temp = 0
        if temp_header2[i] not in new_header:
            if temp_header2[i] not in redundancies:
                new_header.append(temp_header2[i])
                if verbose == True:
                    print("new_header:",new_header)
                temp = temp + temp_terms[i]
                group_label_1 = temp_header2[i][-14:][temp_header2[i][-14:].find('-')+1:temp_header2[i][-14:].find('_')]
                group_label_2 = temp_header2[i][-14:][temp_header2[i][-14:].find('_')+1:]
                if verbose == True:
                    print("Group labels:",group_label_1,group_label_2)
                for j in range(len(temp_header2)):
                    #FOR VC_IQA(A,B) STYLE: if temp_header2[i][:-1*temp_header2[i].find('(')] == temp_header2[j][:-1*temp_header2[j].find('(')]:
                    if temp_header2[i][:2] == temp_header2[j][:2]:
                        if temp_header2[j] not in new_header:
                            if temp_header2[j] not in redundancies:
                                if group_label_1 == temp_header2[j] and group_label_2 == temp_header2[j] and group_label_2 != group_label_1:
                                    temp = temp + temp_terms[j]
                                    redundancies.append(temp_header2[j])
                                    if verbose == True:
                                        print("redundancies",redundancies)
                                        print("temp_header2[j]",temp_header2[j])
                new_terms.append(temp)
    
    if verbose == True:
        print("Inter terms table length:",len(new_terms),len(new_terms[0]))
    print("Inter terms successfully grouped.")
    
    return  new_terms, new_header

"""
###########################################################################################################
FUNCTION: group_vcvx
          group Vc and Vx IQA interatomic terms
          
INPUT: iqa_inter_g, iqa_inter_header_g, verbose=False
    iqa_inter_g         = list of IQA interatomic terms
    iqa_inter_header_g  = header of IQA interatomic terms.
    verbose             = False (default) : Boolean to ask the procedure to say everything (True) or be quiet (False)
    
OUTPUT: [iqa_inter_g2,iqa_inter_header_g2]
    iqa_inter_g2        = list of the grouped IQA terms values
    iqa_inter_header_g2 = list of the grouped terms header
    
ERROR:

Added by Aël Cador (12-04-21)    
###########################################################################################################
"""

def group_vcvx(iqa_inter_g, iqa_inter_header_g, verbose=False):
    print("Grouping Vc and Vx terms...")
    
    len_inter = int(len(iqa_inter_g)/2)
    iqa_inter_g2 = []
    for i in range(len_inter):
        iqa_inter_g2.append(iqa_inter_g[i]+iqa_inter_g[i+len_inter])
    iqa_inter_header_gc = iqa_inter_header_g[:len_inter].copy()
    iqa_inter_g2 = np.array(iqa_inter_g2)
    iqa_inter_header_g2 = []
    for i in range(len(iqa_inter_header_gc)):
        iqa_inter_header_g2.append(iqa_inter_header_gc[i].replace("Vcl","Inter"))
    
    if verbose==True:
        print("Inter term table first element:",iqa_inter_g2[0])
        print("Inter header table first element:",iqa_inter_header_g2[0])
    print("Vc and Vx terms successfully grouped.")
    
    return [iqa_inter_g2, iqa_inter_header_g2]

"""
###########################################################################################################
(Created a procedure to do the .csv and .xlsx files quickly instead of having it done repeatedly in the main script)
FUNCTION: csv
          generates .csv and .xlsx files containing the results of the REG analysis
          
INPUT: reg_inter, iqa_inter_header, reg_intra, iqa_intra_header, reg_disp, iqa_disp_header, SYS, path, g, sign=None, verbose=False
    reg_inter           = array of the interatomic terms with the REG and R² coefficients
    iqa_inter_header    = header (term names) of the interatomic terms
    reg_intra           = array of the intra-atomic terms with the REG and R² coefficients
    iqa_intra_header    = header (term names) of the intra-atomic terms
    reg_disp            = array of the dispersion terms with the REG and R² coefficients
    iqa_disp_header     = header (term names) of the dispersion terms
    SYS                 = name of the system considered
    path                = path to which the .csv file 
    g                   = bool:
                            True = grouped values, 
                            False = ungrouped values
    sign                = array of the sign of terms (optional)
    verbose             = False (default) : Boolean to ask the procedure to say everything (True) or be quiet (False)

OUTPUT: df_final, dataframe_list
    df_final = array of the X most important terms (in terms of REG, both positive and negative) of all segments
    dataframe_list = array of all terms on all segments ordered by REG coefficient

ERROR:
    (just think of creating the folder indicated by "path" before running the script, otherwise you get an error)

Added by Aël Cador (05-03-21)
###########################################################################################################
"""

def csv(reg_inter, iqa_inter_header, reg_intra, iqa_intra_header, reg_disp, iqa_disp_header, SYS, path, g, sign=None, verbose=False):
    print("Starting generating .csv result files...")
    if sign != None:
        sign_intra=sign[0]
        sign_inter=sign[1]
        sign_disp=sign[2]    
    dataframe_list=[]
    df_final = pd.DataFrame()
    
    for i in range(len(reg_inter[0])):
        if verbose==True:
            print("Generating Inter .csv file for segment",i)
        if sign != None:
            temp = [reg_inter[0][i], reg_inter[1][i], sign_inter[i]]                 #inter
            df_inter = pd.DataFrame(temp).transpose()                                #inter
            df_inter.columns = ["REG", "R2","SIGN"]                                  #inter
            #print(df_inter)
        else:                                     
            temp = [reg_inter[0][i], reg_inter[1][i]]                                #inter
            df_inter = pd.DataFrame(temp).transpose()                                #inter
            df_inter.columns = ["REG", "R2"]                                         #inter NB: AC change: R to R2. This will actually be R until squared
        df_inter.index = iqa_inter_header                                            #inter
        df_inter.to_csv(path +'/'+ SYS + "_Inter_seg_" +str(i+1) +".csv", sep=',')   #inter 
        df_inter = df_inter.rename_axis('TERM').reset_index()                        #inter
        
        if verbose==True:
            print("Generating Intra .csv file for segment",i)
        if sign != None:
            temp = [reg_intra[0][i], reg_intra[1][i], sign_intra[i]]                 #intra
            df_intra = pd.DataFrame(temp).transpose()                                #intra
            df_intra.columns = ["REG", "R2","SIGN"]                                  #intra
            #print(df_intra)
        else:
            temp = [reg_intra[0][i], reg_intra[1][i]]                                #intra
            df_intra = pd.DataFrame(temp).transpose()                                #intra
            df_intra.columns = ["REG", "R2"]                                         #intra NB: AC change: R to R2. This will actually be R until squared
        df_intra.index = iqa_intra_header                                            #intra
        df_intra.to_csv(path +'/'+ SYS + "_Intra_seg_" +str(i+1) + ".csv", sep=',')  #intra
        df_intra = df_intra.rename_axis('TERM').reset_index()                        #intra
        
        if verbose==True:
            print("Generating dispersion .csv file for segment",i)
        if sign != None:
            temp = [reg_disp[0][i], reg_disp[1][i], sign_disp[i]]                    #disp
            df_disp = pd.DataFrame(temp).transpose()                                 #disp
            df_disp.columns = ["REG", "R2","SIGN"]                                   #disp
            #print(df_disp)
        else:
            temp = [reg_disp[0][i], reg_disp[1][i]]                                  #disp
            df_disp = pd.DataFrame(temp).transpose()                                 #disp
            df_disp.columns = ["REG", "R2"]                                          #disp NB: AC change: R to R2. This will actually be R until squared
        df_disp.index = iqa_disp_header                                              #disp
        df_disp.to_csv(path +'/'+ SYS + "_Disp_seg_" +str(i+1) +".csv", sep=',')     #disp 
        df_disp = df_disp.rename_axis('TERM').reset_index()
        df_disp.dropna(axis=0, how='any', subset=None, inplace=True) #get rid of "NaN" terms which have a null standard deviation 
        
        df_temp = pd.concat([df_intra, df_inter, df_disp]).sort_values('REG')

        df_temp['R2'] = df_temp['R2']**2 #AC change: R to R2.
        dataframe_list.append(df_temp)
        
        df_temp2=df_temp.copy()
        #print('df_temp2.shape',df_temp2.shape)
        for i in range(df_temp2.shape[0]): #AC change: R² filter (>0.8)
            if df_temp2.iat[i,1]<0.8: #second column is R²
                df_temp2.iat[i,1]=np.nan
        df_temp2=df_temp2.dropna(axis='rows',how='any')
        df_temp3=pd.concat([df_temp2[-5:], df_temp2[:5]], axis=0).sort_values('REG', ascending=False).reset_index()
        #print(df_temp3)        
        df_final=pd.concat([df_final, df_temp3], axis=1)
    print(".csv result files successfully generated.")
    
    print("Starting generating .xlsx result files...")
    if g==False:   
         df_final.to_excel(path + '/' + SYS + '_IQA_segment_analysis.xlsx')
         with pd.ExcelWriter(path + '/' + SYS + '_complete_IQA_segment_analysis.xlsx') as writer: #writer object needed to have the lists written in different sheets
             for i in range(len(dataframe_list)):
                 dli=pd.DataFrame(dataframe_list[i]) #also need the data to be a dataframe to enable writing to XLSX file
                 dli.to_excel(writer, sheet_name='Segment'+str(i+1))
    else:
         df_final.to_excel(path + '/' + SYS + '_IQA_segment_analysis_g.xlsx')
         with pd.ExcelWriter(path + '/' + SYS + '_complete_IQA_segment_analysis_g.xlsx') as writer: #writer object needed to have the lists written in different sheets
             for i in range(len(dataframe_list)):
                 dli=pd.DataFrame(dataframe_list[i]) #also need the data to be a dataframe to enable writing to XLSX file
                 dli.to_excel(writer, sheet_name='Segment'+str(i+1))
    print(".xlsx result files successfully generated.")
    
    return dataframe_list

"""
###########################################################################################################
(Created a procedure to do the .csv and .xlsx files quickly instead of having it done repeatedly in the main script)
FUNCTION: csvnbo
          generates .csv and .xlsx files containing the results of the REG-NBO3 analysis (Gaussian default)

INPUT: reg_inter, iqa_inter_header, SYS, path, g, sign=None, verbose=False
    reg_inter           = array of the NBO interactions with the REG and R² coefficients
    iqa_inter_header    = header (term names) of the NBO interactions
    SYS                 = name of the system considered
    path                = path to which the .csv file 
    g                   = bool:
                            True = grouped segments, 
                            False = ungrouped segments
    sign                = array of the sign of terms (optional)
    verbose             = False (default) : Boolean to ask the procedure to say everything (True) or be quiet (False)

OUTPUT: df_final, dataframe_list
    df_final = array of the X most important terms (in terms of REG, both positive and negative) of all segments
    dataframe_list = array of all terms on all segments ordered by REG coefficient

ERROR:
    (just think of creating the folder indicated by "path" before running the script, otherwise you get an error)

Added by Aël Cador for REGNBO (28-04-21)
###########################################################################################################
"""

def csvnbo(reg_inter, iqa_inter_header, SYS, path, g, sign=None, verbose=False):
    print("Starting generating .csv result files...")
    
    if sign != None:
        sign_inter=sign[0]    
    dataframe_list=[]
    df_final = pd.DataFrame()
    for i in range(len(reg_inter[0])):
        if sign != None:
            temp = [reg_inter[0][i], reg_inter[1][i], sign_inter[i]]                 #inter
            df_inter = pd.DataFrame(temp).transpose()                                #inter
            df_inter.columns = ["REG", "R2","SIGN"]                                  #inter
            if verbose==True:
                print(df_inter)
        else:                                     
            temp = [reg_inter[0][i], reg_inter[1][i]]                                #inter
            if verbose==True:
                print("I",i,len(reg_inter),len(reg_inter[0]),len(reg_inter[0][i]),len(iqa_inter_header),len(iqa_inter_header[i]),
                  len(iqa_inter_header[i]))#[i]))
                print("TEMP",temp,reg_inter[0][i], reg_inter[1][i])
            df_inter = pd.DataFrame(temp).transpose()                                #inter
            df_inter.columns = ["REG", "R2"]                                         #inter NB: AC change: R to R2. This will actually be R until squared
        df_inter.index = iqa_inter_header[i]#[i]#[j] #0/i before                          #inter
        df_inter.to_csv(path +'/'+ SYS + "_NBO_seg_" +str(i+1) +".csv", sep=',')     #inter 
        df_inter = df_inter.rename_axis('TERM').reset_index()                        #inter
        
        df_temp = pd.concat([df_inter]).sort_values('REG')

        df_temp['R2'] = df_temp['R2']**2 #AC change: R to R2.
        dataframe_list.append(df_temp)
        if verbose==True:
            print("DATA")
            print(dataframe_list)
        df_temp2=df_temp.copy()
        if verbose==True:
            print('df_temp2.shape',df_temp2.shape)
        for i in range(df_temp2.shape[0]): #AC change: R² filter (>0.8)
            if df_temp2.iat[i,1]<0.8: #second column is R²
                df_temp2.iat[i,1]=np.nan
        df_temp2=df_temp2.dropna(axis='rows',how='any')
        df_temp3=pd.concat([df_temp2[-5:], df_temp2[:5]], axis=0).sort_values('REG', ascending=False).reset_index()
        if verbose==True:
            print("df_temp3",df_temp3)        
        df_final=pd.concat([df_final, df_temp3], axis=1)
        
    print(".csv result files successfully generated.")
    print("Starting generating .xlsx result files...")
    
    if g==False:   
         df_final.to_excel(path + '/' + SYS + '_NBO_segment_analysis.xlsx')
         with pd.ExcelWriter(path + '/' + SYS + '_complete_NBO_segment_analysis.xlsx') as writer: #writer object needed to have the lists written in different sheets
             for i in range(len(dataframe_list)):
                 dli=pd.DataFrame(dataframe_list[i]) #also need the data to be a dataframe to enable writing to XLSX file
                 dli.to_excel(writer, sheet_name='Segment'+str(i+1))
    else:
         df_final.to_excel(path + '/' + SYS + '_NBO_segment_analysis_3seg.xlsx')
         with pd.ExcelWriter(path + '/' + SYS + '_complete_NBO_segment_analysis_3seg.xlsx') as writer: #writer object needed to have the lists written in different sheets
             for i in range(len(dataframe_list)):
                 dli=pd.DataFrame(dataframe_list[i]) #also need the data to be a dataframe to enable writing to XLSX file
                 dli.to_excel(writer, sheet_name='Segment'+str(i+1))
    print(".xlsx result files successfully generated.")
    return dataframe_list

"""
###########################################################################################################
(Created a procedure to do the .csv and .xlsx files quickly instead of having it done repeatedly in the main script)
FUNCTION: csvnbo7
          generates .csv and .xlsx files containing the results of the REG-NBO7 analysis (most recent version of NBO)

INPUT: reg_inter, iqa_inter_header, SYS, path, g, sign=None, verbose=False
    reg_inter           = array of the NBO interactions with the REG and R² coefficients
    iqa_inter_header    = header (term names) of the NBO interactions
    SYS                 = name of the system considered
    path                = path to which the .csv file 
    g                   = bool:
                            True = grouped segments, 
                            False = ungrouped segments
    sign                = array of the sign of terms (optional)
    verbose             = False (default) : Boolean to ask the procedure to say everything (True) or be quiet (False)

OUTPUT: df_final, dataframe_list
    df_final = array of the X most important terms (in terms of REG, both positive and negative) of all segments
    dataframe_list = array of all terms on all segments ordered by REG coefficient

ERROR:
    (just think of creating the folder indicated by "path" before running the script, otherwise you get an error)

Added by Aël Cador for REGNBO7 (02-06-21)
###########################################################################################################
"""

def csvnbo7(reg_inter1, iqa_inter_header1, reg_inter2, iqa_inter_header2, SYS, path, g, sign=None, verbose=False):
    print("Starting generating .csv result files...")
    
    if sign != None:
        sign_inter1=sign[0]
        sign_inter2=sign[1]#[0] 
    dataframe_list=[]
    df_final = pd.DataFrame()
    for i in range(len(reg_inter1[0])):
        if sign != None:
            temp = [reg_inter1[0][i], reg_inter1[1][i], sign_inter1[i]]                 #inter1
            df_inter1 = pd.DataFrame(temp).transpose()                                #inter1
            df_inter1.columns = ["REG", "R2","SIGN"]                                  #inter1
            if verbose==True:
                print(df_inter1)
        else:                                     
            temp = [reg_inter1[0][i], reg_inter1[1][i]]                                #inter1
            if verbose==True:
                print("I",i,len(reg_inter1),len(reg_inter1[0]),len(reg_inter1[0][i]),len(iqa_inter_header1),len(iqa_inter_header1[i]))#[i]))
                print("TEMP",temp,reg_inter1[0][i], reg_inter1[1][i])
            df_inter1 = pd.DataFrame(temp).transpose()                                #inter1
            df_inter1.columns = ["REG", "R2"]                                         #inter1 NB: AC change: R to R2. This will actually be R until squared
        df_inter1.index = iqa_inter_header1[i]#[i]#[j] #0/i before                          #inter1
        df_inter1.to_csv(path +'/'+ SYS + "_NBO_seg_" +str(i+1) +".csv", sep=',')     #inter1 
        df_inter1 = df_inter1.rename_axis('TERM').reset_index()                        #inter1
           
        if sign != None:
            temp = [reg_inter2[0][i], reg_inter2[1][i], sign_inter2[i]]                 #inter1
            df_inter2 = pd.DataFrame(temp).transpose()                                #inter1
            df_inter2.columns = ["REG", "R2","SIGN"]                                  #inter1
            #print(df_inter)
        else:                                     
            temp = [reg_inter2[0][i], reg_inter2[1][i]]                                #inter2
            if verbose==True:
                print("I",i,len(reg_inter2),len(reg_inter2[0]),len(reg_inter2[0][i]),len(iqa_inter_header2),len(iqa_inter_header2[i]))#[i]))
                print("TEMP",temp,reg_inter2[0][i], reg_inter2[1][i])
            df_inter2 = pd.DataFrame(temp).transpose()                                #inter2
            df_inter2.columns = ["REG", "R2"]                                         #inter2 NB: AC change: R to R2. This will actually be R until squared
        df_inter2.index = iqa_inter_header2[i]#[i]#[j] #0/i before                          #inter2
        df_inter2.to_csv(path +'/'+ SYS + "_NBO_seg_" +str(i+1) +".csv", sep=',')     #inter2 
        df_inter2 = df_inter2.rename_axis('TERM').reset_index()                        #inter2
        
        df_temp = pd.concat([df_inter1,df_inter2]).sort_values('REG')
        
        df_temp['R2'] = df_temp['R2']**2 #AC change: R to R2.
        dataframe_list.append(df_temp)
        if verbose==True:
            print("DATA")
            print(dataframe_list)
        df_temp2=df_temp.copy()
        if verbose==True:
            print('df_temp2.shape',df_temp2.shape)
        for i in range(df_temp2.shape[0]): #AC change: R² filter (>0.8)
            if df_temp2.iat[i,1]<0.8: #second column is R²
                df_temp2.iat[i,1]=np.nan
        df_temp2=df_temp2.dropna(axis='rows',how='any')
        df_temp3=pd.concat([df_temp2[-5:], df_temp2[:5]], axis=0).sort_values('REG', ascending=False).reset_index()
        if verbose==True:
            print("df_temp3",df_temp3)        
        df_final=pd.concat([df_final, df_temp3], axis=1) 
    
    print(".csv result files successfully generated.")
    print("Starting generating .xlsx result files...")
    
    if g==False:   
         df_final.to_excel(path + '/' + SYS + '_NBO7_segment_analysis.xlsx')
         with pd.ExcelWriter(path + '/' + SYS + '_complete_NBO7_segment_analysis.xlsx') as writer: #writer object needed to have the lists written in different sheets
             for i in range(len(dataframe_list)):
                 dli=pd.DataFrame(dataframe_list[i]) #also need the data to be a dataframe to enable writing to XLSX file
                 dli.to_excel(writer, sheet_name='Segment'+str(i+1))
    else:
         df_final.to_excel(path + '/' + SYS + '_NBO7_segment_analysis_3seg.xlsx')
         with pd.ExcelWriter(path + '/' + SYS + '_complete_NBO7_segment_analysis_3seg.xlsx') as writer: #writer object needed to have the lists written in different sheets
             for i in range(len(dataframe_list)):
                 dli=pd.DataFrame(dataframe_list[i]) #also need the data to be a dataframe to enable writing to XLSX file
                 dli.to_excel(writer, sheet_name='Segment'+str(i+1))
    print(".xlsx result files successfully generated.")
    
    return dataframe_list


"""
###########################################################################################################
(Created a procedure to do the .csv and .xlsx files quickly instead of having it done repeatedly in the main script)
FUNCTION: excel
          generates .xlsx files containing the raw data of the REG analysis

INPUT: iqa_inter, iqa_inter_header, iqa_intra, iqa_intra_header, iqa_disp, iqa_disp_header, SYS, path, g
    iqa_inter           = array of the values of the interatomic terms
    iqa_inter_header    = header (term names) of the interatomic terms
    iqa_intra           = array of the values of the intra-atomic terms
    iqa_intra_header    = header (term names) of the intra-atomic terms
    iqa_disp            = array of the values of the dispersion terms
    iqa_disp_header     = header (term names) of the dispersion terms
    SYS                 = name of the system considered
    path                = path to which the .csv file 
    g                   = bool:
                            True = grouped values, 
                            False = ungrouped values

OUTPUT: .xlsx files of the values of the different terms

ERROR:
    (just think of creating the folder indicated by "path" before running the script, otherwise you get an error)

Added by Aël Cador (05-03-21)
###########################################################################################################
""" 

def excel(iqa_inter, iqa_inter_header, iqa_intra, iqa_intra_header, iqa_disp, iqa_disp_header, SYS, path, g):
    print("Starting generating .xlsx raw data files...")
    
    if g==False:
        filename=path + '/' + SYS + '_complete_term_value_list.xlsx'
    else: 
        filename=path + '/' + SYS + '_complete_term_value_list_grouped.xlsx'
    with pd.ExcelWriter(filename) as writer: #writer object needed to have the lists written in different sheets
        dli=pd.DataFrame(iqa_intra)#.transpose() #also need the data to be a dataframe to enable writing to XLSX file
        dli.index = iqa_intra_header
        dli.to_excel(writer, sheet_name='Intra')
        dli=pd.DataFrame(iqa_inter)#.transpose() #also need the data to be a dataframe to enable writing to XLSX file
        dli.index = iqa_inter_header
        dli.to_excel(writer, sheet_name='Inter')
        dli=pd.DataFrame(iqa_disp)#.transpose() #also need the data to be a dataframe to enable writing to XLSX file
        dli.index = iqa_disp_header
        dli.to_excel(writer, sheet_name='Disp')
        
    print(".xlsx raw data files generated.")

"""
###########################################################################################################
(Created a procedure to do the .csv and .xlsx files quickly instead of having it done repeatedly in the main script)
FUNCTION: excelnbo
          generates .xlsx files containing the raw data of the REG-NBO3 analysis

INPUT: nbo_E2,nbo_E2_header,SYS,path,g
    nbo_E2              = array of the values of the E2 terms
    nbo_E2_header       = header (term names) of the E2 terms
    SYS                 = name of the system considered
    path                = path to which the .csv file 
    g                   = bool:
                            True = grouped values, 
                            False = ungrouped values

OUTPUT: .xlsx files of the values of the different terms

ERROR:
    (just think of creating the folder indicated by "path" before running the script, otherwise you get an error)

Added by Aël Cador (May 2021)
###########################################################################################################
""" 

def excelnbo(nbo_E2,nbo_E2_header,SYS,path,g):
    print("Starting generating .xlsx raw data files...")
    
    if g==False:
        filename=path + '/' + SYS + '_complete_NBO_term_value_list.xlsx'
    else: 
        filename=path + '/' + SYS + '_complete_NBO_term_value_list_grouped.xlsx'
    with pd.ExcelWriter(filename) as writer: #writer object needed to have the lists written in different sheets
        dli=pd.DataFrame(nbo_E2)#.transpose() #also need the data to be a dataframe to enable writing to XLSX file
        dli.index = nbo_E2_header
        dli.to_excel(writer, sheet_name='E2')

    print(".xlsx raw data files generated.")
    
"""
###########################################################################################################
(Created a procedure to do the .csv and .xlsx files quickly instead of having it done repeatedly in the main script)
FUNCTION: excelnbo7
          generates .xlsx files containing the raw data of the REG-NBO7 analysis

INPUT: nbo_E2,nbo_E2_header,nbo_E4,nbo_E4_header,SYS,path,g
    nbo_E2              = array of the values of the E2 terms
    nbo_E2_header       = header (term names) of the E2 terms
    nbo_E4              = array of the values of the E4 terms
    nbo_E4_header       = header (term names) of the E4 terms
    SYS                 = name of the system considered
    path                = path to which the .csv file 
    g                   = bool:
                            True = grouped values, 
                            False = ungrouped values

OUTPUT: .xlsx files of the values of the different terms

ERROR:
    (just think of creating the folder indicated by "path" before running the script, otherwise you get an error)

Added by Aël Cador (May 2021)
###########################################################################################################
"""

def excelnbo7(nbo_E2,nbo_E2_header,nbo_E4,nbo_E4_header,SYS,path,g):
    print("Starting generating .xlsx raw data files...")
    
    if g==False:
        filename=path + '/' + SYS + '_complete_NBO7_term_value_list.xlsx'
    else: 
        filename=path + '/' + SYS + '_complete_NBO7_term_value_list_grouped.xlsx'
    with pd.ExcelWriter(filename) as writer: #writer object needed to have the lists written in different sheets
        for i in range(len(nbo_E2)):
            print(len(nbo_E2),len(nbo_E2[0]),len(nbo_E2_header),len(nbo_E2_header[0]))
            dli=pd.DataFrame(nbo_E2[i])#.transpose() #also need the data to be a dataframe to enable writing to XLSX file
            dli.index = nbo_E2_header[i]
            dli.to_excel(writer, sheet_name='E2-seg'+str(i+1))
            dli=pd.DataFrame(nbo_E4[i])#.transpose() #also need the data to be a dataframe to enable writing to XLSX file
            dli.index = nbo_E4_header[i]
            dli.to_excel(writer, sheet_name='E4-seg'+str(i+1))
    
    print(".xlsx raw data files generated.")