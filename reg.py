"""
reg.py v0.0
L. J. Duarte, XXXXXX, P. L. A. Popelier 

Library with function to perform the IQA REG analysis. 
Check for updates at github.com/ljduarte
For details about the method, please see XXXXXXX

Please, report bugs and issues to leo.j.duarte@hotmail.com.br
coded by L. J. Duarte
"""



def regression(A, B, mode=None):
    """
    ###########################################################################################################
    FUNCTION: regression
              Perform a linear regression between A and B (B = slope*A + intercept)

    INPUT: A, B and mode
        A = X values
        B = Y values
        mode = None (default) : for regular linear regression
               "norm"         : use normalized values of A and B
               "std"          : use standardized values for A and B

    OUTPUT: [slope, intercept, person]
        slope     : angular coeficient of the linear regression
        intercept : linear coeficient of the linear regression
        pearson    : Pearson correlation coefficient

    ERROR:
        "Arrays must have the same size" : A and B have different sizes.
        "Mode not recognized. Use None, 'norm' or, 'std'" : invalide value for 'mode'
    ###########################################################################################################
    """
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



def find_critical(Y,X, min_points=3, use_inflex=False):
    """
    ###########################################################################################################
    FUNCTION: find_critical
              Takes the X and Y values of a function and return the critical points
              that fit the desirable criteria.

    INPUT: Y, X, min_points, use_inflex
        Y = ordinate values of a function
        C = abcissa values of a function
        min_point = 5 (default) : minimum amount of points between two critical points
        use_inflex = False (default) : Not search for inflexion points (second derivative = 0)
                     True            : Search for inflexion points (second derivative = 0)

    OUTPUT: critical_point_list
        critical_point_list = Array containing the index of critical points (zero indexed)

    ERROR:
        "Arrays must have the same size" : X and Y have different sizes.
        "Invalid value. Use True or False" : invalid value for use_inflex
        "Too many points between two critical points" : min_points must be lower than the X array length
    ###########################################################################################################
    """
    #ERRORS:  
    if len(Y) != len(X): ## Checks if X and Y have the same size
        raise ValueError("Arrays must have the same size") 
    if use_inflex != True and use_inflex != False: ## Checks if use_inlflex is valid
        raise ValueError("Invalid value. Use True or False")
    if min_points >= len(X):
        raise ValueError("Too many points between two critical points") ## Check if min_points is valid.
   
    #INTERNAL VARIABLES:
    critical_point = [] ## temporary array
    critical_point_list = [] ## final array
    
    Y_prime = [(Y[i+1]-Y[i])/(X[i+1]-X[i]) for i in range(len(Y)-1)] ## First order derivative  

    #FIND CRITICAL POINTS (FIRST DERIVATIVE = 0)    
    if use_inflex == False:
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

##    CHECKS MIN_POINTS CRITERIA
    while critical_point[0] < min_points-1 : critical_point.remove(critical_point[0])
    critical_point_list.append(critical_point[0])
    for i in range(1, len(critical_point)):
        if critical_point[i] - critical_point_list[-1] + 1 >= min_points:
           critical_point_list.append(critical_point[i])
       
    return critical_point_list



def split_segm(A, critical_point_list):
    """
    ###########################################################################################################
    FUNCTION: split_segm
              Takes the A array and divide it into N arrays according with the number of
              critical points. Each array corresponds to a segment of the REG analysis.

    INPUT: A, critical_point_list
        A = Array containing the abcissa values of a function
        critical_points_list = list of the index of the critical points (zero indexed)

    OUTPUT: segm
        segm = Array of arrays containing the values of A for each segment

    ERROR:
        "Invalid critical point index" : index of critical point out of range in the function array.
    ###########################################################################################################
    """
    critical_point_list = list(set(critical_point_list)) ## Remove duplicates and sort in crescent order
    #ERRORS:
    if critical_point_list[-1] >= len(A): ## Checks if index of critical point is in range of function
        raise ValueError("Invalid critical point index") 

    #INTERNAL VARIABLES:
    segm = [] ## array os arrays
    start = 0
    end = len(A)
    
    #DEVIDE FUNCTION INTO SEGMENTS:
    while start < sorted(critical_point_list)[-1]:
        for i in sorted(critical_point_list):
            segm.append([A[j] for j in range(start, i+1)])
            start = i
    segm.append([A[j] for j in range(start, end)])
    return segm

    
def reg(wfn_energy, control_coord, terms, critical=True, np=5, inflex=True, critical_index= [ ], mode="norm"):
    """
    ###########################################################################################################
    FUNCTION: reg
              Perform the REG analysis over all contributions inside "terms" array

    INPUT:  wfn_energy, control_coord, terms, critical, np, inflex, critical_index, mode
        wfn_energy = Array containing the energy values for each point (PES - Energy)
        control_coord = Array containing the control coodinate values for each point (PES - Coordinate)
        terms = Array of arrays, each one corresponding to one IQA contribution.
        critical = True(default): Search for critical points
                 = False        : user must provide the critical point index array
        np = 5(default):        :  minimum amount of points between two critical points (used if critical == True)
        use_inflex = False (default) : Not search for inflexion points (second derivative = 0) (used if critical == True)
                     True            : Search for inflexion points (second derivative = 0 (used if critical == True)
        critical_index = [](default) : list of critical points provided by the user. (used if critical == False)
        mode = None                    : for regular linear regression
               "norm"(default)         : use normalized values of A and B
               "std"                   : use standardized values for A and B

    OUTPUT: [REG_values/pearson_values][segment]
        REG_values = array of REG coefficients for each contribution. Each array correspond to a segment
        pearson_values = array of Pearson coefficients for each contribution. Each array correspond to a segment

    ERROR:
         "PES abcissa and ordenate must have same number of points" : wfn_energy and control_coord have differnt sizes
         "Contributions arrays must have same size": arrays inside terms have diferent sizes
    ###########################################################################################################
    """
    #ERRORS:
    if len(wfn_energy) != len(control_coord): # Checks if abcissa and ordenate of PES have same size
        raise ValueError("PES abcissa and ordenate must have same number of points")  
    for contribution in terms:
        if len(wfn_energy) != len(contribution):  # Checks if all contributions arrays have same size
            raise ValueError("Contributions arrays must have same size")
    
    #INTERNAL VARIABLES:      
    REG_values = []
    pearson_values = []
    split_terms = []  #split_terms[contribution][segment]
    split_wfn = []
    temp_REG = []
    temp_pearson = []

    #FIND CRITICAL POINTS:
    if critical == True: # Search for critical points automatically
        critical_points_list = find_critical(wfn_energy, control_coord, min_points=np, use_inflex=inflex)
    
    else:
        critical_points_list =  critical_index #User provides the critical points index

    #DIVIDE WFN ENERGY INTO SEGMENTS:
    split_wfn = split_segm(wfn_energy, critical_points_list)
    

    #DIVIDE TERMS INTO SEGMENTS:    
    for A in terms:
        split_terms.append(split_segm(A, critical_points_list))
      
    #LINEAR REGREESSION FOR CONTRIBUITON INSIDE EACH SEGMENT:
    for i in range(len(critical_points_list)+1):
        for A in split_terms:
            temp_REG.append(regression(split_wfn[i], A[i], mode=mode)[0])
            temp_pearson.append(regression(split_wfn[i], A[i], mode=mode)[2])
        REG_values.append(temp_REG)
        pearson_values.append(temp_pearson)
        temp_REG = []
        temp_pearson = []
 
    return REG_values, pearson_values



def integration_error(wfn_energies, IQA_energies):
    """
    ###########################################################################################################
    FUNCTION: integration_error
              calculates integration error for each PES point

    INPUT: wfn_energies, IQA_energies
        wfn_energies = list of total energies form wfn files
        IQA_energies = list of total energy obtained from IQA terms

    OUTPUT: [error, RMSE]
        error = list of wfn_energies - IQA_energies
        RMSE = Root mean square error

    ERROR:
        "Arrays must have the same size" : wfn_energies and IQA_energies have different sizes.

    ###########################################################################################################
    """
    #ERRORS:  
    if len(wfn_energies) != len(IQA_energies): ## Checks if X and Y have the same size
        raise ValueError("Arrays must have the same size")
    
    #CREATE LIST FOR EACH POINT ERROR    
    error = [wfn_energies[i] - IQA_energies[i] for i in range(len(wfn_energies))] # error values list
    #CALCULATE RMSE
    temp = [a**2 for a in error]
    RMSE = 2625.5*(sum(temp)/len(temp))**0.5
    return error, RMSE
    



def group_intra(intra_energies, intra_energies_header, groups):
    """
    ###########################################################################################################
    FUNCTION: group_intra
              group IQA intra-atomic terms into the user-defined groups

    INPUT: intra_energies, intra_energies_header, groups
        intra_energies: list of IQA intra-atomic terms
        intra_energies_header: header of IQA intra atomic temrs.
        groups: list of list of atomic labels, each internal list correspond to a diferent group

    OUTPUT: [new_terms, new_header]
        new_terms: list of the grouped IQA terms values
        new_header: list of the grouped terms header

    ERROR:

    ###########################################################################################################
    """
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
                    
    return new_terms, new_header



def group_inter(inter_energies, inter_energies_header, groups):
    """
    ###########################################################################################################
    FUNCTION: group_inter
              group IQA interatomic terms into the user-defined groups

    INPUT: inter_energies, inter_energies_header, groups
        inter_energies: list of IQA interatomic terms
        inter_energies_header: header of IQA interatomic temrs.
        groups: list of list of atomic labels, each internal list correspond to a diferent group

    OUTPUT: [new_terms, new_header]
        new_terms: list of the grouped IQA terms values
        new_header: list of the grouped terms header

    ERROR:

    ###########################################################################################################
    """
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

    #SUM IQA CONTRIBUTIONS OF SAME LABEL
    for term in temp_header:
        temp =0
        if term not in temp_header2:
            temp_header2.append(term)
            for i in range(len(temp_header)):
                if term == temp_header[i]:
                    temp = temp + inter_energies[i]
            temp_terms.append(temp)

    # ORGANIZE AND REMOVE REDUNDANCIES            
   
    for i in range(len(temp_header2)):
        temp = 0
        if temp_header2[i] not in new_header:
            if temp_header2[i] not in redundancies:
                new_header.append(temp_header2[i])        
                temp = temp + temp_terms[i]
                group_label_1 = temp_header2[i][-14:][temp_header2[i][-14:].find('-')+1:temp_header2[i][-14:].find('_')]
                group_label_2 = temp_header2[i][-14:][temp_header2[i][-14:].find('_')+1:]              
                for j in range(len(temp_header2)):
                    if temp_header2[i][:-1*temp_header2[i].find('(')] == temp_header2[j][:-1*temp_header2[j].find('(')]:
                        if temp_header2[j] not in new_header:
                            if temp_header2[j] not in redundancies:
                                if group_label_1 in temp_header2[j] and group_label_2 in temp_header2[j] and group_label_2 != group_label_1:
                                    temp = temp + temp_terms[j]
                                    redundancies.append(temp_header2[j])                   
                new_terms.append(temp)

    return  new_terms, new_header
