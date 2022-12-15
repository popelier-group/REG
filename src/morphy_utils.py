
def read_morphy_outputs(file_list):
    values = []
    for file in file_list:
           
        f = open(file, 'r').readlines()
        interactions_temp = []
        values_temp = []
        for i in range(len(f)):
            if 'Total number of interaction'in f[i]:
                start_line = i+1
    
        for i in range(start_line, len(f), 5):
            if f[i].split()[0] == 'interaction':
                interactions_temp.append('Vcorr(' + f[i].split()[4] + f[i].split()[2] + ',' + f[i].split()[3] + f[i].split()[1] +')')
                if f[i].split()[4] + f[i].split()[2] != f[i].split()[3] + f[i].split()[1]:
                    values_temp.append(float(f[i+1].split()[1]))
                else:
                    values_temp.append(float(f[i+1].split()[1])/2)
        values.append(values_temp)
    
    final_values = []
    #Organize the data
    for i in range(len(values[0])):
        final_values.append([values[j][i] for j in range(len(values))])
        
    return interactions_temp, final_values