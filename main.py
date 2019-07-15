# Import modules
import reg 
import aimall_utils
import numpy as np
import pandas as pd
#Define paths to the files:
path = "./succinic_scan/"
wfn_files = [path + str('{:03d}'.format(i)) + '/' + str('{:03d}'.format(i)) + '.wfn' for i in range(0, 22)]

folders = [path + str('{:03d}'.format(i)) + '/' + str('{:03d}'.format(i)) +'_atomicfiles' for i in range(0, 22)]

#path = './monomer/conformationalAnalysisMonomer/'
#wfn_files = [path + 'dihedral' + str('{:03d}'.format(i)) + '_B3LYP.wfn' for i in range(10, 360, 10)]
#folders = [path + 'dihedral' + str('{:03d}'.format(i)) + '_B3LYP_atomicfiles' for i in range(10, 360, 10)]


groups = [['c10','o3','o15','h21'], #1
           ['o9','c2','o4','h5'], #2
           ['c14','h18','h20','c19','h25','h24'], #3
           ['c1','h6','h8','c7','h11','h12'], #4
           ['c23','o27','o26','h28'], #5
           ['c13','o16','o17','h22']] #6

#groups =[[10,3,15,21], [9,2,4,5],
#          [14,18,20,19,25,24], [1,6,8,7,11,12],
#          [23,27,26,28], [13,16,17,22]]


#Define the IQA terms
intra_prop = ['E_IQA_Intra(A)']
inter_prop = ['VC_IQA(A,B)', 'VX_IQA(A,B)']

#control_coordinate
cc = [i for i in range (180, 70, -5)]

#get wfn energies:
total_energy_wfn = np.array(aimall_utils.get_aimall_wfn_energies(wfn_files))

#get atomslist:
atoms = aimall_utils.get_atom_list(wfn_files[0])

#get intra-atomic IQA terms:
iqa_intra = np.array(aimall_utils.intra_property_from_int_file(folders,intra_prop, atoms)[0])
iqa_intra_header = np.array(aimall_utils.intra_property_from_int_file(folders,intra_prop, atoms)[1])
#get inter-atomic IQA terms:
iqa_inter = np.array(aimall_utils.inter_property_from_int_file(folders, inter_prop, atoms)[0])
iqa_inter_header = aimall_utils.inter_property_from_int_file(folders,inter_prop, atoms)[1]

total_energy_iqa = sum(iqa_inter) + sum(iqa_intra)

#integration error:
rmse_integration = reg.integration_error(total_energy_wfn, total_energy_iqa)[1]

#group intra

reg_intra = reg.reg(total_energy_wfn, cc, iqa_intra)
reg_inter = reg.reg(total_energy_wfn, cc, iqa_inter)

#organize and save 
for i in range(len(reg_inter[0])):
    temp = [reg_inter[0][i], reg_inter[1][i]]
    df = pd.DataFrame(temp).transpose()
    df.columns = ["REG", "R"]
    df.index = iqa_inter_header
    df.to_csv("Inter seg_" +str(i+1) +".csv", sep=',')

for i in range(len(reg_intra[0])):
    temp = [reg_intra[0][i], reg_intra[1][i]]
    df = pd.DataFrame(temp).transpose()
    df.columns = ["REG", "R"]
    df.index = iqa_intra_header
    df.to_csv("Intra seg_" +str(i+1) + ".csv", sep=',')
    
