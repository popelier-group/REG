#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Jun 17 15:20:20 2019

@author: leo
"""
import pandas as pd

#files = ['INTER_AB_energies.ananke','INTRA_energies.ananke',
#         'Intra_KE_energies.ananke', 'Intra_VCL_energies.ananke', 
#         'Intra_VX_energies.ananke',  
#         'VCL_AB_energies.ananke', 'VCL_energies.ananke', 
#         'VX_AB_energies.ananke', 'VX_energies.ananke']

files=['total_energy.ananke']


for i in range(len(files)):
    df = pd.read_csv("./ananke_inputs/" + files[i], sep = '\t', header=None)
    del df[26]
    df=df.transpose()
    headers = ['total_energy']
    new_df  = pd.DataFrame(df.values[:,:], columns=headers)
    print(new_df.shape)
    new_df.to_csv(files[i]+".csv", sep="," , columns=headers)