#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun May 12 22:34:27 2024

@author: leo
"""

import numpy as np
import scan_parameterspace_funcs as fcs
import scan_SPheno_funcs as SPfcs
import pandas as pd
import matplotlib.pyplot as plt

N_points = int(2e3)

tanb_min=fcs.tanb_min
tanb_max=fcs.tanb_max

def find_random_points_coup(N):
    '''Generate N random points in the parameter space. mA, mH, mHpm varying from 125 to 1000, cos(alpha) from -1 to 1, tan(beta) from 1.5 to 5 and M from 1e3 to 1e7.
        Returns: mA, mH, mHpm, cosa, tanb, M.
    '''    
    lmax = 4*np.pi
    
    DataFrame = pd.DataFrame(columns = ['l1', 'l2', 'l3', 'l4', 'l5', 'm12', 'tanb'])
        
    numb = np.random.rand(7,N)
        
    DataFrame['l1'] = -lmax+numb[0]*(2*lmax)
    DataFrame['l2'] = -lmax+numb[1]*(2*lmax)
    DataFrame['l3'] = -lmax+numb[2]*(2*lmax)
    DataFrame['l4'] = -lmax+numb[3]*(2*lmax)
    DataFrame['l5'] = -lmax+numb[4]*(2*lmax)
    DataFrame['m12'] = -(1e3+10**(numb[5]*(6-3)))
    DataFrame['tanb'] = tanb_min+numb[6]*(tanb_max-tanb_min)
        
    return DataFrame

SPheno_input = find_random_points_coup(N_points)

#outp = open("Data_output.dat",'w')
outpt=pd.DataFrame()

for index, row in SPheno_input.iterrows():
    SPfcs.write_spheno_LHA(list(row))
    success = SPfcs.execute_spheno()
    if success:
        outpt = pd.concat([outpt,SPfcs.read_spheno_obs()])
        print(index)
    # else:
    #     proxrow = np.copy(row)
    #     proxrow[5]=-proxrow[5]
    #     SPfcs.write_spheno_LHA(list(proxrow))
    #     success = SPfcs.execute_spheno()
    #     if success:
    #         outpt = pd.concat([outpt,SPfcs.read_spheno_obs()])
    #         print(index)

outpt = outpt.drop_duplicates()
            
#%% 

def plotter_2(param1,param2):

    fig, ax = plt.subplots(figsize=(10, 10))
    
    plt.scatter(outpt[param1],outpt[param2])
    #plt.xlabel(str_to_tex(param1), size=25)
    plt.xticks(size=20)
    #plt.ylabel(str_to_tex(param2), size=25)
    plt.yticks(size=20)
    ax.grid()
    
    plt.show()

    return 0