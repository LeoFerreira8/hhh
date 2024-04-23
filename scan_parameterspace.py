#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Apr 22 11:09:52 2024

@author: leonardoferreira
"""

from anyBSM import anyBSM
import numpy as np
import matplotlib.pyplot as plt
import scan_parameterspace_funcs as fcs
import pandas as pd

def perturbativity_bounds(l):
    '''Returns True if coupling l is above the perturbativity bound, and False otherwise.'''
    lmax = 4*np.pi
    res = np.where(np.abs(l)<lmax,False,True)
    
    return res

def stability_bounds(l_a):
    '''Returns True if array of couplings l_a does not comply with the stability bounds, and False otherwise.'''
    res = np.where(
        (l_a[0] > 0) & 
        (l_a[1] > 0) & 
        (-np.sqrt(l_a[0]*l_a[1]) < l_a[2]) & 
        (l_a[2]+l_a[3]-np.abs(l_a[4]) > -np.sqrt(l_a[0]*l_a[1]))
        ,False,True)
    
    return res
    
#%%

N_points = int(1e5)

Table = fcs.find_random_points(N_points)

sb = np.sin(fcs.beta(Table['tanb']))
cb = np.cos(fcs.beta(Table['tanb']))
sa = np.sin(fcs.alpha(Table['cosa']))
ca = Table['cosa']

mh = fcs.mhSM
v = fcs.vSM
lambda6 = 0
lambda7 = 0

m122 = fcs.m122M(Table['M'], sb, cb)
l1 = fcs.lamb1(Table['mH'], mh, m122, ca, sa, sb, cb, v, lambda6, lambda7)
l2 = fcs.lamb2(Table['mH'], mh, m122, ca, sa, sb, cb, v, lambda6, lambda7)
l3 = fcs.lamb3(Table['mH'], mh, Table['mHpm'], m122, ca, sa, sb, cb, v, lambda6, lambda7)
l4 = fcs.lamb4(Table['mA'], Table['mHpm'], m122, sb, cb, v, lambda6, lambda7)
l5 = fcs.lamb5(Table['mA'], m122, sb, cb, v, lambda6, lambda7)

Table1 = pd.DataFrame(np.array([m122,l1,l2,l3,l4,l5]).T,columns=['m122','l1','l2','l3','l4','l5'])

Table = pd.concat([Table,Table1],axis=1)

#%%                                 Analysis

###     Perturbativity tests

TableP = Table

for i in range(1,6):
    cnd = perturbativity_bounds(TableP['l'+str(i)])
    TableP = TableP.drop(TableP[cnd].index)

###     Stability tests

cnd = stability_bounds([Table['l'+str(i)] for i in range(1,6)])
TableStab = Table.drop(Table[cnd].index)

###     Testing both

cnd = stability_bounds([TableP['l'+str(i)] for i in range(1,6)])
TableTot = TableP.drop(TableP[cnd].index)

#THDM2 = anyBSM('THDMII', scheme_name = 'OS')

#%%                         Calculate lambda

lamb = []
lambtree = []

for i in TableTot.index:
    sino = np.sin(np.arccos(TableTot.at[i,'cosa'])-np.arctan(TableTot.at[i,'tanb']))
    if sino == 1:
        THDM2 = anyBSM('THDMII', scheme_name = 'OSalignment')
    else:
        THDM2 = anyBSM('THDMII', scheme_name = 'OS')
        
    THDM2.setparameters({'Mh2': TableTot.at[i,'mH'], 'MAh2': TableTot.at[i,'mA'], 'MHm2': TableTot.at[i,'mHpm'], 'TanBeta': TableTot.at[i,'tanb'], 'SinBmA': sino}) #Define new mass in anyBSM
    dic = THDM2.lambdahhh()
    lamb.append(-np.real(dic['total'])/fcs.Gammahhh_treelevel(0, 0))  #Recalculate lambda
    lambtree.append(-np.real(dic['treelevel'])/fcs.Gammahhh_treelevel(0, 0))  #Recalculate lambda
    
TableTot = pd.concat([TableTot.reset_index(),pd.DataFrame({'kappa': lamb, 'kappa-tree': lambtree})],axis=1)

#%%                                 Plots

def str_to_tex(strg):
    latex_parameters = [r'$m_A$ [GeV]',r'$m_H$ [GeV]',r'$m_{H^\pm}$ [GeV]',r'$\cos{\alpha}$',r'$\tan{\beta}$',r'$M$ [GeV]',r'$m_{12}^2$ [GeV$^2$]', r'$\lambda_1$', r'$\lambda_2$', r'$\lambda_3$', r'$\lambda_4$', r'$\lambda_5$',r'$\kappa_\lambda$',r'$\kappa_\lambda^{(0)}$']
    
    if strg=='mA':
        return latex_parameters[0]
    elif strg=='mH':
        return latex_parameters[1]
    elif strg=='mHpm':
        return latex_parameters[2]
    elif strg=='cosa':
        return latex_parameters[3]
    elif strg=='tanb':
        return latex_parameters[4]
    elif strg=='M':
        return latex_parameters[5]
    elif strg=='m122':
        return latex_parameters[6]
    elif strg=='l1':
        return latex_parameters[7]
    elif strg=='l2':
        return latex_parameters[8]
    elif strg=='l3':
        return latex_parameters[9]
    elif strg=='l4':
        return latex_parameters[10]
    elif strg=='l5':
        return latex_parameters[11]
    elif strg=='kappa':
        return latex_parameters[12]
    elif strg=='kappa-tree':
        return latex_parameters[13]
    else:
        raise ValueError("Invalid parameter.")
    
def plotter_2(param1,param2):

    fig, ax = plt.subplots(figsize=(10, 10))
    
    plt.scatter(TableTot[param1],TableTot[param2])
    plt.xlabel(str_to_tex(param1), size=25)
    plt.xticks(size=20)
    plt.ylabel(str_to_tex(param2), size=25)
    plt.yticks(size=20)
    ax.grid()
    
    plt.show()

    return 0