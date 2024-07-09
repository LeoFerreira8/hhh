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
import gc
import scan_SPheno_funcs as SPfcs
import scan_higgs_tools_funcs as hggfcs
import bsg
from scipy import stats
from multiprocessing import Pool

N_parameters = 6-1*fcs.alignment
THDM_type = fcs.THDM_type
if THDM_type=='':
    latex_model = '2HDM - Type I'
    THDM_type='THDMI'
elif THDM_type=='II' and not fcs.small_l5:
    latex_model = '2HDM - Type II'
    THDM_type='THDMII'
elif THDM_type=='II' and fcs.small_l5:
    latex_model = '2HDM - 331 EFT'
    THDM_type='THDMII'
    
if fcs.alignment:
    latex_alignment = r'($\beta-\alpha=\frac{\pi}{2}$)'
else:
    latex_alignment = r'($\beta-\alpha\in[\frac{\pi}{2}\pm%.2f]$)' %(fcs.non_alignment_max,)

def perturbativity_bounds(l):
    '''Returns True if coupling l is above the perturbativity bound, and False otherwise.'''
    lmax = np.sqrt(4*np.pi)
    res = np.where(np.abs(l)<lmax,False,True)
    
    return res

def stability_bounds(l_a):
    '''Returns True if array of couplings l_a does not comply with the stability bounds, and False otherwise.'''
    res = np.where(
        (l_a[0] > 0) 
        & (l_a[1] > 0) 
        & (-np.sqrt(l_a[0]*l_a[1]) < l_a[2]) 
        & (l_a[2]+l_a[3]-np.abs(l_a[4]) > -np.sqrt(l_a[0]*l_a[1]))
        ,False,True)
    
    return res

def Spheno_Higgstools_calc_(SPheno_input):
    outpt=pd.DataFrame()
    for index, row in SPheno_input.iterrows():
        SPfcs.write_spheno_LHA(list(row))
        success = SPfcs.execute_spheno()
        if success:
            hggs_SPhen_outpt = pd.concat([SPfcs.read_spheno_obs(),hggfcs.Higgs_tools_scan()],axis=1)
            outpt = pd.concat([outpt,hggs_SPhen_outpt])
        else:
            print('SPheno could not calculate %d' %index)
            proxrow = np.copy(row)
            proxrow[5]=-proxrow[5]
            SPfcs.write_spheno_LHA(list(proxrow))
            success = SPfcs.execute_spheno()
            if success:
                outpt = pd.concat([outpt,SPfcs.read_spheno_obs()])
                print('Calculation of %d successful.' %index)
            else:
                print('Failed')
            
    outpt = outpt.drop_duplicates()
    
    return outpt

def STU_constraint(s,t,u):
    '''Returns True if STU complies with the bounds, and False otherwise.
    ΔS = 0.00+-0.07 
    ΔT = 0.05+-0.06 
    
    ΔS = -0.01+-0.07 
    ΔT = 0.04+-0.06 
    '''
    smax = 0.07
    smin = -0.07
    tmax = 0.11
    tmin = -0.01
    umax = 0.00
    umin = -0.00
    
    res = np.where((smin<s)
                   & (s<smax)
                   & (tmin<t) 
                   & (t<tmax)
                   # & (umin<u)
                   # & (u<umax),
                   ,False,True)
    
    return res

def collider_const(HiggsB):
    return np.where(HiggsB,False,True)

def signals_const(chisq):
    N_observables = hggfcs.signals.observableCount()
    N_d_freedom = N_observables-N_parameters
    
    Δχ = chisq-np.min(chisq)
    CL = 1-stats.chi2.cdf(Δχ, N_d_freedom)
    
    return np.where(CL>0.95,False,True)
    
#%%

N_points = int(5e4)

Table = fcs.find_random_points(N_points)

sb = np.sin(fcs.beta(Table['tanb']))
cb = np.cos(fcs.beta(Table['tanb']))
sa = np.sin(fcs.alpha(Table['cosa']))
ca = Table['cosa']

mh = fcs.mhSM
v = fcs.v
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
Table1 = None

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
TableTot.reset_index()

Table = None
TableP = None
TableStab = None
l1 = None
l2 = None
l3 = None
l4 = None
l5 = None
sa = None
ca = None
sb = None
cb = None

gc.collect()

#%%                         Calculate lambda

lamb = []
lambtree = []

sino = np.sin(fcs.beta(TableTot['tanb'])-fcs.alpha(TableTot['cosa']))
if fcs.alignment:
    THDM2 = anyBSM(THDM_type, scheme_name = 'OSalignment')
else:
    THDM2 = anyBSM(THDM_type, scheme_name = 'OS')

for i in TableTot.index:
    # sino = np.sin(fcs.beta(TableTot.at[i,'tanb'])-fcs.alpha(TableTot.at[i,'cosa']))
    # if sino == 1:
    #     THDM2 = anyBSM('THDMI', scheme_name = 'OSalignment')
    # else:
    #     THDM2 = anyBSM('THDMI', scheme_name = 'OS')
    if not fcs.alignment and sino[i]==1.0:
        THDM2.load_renormalization_scheme('OSalignment')
        
    THDM2.setparameters({'Mh2': TableTot.at[i,'mH'], 'MAh2': TableTot.at[i,'mA'], 'MHm2': TableTot.at[i,'mHpm'], 'TanBeta': TableTot.at[i,'tanb'], 'SinBmA': sino[i],'M': TableTot.at[i,'M']}) #Define new mass in anyBSM
    dic = THDM2.lambdahhh()
    lamb.append(-np.real(dic['total'])/fcs.Gammahhh_treelevel(0, 0))  #Recalculate lambda
    lambtree.append(-np.real(dic['treelevel'])/fcs.Gammahhh_treelevel(0, 0))  #Recalculate lambda
    
sino = np.sin(fcs.beta(TableTot['tanb'])-fcs.alpha(TableTot['cosa']))
# Using that sinBmA ~ 1-x²/2
kappa_kan_x = fcs.Gammahhh_oneloop(np.sqrt(2*(1-sino)), TableTot['M'], TableTot['mH'], TableTot['mA'], TableTot['mHpm'])/fcs.Gammahhh_treelevel(0, 0)
kappa_kan = fcs.Gammahhh_oneloop_cos(fcs.beta(TableTot['tanb'])-fcs.alpha(TableTot['cosa']), fcs.beta(TableTot['tanb']), TableTot['M'], TableTot['mH'], TableTot['mA'], TableTot['mHpm'])/fcs.Gammahhh_treelevel(0, 0)
    
TableTot = pd.concat([TableTot,pd.DataFrame({'kappa': lamb, 'kappa-tree': lambtree, 'kappa-kan-x': kappa_kan_x, 'kappa-kan': kappa_kan})],axis=1)

#%%                                 Calculate S,T,U & collider constraints

Sp_in = TableTot.T.loc[['l1','l2','l3','l4','l5','m122','tanb']].T
Sp_in['m122'] = -Sp_in['m122'] # Different convention from SPheno.
Sp_in['l1'] = Sp_in['l1']/2 # Different convention from SPheno.
Sp_in['l2'] = Sp_in['l2']/2 # Different convention from SPheno.
Sp_in = Sp_in.rename({'m122': 'm12'},axis=1) # Fixed name with SPheno.

TotalSP = Spheno_Higgstools_calc_(Sp_in)

#%%                                 Impose bounds from STU and Collider

STU = TotalSP.T.loc[['S-parameter (1-loop BSM)','T-parameter (1-loop BSM)','U-parameter (1-loop BSM)','HiggsB','HiggsS']].T
STU=STU.drop(index=0).set_index(TableTot.index)
TableTot_STU = pd.concat([TableTot,STU],axis=1)

cnd = STU_constraint(TableTot_STU['S-parameter (1-loop BSM)'],TableTot_STU['T-parameter (1-loop BSM)'],TableTot_STU['U-parameter (1-loop BSM)'])
TableTot_STU = TableTot_STU.drop(TableTot_STU[cnd].index)

cnd = collider_const(TableTot_STU['HiggsB'])
TableTot_STU_Collid = TableTot_STU.drop(TableTot_STU[cnd].index)

#%%                                 Impose bounds from HiggsSignals

cnd = signals_const(np.array(TableTot_STU_Collid['HiggsS'],dtype=float))
TableTot_STU_Collid = TableTot_STU_Collid.drop(TableTot_STU_Collid[cnd].index)

#%%                                 Impose bounds from BSG

cnd = bsg.Constraints_BSG(np.array(TableTot_STU_Collid['tanb'],dtype=float), np.array(TableTot_STU_Collid['mHpm'],dtype=float))
TableTot_STU_Collid_BSG = TableTot_STU_Collid.drop(TableTot_STU_Collid[cnd].index)

#%%                             Save data

Saves=False
if fcs.alignment:
    strga = 'A'
else:
    strg = str(fcs.non_alignment_max)
    
if fcs.small_l5:
    strgl5 = '331lk'
else:
    strgl5 = ''

if Saves:
    TableTot.to_csv('THDM'+fcs.THDM_type+strgl5+'-'+strga+'-Theo.csv')
    TableTot_STU.to_csv('THDM'+fcs.THDM_type+strgl5+'-'+strga+'-STU.csv')
    TableTot_STU_Collid.to_csv('THDM'+fcs.THDM_type+strgl5+'-'+strga+'-Collid.csv')

#%%                                 Plots

def str_to_tex(strg):
    latex_parameters = [r'$m_A$ [GeV]',r'$m_H$ [GeV]',r'$m_{H^\pm}$ [GeV]',r'$\cos{\alpha}$',r'$\tan{\beta}$',r'$M$ [GeV]',r'$m_{12}^2$ [GeV$^2$]', r'$\lambda_1$', r'$\lambda_2$', r'$\lambda_3$', r'$\lambda_4$', r'$\lambda_5$',r'$\kappa_\lambda$',r'$\kappa_\lambda^{(0)}$',r'$\kappa_\lambda^{\text{Kan. Aprox}}$',r'$\kappa_\lambda^{\text{Kan}}$',r'$\sin{(\beta-\alpha)}$']
    
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
    elif strg=='kappa-kan-x':
        return latex_parameters[14]
    elif strg=='kappa-kan':
        return latex_parameters[15]
    elif strg=='sino':
        return latex_parameters[16]
    else:
        raise ValueError("Invalid parameter.")
    
def plotter(param1,param2,*dataset):

    fig, ax = plt.subplots(figsize=(10, 10))
    
    plt.title(latex_model +' '+ latex_alignment,size=20)
    
    for tbl in dataset:
        if param1=='sino':
            sino = np.sin(fcs.beta(tbl['tanb'])-fcs.alpha(tbl['cosa']))
            plt.scatter(sino,tbl[param2])
        elif param2=='sino':
            sino = np.sin(fcs.beta(tbl['tanb'])-fcs.alpha(tbl['cosa']))
            plt.scatter(tbl[param1],sino)
        else:
            plt.scatter(tbl[param1],tbl[param2])
    plt.xlabel(str_to_tex(param1), size=25)
    plt.xticks(size=20)
    #plt.xlim(125,900)
    plt.ylabel(str_to_tex(param2), size=25)
    plt.yticks(size=20)
    #plt.ylim(125,900)
    ax.grid()
    
    plt.show()

    return 0

def plotter_3(param1,param2,param3,*dataset):

    fig, ax = plt.subplots(figsize=(10, 10))
    
    plt.title(latex_model +' '+ latex_alignment,size=20)
    
    for tbl in dataset:
        if param1=='sino':
            sino = np.sin(fcs.beta(tbl['tanb'])-fcs.alpha(tbl['cosa']))
            plt.scatter(sino,tbl[param2],c=tbl[param3])
        elif param2=='sino':
            sino = np.sin(fcs.beta(tbl['tanb'])-fcs.alpha(tbl['cosa']))
            plt.scatter(tbl[param1],sino,c=tbl[param3])
        elif param3=='sino':
            sino = np.sin(fcs.beta(tbl['tanb'])-fcs.alpha(tbl['cosa']))
            plt.scatter(tbl[param1],tbl[param2],c=sino)
        else:
            plt.scatter(tbl[param1],tbl[param2],c=tbl[param3])
        
    cbar=plt.colorbar(label=str_to_tex(param3))
    plt.xlabel(str_to_tex(param1), size=25)
    plt.xticks(size=20)
    #plt.xlim(125,900)
    plt.ylabel(str_to_tex(param2), size=25)
    plt.yticks(size=20)
    #plt.ylim(125,900)
    ax.grid()
    cbar.ax.tick_params(labelsize=20)  # Change size of tick labels
    cbar.ax.yaxis.label.set_size(20)   # Change size of colorbar label
    
    plt.show()

    return 0

def plotter_comp(param1,param2,param3,param4,*dataset):    

    fig, ax = plt.subplots(figsize=(10, 10))
    
    plt.title(latex_model +' '+ latex_alignment,size=20)
    
    for tbl in dataset:
        ratio = tbl[param3]/tbl[param4]
        if param1=='sino':
            sino = np.sin(fcs.beta(tbl['tanb'])-fcs.alpha(tbl['cosa']))
            plt.scatter(sino,tbl[param2],c=ratio)
        elif param2=='sino':
            sino = np.sin(fcs.beta(tbl['tanb'])-fcs.alpha(tbl['cosa']))
            plt.scatter(tbl[param1],sino,c=ratio)
        else:
            plt.scatter(tbl[param1],tbl[param2],c=ratio)
    cbar=plt.colorbar(label=str_to_tex(param3)+'/'+str_to_tex(param4))
    plt.xlabel(str_to_tex(param1), size=25)
    plt.xticks(size=20)
    #plt.xlim(125,900)
    plt.ylabel(str_to_tex(param2), size=25)
    plt.yticks(size=20)
    #plt.ylim(125,900)
    ax.grid()
    cbar.ax.tick_params(labelsize=20)  # Change size of tick labels
    cbar.ax.yaxis.label.set_size(20)   # Change size of colorbar label
    
    plt.show()

    return 0

def plotter_diff(param1,param2,param3,param4,*dataset):

    fig, ax = plt.subplots(figsize=(10, 10))
    
    plt.title(latex_model +' '+ latex_alignment,size=20)
    
    for tbl in dataset:
        plt.scatter(tbl[param1]-tbl[param2],tbl[param3]-tbl[param4])
        
    plt.xlabel(str_to_tex(param1).replace(' [GeV]','')+'-'+str_to_tex(param2), size=25)
    plt.xticks(size=20)
    plt.xlim(-200,200)
    plt.ylabel(str_to_tex(param3).replace(' [GeV]','')+'-'+str_to_tex(param4), size=25)
    plt.yticks(size=20)
    plt.ylim(-200,200)
    ax.grid()
    
    plt.show()

    return 0

def plotter_diff_color(param1,param2,param3,param4,param5,*dataset):

    fig, ax = plt.subplots(figsize=(10, 10))
    
    plt.title(latex_model +' '+ latex_alignment,size=20)
    
    for tbl in dataset:
        if param5=='sino':
            sino = np.sin(fcs.beta(tbl['tanb'])-fcs.alpha(tbl['cosa']))
            plt.scatter(tbl[param1]-tbl[param2],tbl[param3]-tbl[param4],c=sino)
        else:
            plt.scatter(tbl[param1]-tbl[param2],tbl[param3]-tbl[param4],c=tbl[param5])
    
    cbar=plt.colorbar(label=str_to_tex(param5))
    plt.xlabel(str_to_tex(param1).replace(' [GeV]','')+'-'+str_to_tex(param2), size=25)
    plt.xticks(size=20)
    plt.xlim(-450,200)
    plt.ylabel(str_to_tex(param3).replace(' [GeV]','')+'-'+str_to_tex(param4), size=25)
    plt.yticks(size=20)
    plt.ylim(-320,350)
    ax.grid()
    cbar.ax.tick_params(labelsize=20)  # Change size of tick labels
    cbar.ax.yaxis.label.set_size(20)   # Change size of colorbar label
    
    plt.show()

    return 0

def plotter_diff_comp(param1,param2,param3,param4,param5,param6,*dataset):    

    fig, ax = plt.subplots(figsize=(10, 10))
    
    plt.title(latex_model +' '+ latex_alignment,size=20)
    
    for tbl in dataset:
        ratio = tbl[param5]/tbl[param6]
        plt.scatter(tbl[param1]-tbl[param2],tbl[param3]-tbl[param4],c=ratio)
    cbar=plt.colorbar(label=str_to_tex(param5)+'/'+str_to_tex(param6))
    plt.xlabel(str_to_tex(param1).replace(' [GeV]','')+'-'+str_to_tex(param2), size=25)
    plt.xticks(size=20)
    plt.xlim(-450,200)
    plt.ylabel(str_to_tex(param3).replace(' [GeV]','')+'-'+str_to_tex(param4), size=25)
    plt.yticks(size=20)
    plt.ylim(-320,350)
    ax.grid()
    cbar.ax.tick_params(labelsize=20)  # Change size of tick labels
    cbar.ax.yaxis.label.set_size(20)   # Change size of colorbar label
    
    plt.show()

    return 0