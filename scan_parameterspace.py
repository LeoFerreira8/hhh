#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Apr 22 11:09:52 2024

@author: leonardoferreira
"""

from anyBSM import anyBSM
import numpy as np
import scan_parameterspace_funcs as fcs
import pandas as pd
import gc
import scan_SPheno_funcs as SPfcs
import scan_higgs_tools_funcs as hggfcs
import bsg
import quartic_couplings as qtcp
from scipy import stats
from time import time

s_mean = 0.00
δs = 0.07

t_mean = 0.05
δt = 0.06

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
    lmax = 4*np.pi
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
    smax = s_mean+δs
    smin = s_mean-δs
    tmax = t_mean+δt
    tmin = t_mean-δt
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

def calculate_lambda(DtFrame):
    lamb = []
    lambtree = []

    sino = np.sin(fcs.beta(DtFrame['tanb'])-fcs.alpha(DtFrame['cosa']))
    
    if fcs.alignment:
        THDM2 = anyBSM(THDM_type, scheme_name = 'OSalignment')
    else:
        THDM2 = anyBSM(THDM_type, scheme_name = 'OS')

    for i in DtFrame.index:
        if not fcs.alignment and sino[i]==1.0:
            THDM2.load_renormalization_scheme('OSalignment')
            
        THDM2.setparameters({'Mh2': DtFrame.at[i,'mH'], 'MAh2': DtFrame.at[i,'mA'], 'MHm2': DtFrame.at[i,'mHpm'], 'TanBeta': DtFrame.at[i,'tanb'], 'SinBmA': sino[i],'M': DtFrame.at[i,'M'], 'MWm': fcs.MW, 'MWp': fcs.MW}) #Define new mass in anyBSM
        THDM2.progress=False
        THDM2.warnSSSS=False
        dic = THDM2.lambdahhh()
        lamb.append(-np.real(dic['total'])/fcs.Gammahhh_treelevel(0, 0))  #Recalculate lambda
        lambtree.append(-np.real(dic['treelevel'])/fcs.Gammahhh_treelevel(0, 0))  #Recalculate lambda
        
    return lamb, lambtree

def calculate_quartics(DtFrame):
    c93 = qtcp.C_93(fcs.alpha(DtFrame['cosa']), fcs.beta(DtFrame['tanb']), DtFrame['mH'], DtFrame['l5'])
    c94 = qtcp.C_94(fcs.alpha(DtFrame['cosa']), fcs.beta(DtFrame['tanb']), DtFrame['mH'], DtFrame['mA'], DtFrame['l5'])
    c102 = qtcp.C_102(fcs.alpha(DtFrame['cosa']), fcs.beta(DtFrame['tanb']), DtFrame['mH'], DtFrame['mA'], DtFrame['l5'])
    c123 = qtcp.C_123(fcs.alpha(DtFrame['cosa']), fcs.beta(DtFrame['tanb']), DtFrame['mH'], DtFrame['l5'])
    c140 = qtcp.C_140(fcs.alpha(DtFrame['cosa']), fcs.beta(DtFrame['tanb']), DtFrame['mH'], DtFrame['l5'])
    
    return [c93,c94,c102,c123,c140]

def perturbative_unitarity_const(c):
    cmax = 8*np.pi
    res = np.where(np.abs(c)<cmax,False,True)
    
    return res

def calculate_eigenvalues(DtFrame):
    if fcs.alignment:
        THDM2 = anyBSM(THDM_type, scheme_name = 'OSalignment')
    else:
        THDM2 = anyBSM(THDM_type, scheme_name = 'OS')
        
    sino = np.sin(fcs.beta(DtFrame['tanb'])-fcs.alpha(DtFrame['cosa']))
    
    a0=[]
        
    THDM2.progress=False
    THDM2.warnSSSS=False
    
    for i in DtFrame.index:
        if not fcs.alignment and sino[i]==1.0:
            THDM2.load_renormalization_scheme('OSalignment')
            
        a0.append(THDM2.eigSSSS(parameters={'Mh2': DtFrame.at[i,'mH'], 'MAh2': DtFrame.at[i,'mA'], 'MHm2': DtFrame.at[i,'mHpm'], 'TanBeta': DtFrame.at[i,'tanb'], 'SinBmA': sino[i],'M': DtFrame.at[i,'M'], 'MWm': fcs.MW, 'MWp': fcs.MW}))
    
    return a0

def perturbative_unitarity_const_a0(a0):
    a0max = 0.5
    res = np.where(np.abs(a0)<a0max,False,True)
    
    return res

#%%                         Function to call

def main_module(N_points):
    s1 = time()

    Table = fcs.find_random_points(int(N_points))

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
    e1 = time()
    print('Duration - Setting up the points: %f' %(e1-s1))

    #%%                                 Analysis

    ###     Perturbativity tests

    s2 = time()
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
    e2 = time()

    print('Duration - Theoretical bounds: %f' %(e2-s2))

    #%%                         Calculate lambda

    s3 = time()
    # lamb = []
    # lambtree = []

    # sino = np.sin(fcs.beta(TableTot['tanb'])-fcs.alpha(TableTot['cosa']))
    # if fcs.alignment:
    #     THDM2 = anyBSM(THDM_type, scheme_name = 'OSalignment')
    # else:
    #     THDM2 = anyBSM(THDM_type, scheme_name = 'OS')

    # for i in TableTot.index:
    #     if not fcs.alignment and sino[i]==1.0:
    #         THDM2.load_renormalization_scheme('OSalignment')
            
    #     THDM2.setparameters({'Mh2': TableTot.at[i,'mH'], 'MAh2': TableTot.at[i,'mA'], 'MHm2': TableTot.at[i,'mHpm'], 'TanBeta': TableTot.at[i,'tanb'], 'SinBmA': sino[i],'M': TableTot.at[i,'M']}) #Define new mass in anyBSM
    #     dic = THDM2.lambdahhh()
    #     lamb.append(-np.real(dic['total'])/fcs.Gammahhh_treelevel(0, 0))  #Recalculate lambda
    #     lambtree.append(-np.real(dic['treelevel'])/fcs.Gammahhh_treelevel(0, 0))  #Recalculate lambda
    
    lamb, lambtree = calculate_lambda(TableTot)
    sino = np.sin(fcs.beta(TableTot['tanb'])-fcs.alpha(TableTot['cosa']))
    
    # Using that sinBmA ~ 1-x²/2
    kappa_kan_x = fcs.Gammahhh_oneloop(np.sqrt(2*(1-sino)), TableTot['M'], TableTot['mH'], TableTot['mA'], TableTot['mHpm'])/fcs.Gammahhh_treelevel(0, 0)
    kappa_kan = fcs.Gammahhh_oneloop_cos(fcs.beta(TableTot['tanb'])-fcs.alpha(TableTot['cosa']), fcs.beta(TableTot['tanb']), TableTot['M'], TableTot['mH'], TableTot['mA'], TableTot['mHpm'])/fcs.Gammahhh_treelevel(0, 0)
        
    TableTot = pd.concat([TableTot,pd.DataFrame({'kappa': lamb, 'kappa-tree': lambtree, 'kappa-kan-x': kappa_kan_x, 'kappa-kan': kappa_kan})],axis=1)

    e3 = time()
    print('Duration - Calculating kappa: %f' %(e3-s3))

    #%%                                 Calculate S,T,U & collider constraints

    s4 = time()
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
    e4 = time()
    print('Duration - STU & Collider bounds: %f' %(e4-s4))

    #%%                                 Impose bounds from HiggsSignals

    cnd = signals_const(np.array(TableTot_STU_Collid['HiggsS'],dtype=float))
    TableTot_STU_Collid = TableTot_STU_Collid.drop(TableTot_STU_Collid[cnd].index)

    #%%                                 Impose bounds from BSG

    s5 = time()
    cnd = bsg.Constraints_BSG(np.array(TableTot_STU_Collid['tanb'],dtype=float), np.array(TableTot_STU_Collid['mHpm'],dtype=float))
    TableTot_STU_Collid_BSG = TableTot_STU_Collid.drop(TableTot_STU_Collid[cnd].index)
    e5 = time()
    print('Duration - BSG bounds: %f' %(e5-s5))

    return TableTot, TableTot_STU, TableTot_STU_Collid, TableTot_STU_Collid_BSG