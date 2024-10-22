#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""

MW recalculation

Created on Thu Aug  8 10:32:12 2024

@author: leo
"""

import scan_parameterspace as spr
import scan_parameterspace_funcs as fcs
import bsg
import pandas as pd
import numpy as np
from time import time
from pathlib import Path

CDF = {'Name': 'CDF',       #PDG 2024 CDF 2022
       'MW': 80.4335,
       'ΔS': 0.15,
       'δS': 0.08,
       'ΔT': 0.27,
       'δT': 0.06,
       'corr': 0.93
       }
    
Combined = {'Name': 'Combined',        #J. de Blas, M. Pierini, L. Reina, and L. Silvestrini, Phys. Rev. Lett. 129, 271801 (2022),29 2204.04204.
           'MW': 80.4133,
           'ΔS': 0.086,
           'δS': 0.077,
           'ΔT': 0.177,
           'δT': 0.070,
           'corr': 0.89
           } 

PDG = {'Name': 'PDG',        #PDG 2024
           'MW': 80.3692,
           'ΔS': -0.05,
           'δS': 0.07,
           'ΔT': 0.00,
           'δT': 0.06,
           'corr': 0.93
           }

###                     Model and parameters

THDM_type = fcs.THDM_type
small_l5 = fcs.small_l5
alignment = fcs.alignment
non_alignment_max = fcs.non_alignment_max

if THDM_type=='':
    latex_model = '2HDM - Type I'
    THDM_type=''
elif THDM_type=='II' and not small_l5:
    latex_model = '2HDM - Type II'
elif THDM_type=='II' and small_l5:
    latex_model = '2HDM - 331 EFT'
    
if alignment:
    latex_alignment = r'($\beta-\alpha=\frac{\pi}{2}$)'
else:
    latex_alignment = r'($\beta-\alpha\in[\frac{\pi}{2}\pm%.2f]$)' %(non_alignment_max,)

if alignment:
    strga = 'A'
else:
    strg = str(non_alignment_max)
    strga = 'NA'+strg
    
if small_l5:
    strgl5 = '331lk'
else:
    strgl5 = ''

set_dir = 'data_'+'THDM'+THDM_type+strgl5+'-'+strga+'/'

TableTot = pd.read_csv('./'+set_dir+'/THDM'+THDM_type+strgl5+'-'+strga+'-Theo_PDG.csv')

#TableTot = TableTot.loc[:10e3]

PDG['STUdata'] = TableTot
CDF['STUdata'] = TableTot
Combined['STUdata'] = TableTot

proxies = []
proxies1 = []
proxies2 = []
proxies3 = []

dicts = [CDF,Combined]

for i,d in enumerate(dicts):
    fcs.MW = d['MW']
    spr.s_mean=d['ΔS']
    spr.δs = d['δS']
    spr.t_mean=d['ΔT']
    spr.δt = d['δT']
    spr.corr = d['corr']

    cnd = spr.STU_constraint(d['STUdata']['S-parameter (1-loop BSM)'],d['STUdata']['T-parameter (1-loop BSM)'],d['STUdata']['U-parameter (1-loop BSM)'])
    proxies.append(d['STUdata'].drop(d['STUdata'][cnd].index))
    
    if not proxies[i].empty:
        cnd = spr.collider_const(proxies[i]['HiggsB'])
        proxies1.append(proxies[i].drop(proxies[i][cnd].index))
    else:
        proxies1.append(proxies[i])

    #%%                                 Impose bounds from HiggsSignals
    
    if not proxies1[i].empty:
        cnd = spr.signals_const(np.array(proxies1[i]['HiggsS'],dtype=float))
        proxies1[i] = proxies1[i].drop(proxies1[i][cnd].index)
    
    #%%                                 Impose bounds from BSG
    
    if not proxies1[i].empty:
        cnd = bsg.Constraints_BSG(np.array(proxies1[i]['tanb'],dtype=float), np.array(proxies1[i]['mHpm'],dtype=float))
        proxies2.append(proxies1[i].drop(proxies1[i][cnd].index))
    else:
        proxies2.append(proxies1[i])
    
    #%%                                 Calculate lambda
    
    if not proxies2[i].empty:
        lamb, lambtree = spr.calculate_lambda(proxies2[i])
        sino = np.sin(fcs.beta(proxies2[i]['tanb'])-fcs.alpha(proxies2[i]['cosa']))
    
        # Using that sinBmA ~ 1-x²/2
        kappa_kan_x = fcs.Gammahhh_oneloop(np.sqrt(2*(1-sino)), proxies2[i]['M'], proxies2[i]['mH'], proxies2[i]['mA'], proxies2[i]['mHpm'])/fcs.Gammahhh_treelevel(0, 0)
        kappa_kan = fcs.Gammahhh_oneloop_cos(fcs.beta(proxies2[i]['tanb'])-fcs.alpha(proxies2[i]['cosa']), fcs.beta(proxies2[i]['tanb']), proxies2[i]['M'], proxies2[i]['mH'], proxies2[i]['mA'], proxies2[i]['mHpm'])/fcs.Gammahhh_treelevel(0, 0)
        
        proxies2[i] = pd.concat([proxies2[i],pd.DataFrame({'kappa(new)': lamb, 'kappa-tree(new)': lambtree, 'kappa-kan-x(new)': kappa_kan_x, 'kappa-kan(new)': kappa_kan})],axis=1)
    
    #%%                                 PU bounds
    
    if not proxies2[i].empty:
        proxies2[i]=proxies2[i].reset_index(drop=True)
        TableTot_STU_Collid_BSG_unit = pd.concat([proxies2[i],pd.DataFrame(np.abs(np.array(spr.calculate_quartics(proxies2[i]))).T,columns=['c93','c94','c102','c123','c140'])],axis=1)
    
        TableTot_STU_Collid_BSG_unit = pd.concat([TableTot_STU_Collid_BSG_unit,pd.DataFrame(np.array(spr.calculate_eigenvalues(TableTot_STU_Collid_BSG_unit)).T,columns=['a0'])],axis=1)
    
        cnd = spr.perturbative_unitarity_const_a0(TableTot_STU_Collid_BSG_unit['a0'])
        proxies3.append(TableTot_STU_Collid_BSG_unit.drop(TableTot_STU_Collid_BSG_unit[cnd].index))
        

#%%

CDF['STUdata'],Combined['STUdata'] = proxies#,PDG['STUdata'] = proxies    
CDF['Colliderdata'],Combined['Colliderdata'] = proxies1#,PDG['Colliderdata'] = proxies1
CDF['BSGdata'],Combined['BSGdata'] = proxies2#,PDG['BSGdata'] = proxies2

#####       Save CDF

set_dir_CDF = Path('./data_'+'THDM'+THDM_type+strgl5+'-'+strga+'_CDF_'+'/')
set_dir_CDF.mkdir(parents=True, exist_ok=True)

path_files = [Path('./'+set_dir_CDF.parts[0]+'/THDM'+fcs.THDM_type+strgl5+'-'+strga+'-Theo.csv'),Path('./'+set_dir_CDF.parts[0]+'/THDM'+fcs.THDM_type+strgl5+'-'+strga+'-STU.csv'),Path('./'+set_dir_CDF.parts[0]+'/THDM'+fcs.THDM_type+strgl5+'-'+strga+'-Collid.csv'),Path('./'+set_dir_CDF.parts[0]+'/THDM'+fcs.THDM_type+strgl5+'-'+strga+'-BSG.csv')]

# TableTot.to_csv(path_files[0],index=False)
# CDF['STUdata'].to_csv(path_files[1],index=False)
# CDF['Colliderdata'].to_csv(path_files[2],index=False)
# CDF['BSGdata'].to_csv(path_files[3],index=False)

#####       Save Combined

set_dir_Comb = Path('./data_'+'THDM'+THDM_type+strgl5+'-'+strga+'_Comb_'+'/')
set_dir_Comb.mkdir(parents=True, exist_ok=True)

path_files = [Path('./'+set_dir_Comb.parts[0]+'/THDM'+fcs.THDM_type+strgl5+'-'+strga+'-Theo.csv'),Path('./'+set_dir_Comb.parts[0]+'/THDM'+fcs.THDM_type+strgl5+'-'+strga+'-STU.csv'),Path('./'+set_dir_Comb.parts[0]+'/THDM'+fcs.THDM_type+strgl5+'-'+strga+'-Collid.csv'),Path('./'+set_dir_Comb.parts[0]+'/THDM'+fcs.THDM_type+strgl5+'-'+strga+'-BSG.csv')]

# TableTot.to_csv(path_files[0],index=False)
# Combined['STUdata'].to_csv(path_files[1],index=False)
# Combined['Colliderdata'].to_csv(path_files[2],index=False)
# Combined['BSGdata'].to_csv(path_files[3],index=False)

#####       Save PDG

# set_dir_PDG = Path('./data_'+'THDM'+THDM_type+strgl5+'-'+strga+'/')
# set_dir_PDG.mkdir(parents=True, exist_ok=True)

# path_files = [Path('./'+set_dir_PDG.parts[0]+'/THDM'+fcs.THDM_type+strgl5+'-'+strga+'-Theo_PDG.csv'),Path('./'+set_dir_PDG.parts[0]+'/THDM'+fcs.THDM_type+strgl5+'-'+strga+'-STU_PDG.csv'),Path('./'+set_dir_PDG.parts[0]+'/THDM'+fcs.THDM_type+strgl5+'-'+strga+'-Collid_PDG.csv'),Path('./'+set_dir_PDG.parts[0]+'/THDM'+fcs.THDM_type+strgl5+'-'+strga+'-BSG_PDG.csv')]

# TableTot.to_csv(path_files[0],index=False)
# PDG['STUdata'].to_csv(path_files[1],index=False)
# PDG['Colliderdata'].to_csv(path_files[2],index=False)
# PDG['BSGdata'].to_csv(path_files[3],index=False)