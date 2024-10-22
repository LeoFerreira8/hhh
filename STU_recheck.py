#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Oct 18 18:37:24 2024

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

set_dir = 'data_'+'THDM'+THDM_type+strgl5+'-'+strga+'-D/'

#TableTot = pd.read_csv('./'+set_dir+'/THDM'+THDM_type+strgl5+'-'+strga+'-Theo_PDG.csv')
TableTot_STU = pd.read_csv('./'+set_dir+'/THDM'+THDM_type+strgl5+'-'+strga+'-STU_PDG.csv')
TableTot_STU_Collid = pd.read_csv('./'+set_dir+'/THDM'+THDM_type+strgl5+'-'+strga+'-Collid_PDG.csv')
TableTot_STU_Collid_BSG = pd.read_csv('./'+set_dir+'/THDM'+THDM_type+strgl5+'-'+strga+'-BSG_PDG.csv')
TableTot_STU_Collid_BSG_unit = pd.read_csv('./'+set_dir+'/THDM'+THDM_type+strgl5+'-'+strga+'-PU_PDG.csv')

#TableTot = TableTot.loc[:10e3]

proxies = []
proxies1 = []
proxies2 = []
proxies3 = []

cnd = spr.STU_constraint(TableTot_STU['S-parameter (1-loop BSM)'],TableTot_STU['T-parameter (1-loop BSM)'],TableTot_STU['U-parameter (1-loop BSM)'])
proxies = TableTot_STU.drop(TableTot_STU[cnd].index)

# if not proxies.empty:
#     cnd = spr.collider_const(proxies[i]['HiggsB'])
#     proxies1.append(proxies[i].drop(proxies[i][cnd].index))
# else:
#     proxies1.append(proxies[i])

#     #%%                                 Impose bounds from HiggsSignals
    
# if not proxies1[i].empty:
#     cnd = spr.signals_const(np.array(proxies1[i]['HiggsS'],dtype=float))
#     proxies1[i] = proxies1[i].drop(proxies1[i][cnd].index)