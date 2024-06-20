#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon May 20 15:18:31 2024

@author: leo
"""

import Higgs.predictions as HP
import Higgs.bounds as HB
import Higgs.signals as HS
import Higgs.tools.Input as hinput
import scan_SPheno_funcs as SPfcs
from contextlib import chdir
import datetime
import pandas as pd
import scan_parameterspace_funcs as fcs

keep_log = False
Higgs_tools_database_path = "/home/leo/Documents/Unicamp/HEPTools/higgstools/Database/"

bounds = HB.Bounds(Higgs_tools_database_path+"hbdataset/") # load HB dataset
signals = HS.Signals(Higgs_tools_database_path+"hsdataset/") # load HS dataset

neutralIds = [25,35,36]
chargedIds = [37]
neutralIdStrings = [str(id) for id in neutralIds]
chargedIdStrings = [str(id) for id in chargedIds]

def Higgs_tools_scan():
    dc = hinput.readHB5SLHA(SPfcs.Spheno_path+"SPheno.spc.THDM"+fcs.THDM_type, neutralIds, chargedIds)
    pred = hinput.predictionsFromDict(dc, neutralIdStrings, chargedIdStrings, [])
    h = pred.particle('25')

    # evaluate HiggsBounds
    hbresult = bounds(pred)
    
    # evaluate HiggsSignals
    chisq = signals(pred)
    
    with chdir(Higgs_tools_database_path+'Logs'):
        
        if keep_log and not bool(hbresult):
            fil = open("Higgs_tools_log"+str(datetime.datetime.now()),'w')
            
            fil.write(hbresult)
            fil.write(print(f"HiggsSignals chisq: {chisq}"))
            
            fil.close()
    
    outpt=pd.DataFrame({'HiggsB': bool(hbresult), 'HiggsS': chisq},index=[1])

    return outpt