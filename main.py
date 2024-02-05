#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Oct 23 17:53:30 2023

@author: leonardoferreira
"""

from anyBSM import anyBSM
import numpy as np

SM = anyBSM('SM',scheme_name = 'MS')

THDM1 = anyBSM('THDMI', scheme_name = 'MS')
THDM2 = anyBSM('THDMII', scheme_name = 'MS')

lamb = SM.lambdahhh()
lamb1 = THDM1.lambdahhh()
lamb2 = THDM2.lambdahhh()

#%% Cross-check

mh = 125.1 #GeV

### ---------      SM param
MW = 80.379
MZ = 91.187
alphaQED = 137.035999679

e = 2*np.sqrt(np.pi/(alphaQED))
v = 2*MW*np.sqrt(1-MW**2/MZ**2)/e
### ----------

def Gammahhh_treelevel(mH,x,M):
    G = -(3*mh**2/v)*(1+(3/2)*(1-4*M**2/(3*mh**2))*x**2)
    
    return G

x = np.arcsin(0.9)-np.pi/2

#lambdahhh is defined as -Gammahhh

print(-Gammahhh_treelevel(300, x, 0))