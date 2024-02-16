#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Oct 23 17:53:30 2023

@author: leonardoferreira
"""

from anyBSM import anyBSM
import numpy as np

SM = anyBSM('SM',scheme_name = 'OS')
SM.set_evaluation_mode('analytical')
hself = SM.Sigma('h', simplify = False,draw=True) # this returns a string
# just convert the string to sympy
SM.sympify(hself, simplify = False)
# convert to sympy and replace all internal UFO parameters with

THDM1 = anyBSM('THDMI', scheme_name = 'OS')
THDM2 = anyBSM('THDMII', scheme_name = 'OS')
THDM2.set_evaluation_mode('analytical')
hself2 = THDM2.Sigma('h', simplify = False,draw=True) # this returns a string
# just convert the string to sympy
THDM2.sympify(hself2, simplify = True)
# convert to sympy and replace all internal UFO parameters with

lamb = SM.lambdahhh()
lamb1 = THDM1.lambdahhh()
lamb2 = THDM2.lambdahhh()

#%% Cross-check

mh = 125.1 #GeV

### ---------      SM param
MW = 80.379
MZ = 91.187
mt = 172.5
alphaQED = 137.035999679

e = 2*np.sqrt(np.pi/(alphaQED))
v = 2*MW*np.sqrt(1-MW**2/MZ**2)/e
### ----------

def Gammahhh_treelevel(x,M):
    G = -(3*mh**2/v)*(1+(3/2)*(1-4*M**2/(3*mh**2))*x**2)
    
    return G

def Gammahhh_oneloop_SM_like(x,M):
    G = -(3*mh**2/v)*(1+(3/2)*(1-4*M**2/(3*mh**2))*x**2-3*mt**4/(3*np.pi**2*mh**2*v**2))
    
    return G

x = np.arcsin(0.98)-np.pi/2

#lambdahhh is defined as -Gammahhh

print(-Gammahhh_treelevel(0, 0))
print(-Gammahhh_oneloop_SM_like(0, 0))