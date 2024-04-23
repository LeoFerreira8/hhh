#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Mar 20 14:46:05 2024

@author: leonardoferreira
"""

from anyBSM import anyBSM
import numpy as np
import matplotlib.pyplot as plt

SM = anyBSM('SM',scheme_name = 'OS')
THDM1 = anyBSM('THDMI', scheme_name = 'OS')
THDM2 = anyBSM('THDMII', scheme_name = 'OSalignment')

lamb = SM.lambdahhh()
lamb1 = THDM1.lambdahhh()
lamb2 = THDM2.lambdahhh()

#%% Cross-check

mhSM = 125.1 #GeV
mh = mhSM

### ---------      SM param
MW = 80.379
MZ = 91.187
mtSM = 172.5
mt = mtSM
alphaQED = 137.035999679

e = 2*np.sqrt(np.pi/(alphaQED))
vSM = 2*MW*np.sqrt(1-MW**2/MZ**2)/e
v = vSM
### ----------    Standard definitions in anyBSM
mA = 300 #Pseudoscalar mass
#mHpm = 550 #Charged scalar mass
mHpm = 300 #Charged scalar mass
mH = 300  #Heavy Higgs mass
sBmA = 1 #Standard definition of sin(β-α) in anyBSM
x = np.arcsin(sBmA)-np.pi/2
M = 300 #Standard definition of M in anyBSM
BmA = np.arcsin(sBmA)
Beta = np.arctan(2)

### ----------

def Gammahhh_treelevel(x,M):
    G = -(3*mh**2/v)*(1+(3/2)*(1-4*M**2/(3*mh**2))*x**2)
    
    return G

def Gammahhh_oneloop_SM_like(x,M):
    G = -(3*mh**2/v)*(1+(3/2)*(1-4*M**2/(3*mh**2))*x**2-3*mt**4/(3*np.pi**2*mh**2*v**2))
    
    return G

def Gammahhh_oneloop(x,M,mH,mA,mHpm):
    G = -(3*mh**2/v)*(1+(3/2)*(1-4*M**2/(3*mh**2))*x**2-3*mt**4/(3*np.pi**2*mh**2*v**2)+(mH**4/(12*np.pi**2*mh**2*v**2))*(1-M**2/(mH**2))**3+(mA**4/(12*np.pi**2*mh**2*v**2))*(1-M**2/(mA**2))**3+(mHpm**4/(6*np.pi**2*mh**2*v**2))*(1-M**2/(mHpm**2))**3+(mh**4/(2*np.pi**2*mh**2*v**2)))
    
    return G

#lambdahhh is defined as -Gammahhh

#%% Plots para checar renormalização

###                         Plot 1 - SM comparison with 1-loop as function of mh

mA_array = np.linspace(-150,150,101)+mA #Plot mhSM +-50 GeV
Γ = np.copy(mA_array)
Γ1 = np.copy(mA_array)

for i,mass in enumerate(mA_array):
    THDM2.setparameters({'MAh2': mass, 'MHm2': mHpm, 'Mh2': mH, 'M': M, 'SinBmA': sBmA}) #Define new mass in anyBSM
    lamb = THDM2.lambdahhh()  #Recalculate lambda
    Γ[i] = np.real(lamb["total"])
    Γ1[i] = -np.real(Gammahhh_oneloop(x, M, mH, mass, mHpm))

#%%

fig, ax = plt.subplots(figsize=(10, 10))

plt.plot(mA_array, Γ,label="AnyBSM")
plt.plot(mA_array, Γ1,label="Kanemura")
plt.yticks(size=20)
plt.ylabel(r'$\Gamma$', size=25)
plt.xticks(size=20)
plt.xlabel(r'Mass $m_\phi$ [GeV]', size=25)
#plt.ylim(150,330)
ax.grid()

plt.legend(fontsize=20)

plt.show()