#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Apr 12 14:25:54 2024

@author: leo
"""

from anyBSM import anyBSM
import numpy as np
import matplotlib.pyplot as plt

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
### ----------

def Gammahhh_treelevel(x,M):
    G = -(3*mh**2/v)*(1+(3/2)*(1-4*M**2/(3*mh**2))*x**2)
    
    return G

#%%

THDM2 = anyBSM('THDMII', scheme_name = 'OSalignment')
THDM2.setparameters({'SinBmA': 1}) #Define new mass in anyBSM

mA_array = np.linspace(200,800,35) #Plot mA GeV
mH_array = np.linspace(200,800,36)
m1, m2 = np.meshgrid(mA_array,mH_array,indexing='ij')

δΓ = np.copy(m1)

for i,x in enumerate(mA_array):
    for j, y in enumerate(mH_array):
        THDM2.setparameters({'Mh2': x, 'MAh2': y}) #Define new mass in anyBSM
        lamb = THDM2.lambdahhh()  #Recalculate lambda
        δΓ[i,j] = (np.real(lamb["total"]))/(-Gammahhh_treelevel(0, 0))


fig, ax = plt.subplots(figsize=(10, 10))

ctrf = plt.contourf(mH_array,mA_array, δΓ)
plt.yticks(size=20)
plt.xlabel(r'$m_H$ [GeV]', size=25)
plt.xticks(size=20)
plt.ylabel(r'$m_A$ [GeV]', size=25)
ax.grid()

fig.colorbar(ctrf, ax=ax)

plt.show()

#%%

THDM2 = anyBSM('THDMII', scheme_name = 'OS')

mH_array = np.linspace(200,800,35) #Plot mA GeV
mHpm_array = np.linspace(200,800,36)
m1, m2 = np.meshgrid(mH_array,mHpm_array,indexing='ij')

δΓ = np.copy(m1)

for i,x in enumerate(mH_array):
    for j, y in enumerate(mHpm_array):
        THDM2.setparameters({'Mh2': x, 'MHm2': y}) #Define new mass in anyBSM
        lamb = THDM2.lambdahhh()  #Recalculate lambda
        δΓ[i,j] = (np.real(lamb["total"]))/(-Gammahhh_treelevel(0, 0))


fig, ax = plt.subplots(figsize=(10, 10))

ctrf = plt.contourf(mH_array,mHpm_array, δΓ.T)
plt.yticks(size=20)
plt.xticks(size=20)
plt.xlabel(r'$m_H$ [GeV]', size=25)
plt.ylabel(r'$m_{H^\pm}$ [GeV]', size=25)
ax.grid()

fig.colorbar(ctrf, ax=ax)

plt.show()

#%%

THDM2 = anyBSM('THDMII', scheme_name = 'OS')

mH_array = np.linspace(200,800,35) #Plot mA GeV
M_array = np.linspace(200,1600,36)
m1, m2 = np.meshgrid(mH_array,M_array,indexing='ij')

δΓ = np.copy(m1)

for i,x in enumerate(mH_array):
    for j, y in enumerate(M_array):
        THDM2.setparameters({'Mh2': x, 'M': y}) #Define new mass in anyBSM
        lamb = THDM2.lambdahhh()  #Recalculate lambda
        δΓ[i,j] = (np.real(lamb["total"]))/(-Gammahhh_treelevel(0, 0))


fig, ax = plt.subplots(figsize=(10, 10))

ctrf = plt.contourf(mH_array,M_array, δΓ.T)
plt.yticks(size=20)
plt.xticks(size=20)
plt.xlabel(r'$m_H$ [GeV]', size=25)
plt.ylabel(r'$M$ [GeV]', size=25)
ax.grid()

fig.colorbar(ctrf, ax=ax)

plt.show()

#%%

THDM2 = anyBSM('THDMII', scheme_name = 'OS')

sBA_array = np.linspace(-0.999,0.999,35) #Plot mA GeV
M_array = np.linspace(200,800,36)
m1, m2 = np.meshgrid(M_array,sBA_array,indexing='ij')

δΓ = np.copy(m1)

for i,x in enumerate(M_array):
    for j, y in enumerate(sBA_array):
        THDM2.setparameters({'M': x, 'SinBmA': y}) #Define new mass in anyBSM
        lamb = THDM2.lambdahhh()  #Recalculate lambda
        δΓ[i,j] = (np.real(lamb["total"]))/(-Gammahhh_treelevel(0, 0))


fig, ax = plt.subplots(figsize=(10, 10))

ctrf = plt.contourf(M_array,sBA_array, δΓ.T)
plt.yticks(size=20)
plt.xticks(size=20)
plt.xlabel(r'$M$ [GeV]', size=25)
plt.ylabel(r'$\sin_{\beta-\alpha}$', size=25)
ax.grid()

fig.colorbar(ctrf, ax=ax)

plt.show()
