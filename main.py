#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Oct 23 17:53:30 2023

@author: leonardoferreira
"""

from anyBSM import anyBSM
import numpy as np
import matplotlib.pyplot as plt

SM = anyBSM('SM',scheme_name = 'OS')
THDM1 = anyBSM('THDMI', scheme_name = 'OS')
THDM2 = anyBSM('THDMII', scheme_name = 'OS')

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
### ----------

def Gammahhh_treelevel(x,M):
    G = -(3*mh**2/v)*(1+(3/2)*(1-4*M**2/(3*mh**2))*x**2)
    
    return G

def Gammahhh_oneloop_SM_like(x,M):
    G = -(3*mh**2/v)*(1+(3/2)*(1-4*M**2/(3*mh**2))*x**2-3*mt**4/(3*np.pi**2*mh**2*v**2))
    
    return G

def Gammahhh_oneloop(x,M,mH,mA,mHpm):
    G = -(3*mh**2/v)*(1+(3/2)*(1-4*M**2/(3*mh**2))*x**2-3*mt**4/(3*np.pi**2*mh**2*v**2)+(mH**4/(12*np.pi**2*mh**2*v**2))*(1-M**2/(mH**2))**3+(mA**4/(12*np.pi**2*mh**2*v**2))*(1-M**2/(mA**2))**3+(mHpm**4/(6*np.pi**2*mh**2*v**2))*(1-M**2/(mHpm**2))**3)
    
    return G

x = np.arcsin(0.999)-np.pi/2

#lambdahhh is defined as -Gammahhh

print(-Gammahhh_treelevel(0, 0))
print(-Gammahhh_oneloop_SM_like(0, 0))

print((-lamb["total"]-Gammahhh_oneloop_SM_like(0, 0))/Gammahhh_treelevel(0, 0))

#%% Plots para checar renormalização

###                         Plot 1

yarray = np.array([0.037,0.0333,0.0314,0.0306,0.0304,0.0308,0.0316,0.0319,0.0328,0.0344,0.0365,0.0397])
xarray = np.array([60,70,80,90,100,110,120,125.1,130,140,150,160])

plt.rc('text', usetex=True)
# plt.rc('xtick',labelsize=20)
# plt.rc('ytick',labelsize=20)

fig, ax = plt.subplots(figsize=(10, 10))

plt.plot(xarray, yarray)
plt.yticks(size=20)
plt.ylabel(r'$\delta \Gamma_{(1)} / \Gamma_0$', size=25)
#plt.ylim((-1,1))
plt.xticks(size=20)
plt.xlabel(r'SM-Higgs mass $m_h$ [GeV]', size=25)
ax.grid()

plt.show()

###                         Plot 1.5

mh_array = np.linspace(-50,50,101)+mh
δΓ = np.copy(mh_array)

for i,x in enumerate(mh_array):
    SM.setparameters({'Mh': x})
    lamb = SM.lambdahhh()
    mh = x
    δΓ[i] = (np.abs(lamb["total"])-(-Gammahhh_oneloop_SM_like(0, 0)))/(-Gammahhh_treelevel(0, 0))


fig, ax = plt.subplots(figsize=(10, 10))

plt.plot(mh_array, δΓ)
plt.yticks(size=20)
plt.ylabel(r'$\delta \Gamma_{(1)} / \Gamma_0$', size=25)
plt.xticks(size=20)
plt.xlabel(r'SM-Higgs mass $m_h$ [GeV]', size=25)
ax.grid()

plt.show()

mh = mhSM
SM = anyBSM('SM',scheme_name = 'OS')

###                         Plot 2

mt_array = np.linspace(-50,50,101)+mt
δΓ = np.copy(mt_array)

for i,x in enumerate(mt_array):
    SM.setparameters({'Mu3': x})
    lamb = SM.lambdahhh()
    mt = x
    δΓ[i] = (np.abs(lamb["total"])-(-Gammahhh_oneloop_SM_like(0, 0)))/(-Gammahhh_treelevel(0, 0))


fig, ax = plt.subplots(figsize=(10, 10))

plt.plot(mt_array, δΓ)
plt.yticks(size=20)
plt.ylabel(r'$\delta \Gamma_{(1)} / \Gamma_0$', size=25)
plt.xticks(size=20)
plt.xlabel(r'Top-quark mass $m_t$ [GeV]', size=25)
ax.grid()

plt.show()

mt = mtSM
SM = anyBSM('SM',scheme_name = 'OS')

#%%                         Plot THDM

mA = 300
mHpm = 550
mH = 300

###                         Plot 3

sBmA_std = 0.8

M_array = np.linspace(0,500,100)
Γ = np.copy(M_array)
THDM2 = anyBSM('THDMII', scheme_name = 'OSalignment')
THDM2.setparameters({'SinBmA': sBmA_std})

for i,x in enumerate(M_array):
    THDM2.setparameters({'M': x})
    lamb = THDM2.lambdahhh()
    Γ[i] = np.abs(lamb["total"])

fig, ax = plt.subplots(figsize=(10, 10))

plt.plot(M_array, Γ,c='k',label='AnyBSM')
plt.plot(M_array, -Gammahhh_treelevel(np.arcsin(sBmA_std)-np.pi/2, M_array),c='y',label='Kanemura (tree-level)')
plt.plot(M_array, -Gammahhh_oneloop(np.arcsin(sBmA_std)-np.pi/2, M_array,mH,mA,mHpm),c='orange',label='Kanemura (one-loop)')
plt.yticks(size=20)
plt.ylabel(r'$\Gamma$ [GeV]', size=25)
plt.xticks(size=20)
plt.xlabel(r'$M$ [GeV]', size=25)
ax.grid()

plt.legend(fontsize=20)

plt.show()

###                         Plot 4

M_std = 200

sBmA_array = np.linspace(0,1,20)
Γ = np.copy(sBmA_array)
THDM2 = anyBSM('THDMII', scheme_name = 'OSalignment')
THDM2.setparameters({'M': M_std})

for i,x in enumerate(sBmA_array):
    THDM2.setparameters({'SinBmA': x})
    lamb = THDM2.lambdahhh()
    Γ[i] = np.abs(lamb["total"])

fig, ax = plt.subplots(figsize=(10, 10))

plt.plot(sBmA_array, Γ,c='k',label='AnyBSM')
plt.plot(sBmA_array, -Gammahhh_treelevel(np.arcsin(sBmA_array)-np.pi/2, M_std),c='y',label='Kanemura (tree-level)')
plt.plot(sBmA_array, -Gammahhh_oneloop(np.arcsin(sBmA_array)-np.pi/2, M_std,mH,mA,mHpm),c='orange',label='Kanemura (one-loop)')
plt.yticks(size=20)
plt.ylabel(r'$\Gamma$ [GeV]', size=25)
plt.xticks(size=20)
plt.xlabel(r'$s_{\beta-\alpha}$', size=25)
ax.grid()

plt.legend(fontsize=20)

plt.show()

#%%                         Plot SM

###                         Plot 5

mh_array = np.linspace(-50,50,101)+mh
Γ = np.copy(mh_array)
Γ0 = np.copy(mh_array)
ΓK = np.copy(mh_array)
Γ0K = np.copy(mh_array)

for i,x in enumerate(mh_array):
    SM.setparameters({'Mh': x})
    lamb = SM.lambdahhh()
    mh = x
    Γ[i] = np.abs(lamb["total"])
    Γ0[i] = np.abs(lamb["treelevel"])
    ΓK[i] = -Gammahhh_oneloop_SM_like(0, 0)
    Γ0K[i] = -Gammahhh_treelevel(0, 0)

#%%

fig, ax = plt.subplots(figsize=(10, 10))

plt.plot(mh_array, Γ0,label='AnyBSM (tree-level)',linestyle='dashed')
plt.plot(mh_array, Γ0K,label='Kanemura (tree-level)',linestyle='dashdot')
plt.plot(mh_array, Γ,label='AnyBSM (one-loop)')
plt.plot(mh_array, ΓK,label='Kanemura (one-loop)')
plt.yticks(size=20)
plt.ylabel(r'$\Gamma$', size=25)
plt.xticks(size=20)
plt.xlabel(r'SM-Higgs mass $m_h$ [GeV]', size=25)
ax.grid()
plt.title('SM - OS Scheme', size=25)

plt.legend(fontsize=20)

plt.show()

mh = mhSM
SM = anyBSM('SM',scheme_name = 'OS')

#%%                         Plot SM

###                         Plot 6

mt_array = np.linspace(-50,50,101)+mt
Γ = np.copy(mt_array)
Γ0 = np.copy(mt_array)
ΓK = np.copy(mt_array)
Γ0K = np.copy(mt_array)

for i,x in enumerate(mt_array):
    SM.setparameters({'Mu3': x})
    lamb = SM.lambdahhh()
    mt = x
    Γ[i] = np.abs(lamb["total"])
    Γ0[i] = np.abs(lamb["treelevel"])
    ΓK[i] = -Gammahhh_oneloop_SM_like(0, 0)
    Γ0K[i] = -Gammahhh_treelevel(0, 0)
    
#%%

fig, ax = plt.subplots(figsize=(10, 10))

plt.plot(mt_array, Γ0,label='AnyBSM (tree-level)',linestyle='dashed')
plt.plot(mt_array, Γ0K,label='Kanemura (tree-level)',linestyle='dashdot')
plt.plot(mt_array, Γ,label='AnyBSM (one-loop)')
plt.plot(mt_array, ΓK,label='Kanemura (one-loop)')
plt.yticks(size=20)
plt.ylabel(r'$\Gamma$', size=25)
plt.xticks(size=20)
plt.xlabel(r'Top-quark mass $m_t$ [GeV]', size=25)
ax.grid()
plt.title('SM - OS Scheme', size=25)

plt.legend(fontsize=20)

plt.show()

mt = mtSM
SM = anyBSM('SM',scheme_name = 'OS')

#%%                         Plot SM

SM = anyBSM('SM',scheme_name = 'MS')

###                         Plot 5

mh_array = np.linspace(-50,50,101)+mh
Γ = np.copy(mh_array)
Γ0 = np.copy(mh_array)
ΓK = np.copy(mh_array)
Γ0K = np.copy(mh_array)

for i,x in enumerate(mh_array):
    SM.setparameters({'Mh': x})
    lamb = SM.lambdahhh()
    mh = x
    Γ[i] = np.abs(lamb["total"])
    Γ0[i] = np.abs(lamb["treelevel"])
    ΓK[i] = -Gammahhh_oneloop_SM_like(0, 0)
    Γ0K[i] = -Gammahhh_treelevel(0, 0)

#%%

fig, ax = plt.subplots(figsize=(10, 10))

plt.plot(mh_array, Γ0,label='AnyBSM (tree-level)',linestyle='dashed')
plt.plot(mh_array, Γ0K,label='Kanemura (tree-level)',linestyle='dashdot')
plt.plot(mh_array, Γ,label='AnyBSM (one-loop)')
plt.plot(mh_array, ΓK,label='Kanemura (one-loop)')
plt.yticks(size=20)
plt.ylabel(r'$\Gamma$', size=25)
plt.xticks(size=20)
plt.xlabel(r'SM-Higgs mass $m_h$ [GeV]', size=25)
ax.grid()
plt.title('SM - MS Scheme', size=25)

plt.legend(fontsize=20)

plt.show()

mh = mhSM
SM = anyBSM('SM',scheme_name = 'OS')

#%%                         Plot SM

SM = anyBSM('SM',scheme_name = 'MS')

###                         Plot 6

mt_array = np.linspace(-50,50,101)+mt
Γ = np.copy(mt_array)
Γ0 = np.copy(mt_array)
ΓK = np.copy(mt_array)
Γ0K = np.copy(mt_array)

for i,x in enumerate(mt_array):
    SM.setparameters({'Mu3': x})
    lamb = SM.lambdahhh()
    mt = x
    Γ[i] = np.abs(lamb["total"])
    Γ0[i] = np.abs(lamb["treelevel"])
    ΓK[i] = -Gammahhh_oneloop_SM_like(0, 0)
    Γ0K[i] = -Gammahhh_treelevel(0, 0)
    
#%%

fig, ax = plt.subplots(figsize=(10, 10))

plt.plot(mt_array, Γ0,label='AnyBSM (tree-level)',linestyle='dashed')
plt.plot(mt_array, Γ0K,label='Kanemura (tree-level)',linestyle='dashdot')
plt.plot(mt_array, Γ,label='AnyBSM (one-loop)')
plt.plot(mt_array, ΓK,label='Kanemura (one-loop)')
plt.yticks(size=20)
plt.ylabel(r'$\Gamma$', size=25)
plt.xticks(size=20)
plt.xlabel(r'Top-quark mass $m_t$ [GeV]', size=25)
ax.grid()
plt.title('SM - MS Scheme', size=25)

plt.legend(fontsize=20)

plt.show()

mt = mtSM
SM = anyBSM('SM',scheme_name = 'OS')