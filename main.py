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

###                         Plot 1 - SM comparison with 1-loop as function of mh

mh_array = np.linspace(-50,50,101)+mhSM #Plot mhSM +-50 GeV
δΓ = np.copy(mh_array)

for i,x in enumerate(mh_array):
    SM.setparameters({'Mh': x}) #Define new mass in anyBSM
    lamb = SM.lambdahhh()  #Recalculate lambda
    mh = x  #Redefine mh globally to recalculate Gamma
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

###                         Plot 2 - SM comparison with 1-loop as function of mt

mt_array = np.linspace(-50,50,101)+mtSM
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

#%%                         Plots 2HDM

#Standard definitions in anyBSM
mA = 300 #Pseudoscalar mass
mHpm = 550 #Charged scalar mass
mH = 300  #Heavy Higgs mass
sBmA_std = 0.9 #Standard definition of sin(β-α) in anyBSM
M_std = 200 #Standard definition of M in anyBSM

###                         Plot 3 - 2HDM 1-loop comparison as function of M

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

###                         Plot 4 - 2HDM 1-loop comparison as function of sβα

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

#%%                         Plots SM - Tree and 1-loop level comparison in OS scheme

###                         Plot 5 - As function of mh

mh_array = np.linspace(-50,50,101)+mhSM
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

#%%                         Plot SM - Tree and 1-loop level comparison in OS scheme

###                         Plot 6 - As function of mt

mt_array = np.linspace(-50,50,101)+mtSM
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

#%%                         Plot SM - Tree and 1-loop level comparison in MS scheme

SM = anyBSM('SM',scheme_name = 'MS')

###                         Plot 5 - As function of mh

mh_array = np.linspace(-50,50,101)+mhSM
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

#%%                         Plot SM - Tree and 1-loop level comparison in MS scheme

SM = anyBSM('SM',scheme_name = 'MS')

###                         Plot 6 - As function of mt

mt_array = np.linspace(-50,50,101)+mtSM
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

#%%                         Plot SM - Tree and 1-loop level comparison in MS scheme (varying the scale as well)

SM = anyBSM('SM',scheme_name = 'MS')

###                         Plot 7 - As function of mt

mt_array = np.linspace(-50,50,101)+mtSM
Γ = np.copy(mt_array)
Γ0 = np.copy(mt_array)
ΓK = np.copy(mt_array)
Γ0K = np.copy(mt_array)

for i,x in enumerate(mt_array):
    SM.setparameters({'Mu3': x})
    SM.setparameters(params={'Qren': x}) #Varying the RS
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
plt.title('SM - MS Scheme (scale variation)', size=25)

plt.legend(fontsize=20)

plt.show()

mt = mtSM

#%%                         Plot SM - Tree and 1-loop level comparison in MS scheme (varying the scale as well)

SM = anyBSM('SM',scheme_name = 'MS')

###                         Plot 8 - As function of mh

mh_array = np.linspace(-50,150,101)+mhSM
Γ = np.copy(mh_array)
Γ0 = np.copy(mh_array)
ΓK = np.copy(mh_array)
Γ0K = np.copy(mh_array)

for i,x in enumerate(mh_array):
    SM.setparameters({'Mh': x})
    SM.setparameters(params={'Qren': x})
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
plt.title('SM - MS Scheme (scale variation)', size=25)

plt.legend(fontsize=20)

plt.show()

mh = mhSM

#%%                         Plot SM - Tree and 1-loop level comparison in MS scheme (varying the scale as well)

SM = anyBSM('SM',scheme_name = 'MS')

###                         Plot 9 - As function of the scale

mh_array = np.linspace(-50,200,101)+mh
Γ = np.copy(mh_array)
Γ0 = np.copy(mh_array)
ΓK = np.copy(mh_array)
Γ0K = np.copy(mh_array)

for i,x in enumerate(mh_array):
    SM.setparameters(params={'Qren': x})
    lamb = SM.lambdahhh()
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
plt.xlabel(r'$\mu$ [GeV]', size=25)
ax.grid()
plt.title('SM - MS Scheme', size=25)

plt.legend(fontsize=20)

plt.show()

mh = mhSM