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

THDM2.setparameters({'MAh2': 200, 'MHm2': 800, 'Mh2': 800, 'M': 200, 'SinBmA': 0.8, 'TanBeta': 40, 'Mh1': 125.6}) #Define new mass in anyBSM

#lamb = SM.lambdahhh()
#lamb1 = THDM1.lambdahhh()
#THDM2.set_evaluation_mode('analytical')
#lamb2 = THDM2.lambdahhh(draw=0.50,simplify=True)
#expre = THDM2.SolveDependencies(lamb2["total"])
#lamb2 = THDM2.lambdahhh()

plt.rc('text', usetex=True)
plt.rc('xtick',labelsize=40)
plt.rc('ytick',labelsize=40)

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
lambdahhhSM = 3*mh**2/v
### ----------

def Gammahhh_treelevel(x,M):
    G = -(3*mh**2/v)*(1+(3/2)*(1-4*M**2/(3*mh**2))*x**2)
    
    return G

def Gammahhh_oneloop_SM_like_full(x,M):
    G = -(3*mh**2/v)*(1+(3/2)*(1-4*M**2/(3*mh**2))*x**2-3*mt**4/(3*np.pi**2*mh**2*v**2)+mh**4/(2*np.pi**2*mh**2*v**2))
    
    return G

def Gammahhh_oneloop_SM_like(x,M):
    G = -(3*mh**2/v)*(1+(3/2)*(1-4*M**2/(3*mh**2))*x**2-3*mt**4/(3*np.pi**2*mh**2*v**2))
    
    return G

def Gammahhh_oneloop(x,M,mH,mA,mHpm):
    G = -(3*mh**2/v)*(1+(3/2)*(1-4*M**2/(3*mh**2))*x**2-3*mt**4/(3*np.pi**2*mh**2*v**2)+(mH**4/(12*np.pi**2*mh**2*v**2))*(1-M**2/(mH**2))**3+(mA**4/(12*np.pi**2*mh**2*v**2))*(1-M**2/(mA**2))**3+(mHpm**4/(6*np.pi**2*mh**2*v**2))*(1-M**2/(mHpm**2))**3)
    
    return G

def Gammahhh_treelevel_cos(ba,b,M):
    G = -(3*mh**2/(2*v*np.sin(2*b)))*(np.cos(-3*ba+2*b)+3*np.cos(-ba+2*b)-4*np.cos(ba)**2*np.cos(-ba+2*b)*(M**2/mh**2))
    
    return G

def Gammahhh_oneloop_cos_kan(ba,b,M,mH,mA,mHpm):
    G = Gammahhh_treelevel_cos(ba, b, M)-(3*mh**2/v)*(-3*mt**4/(3*np.pi**2*mh**2*v**2)+(mH**4/(12*np.pi**2*mh**2*v**2))*(1-M**2/(mH**2))**3+(mA**4/(12*np.pi**2*mh**2*v**2))*(1-M**2/(mA**2))**3+(mHpm**4/(6*np.pi**2*mh**2*v**2))*(1-M**2/(mHpm**2))**3)
    #G = Gammahhh_treelevel_cos(ba, b, M)-(3*mh**2/v)*(-3*mt**4/(3*np.pi**2*mh**2*v**2)+(mH**4/(12*np.pi**2*mh**2*v**2))*(1-M**2/(mH**2))**3+(mA**4/(12*np.pi**2*mh**2*v**2))*(1-M**2/(mA**2))**3+(mHpm**4/(6*np.pi**2*mh**2*v**2))*(1-M**2/(mHpm**2))**3)
    
    return G

def Gammahhh_oneloop_cos(ba,b,M,mH,mA,mHpm):
    G = Gammahhh_treelevel_cos(ba, b, M)-(3*mh**2/v)*(-(np.cos(ba)/np.tan(b)+np.sin(ba))**3*3*mt**4/(3*np.pi**2*mh**2*v**2)+(mH**4/(12*np.pi**2*mh**2*v**2))*(1-M**2/(mH**2))**3+(mA**4/(12*np.pi**2*mh**2*v**2))*(1-M**2/(mA**2))**3+(mHpm**4/(6*np.pi**2*mh**2*v**2))*(1-M**2/(mHpm**2))**3)
    #G = Gammahhh_treelevel_cos(ba, b, M)-(3*mh**2/v)*(-3*mt**4/(3*np.pi**2*mh**2*v**2)+(mH**4/(12*np.pi**2*mh**2*v**2))*(1-M**2/(mH**2))**3+(mA**4/(12*np.pi**2*mh**2*v**2))*(1-M**2/(mA**2))**3+(mHpm**4/(6*np.pi**2*mh**2*v**2))*(1-M**2/(mHpm**2))**3)
    
    return G

def Gammahhh_oneloop_cos_correc(ba,b,M,mH,mA,mHpm):
    G = Gammahhh_treelevel_cos(ba, b, M)-(3*mh**2/v)*(-(np.cos(ba)/np.tan(b)+np.sin(ba))**3*3*mt**4/(3*np.pi**2*mh**2*v**2)+mh**4/(2*np.pi**2*mh**2*v**2)+(mH**4/(12*np.pi**2*mh**2*v**2))*(1-M**2/(mH**2))**3+(mA**4/(12*np.pi**2*mh**2*v**2))*(1-M**2/(mA**2))**3+(mHpm**4/(6*np.pi**2*mh**2*v**2))*(1-M**2/(mHpm**2))**3)
    #G = Gammahhh_treelevel_cos(ba, b, M)-(3*mh**2/v)*(-3*mt**4/(3*np.pi**2*mh**2*v**2)+(mH**4/(12*np.pi**2*mh**2*v**2))*(1-M**2/(mH**2))**3+(mA**4/(12*np.pi**2*mh**2*v**2))*(1-M**2/(mA**2))**3+(mHpm**4/(6*np.pi**2*mh**2*v**2))*(1-M**2/(mHpm**2))**3)
    
    return G

#%% Plots para checar renormalização

###                         Plot 1 - SM comparison with 1-loop as function of mh

mh_array = np.linspace(-50,100,101)+mhSM #Plot mhSM +-50 GeV
δΓ = np.copy(mh_array)
δΓ2 = np.copy(mh_array)

for i,x in enumerate(mh_array):
    SM.setparameters({'Mh': x}) #Define new mass in anyBSM
    lamb = SM.lambdahhh()  #Recalculate lambda
    mh = x  #Redefine mh globally to recalculate Gamma
    kappa_kan = (-Gammahhh_oneloop_SM_like(0, 0))/(-Gammahhh_treelevel(0, 0))
    kappa_kan_full = (-Gammahhh_oneloop_SM_like_full(0, 0))/(-Gammahhh_treelevel(0, 0))
    kappa_any = np.real(lamb["total"])/(-Gammahhh_treelevel(0, 0))
    δΓ[i] = (kappa_any-kappa_kan)/kappa_kan
    δΓ2[i] = (kappa_any-kappa_kan_full)/kappa_kan_full


fig, ax = plt.subplots(figsize=(10, 10))

plt.plot(mh_array, δΓ,label='Without Higgs loop correction')
plt.plot(mh_array, δΓ2,label='With Higgs loop correction')
plt.yticks(size=20)
plt.ylabel(r'$(\kappa_\lambda^{{Kan}}-\kappa_\lambda)/\kappa_\lambda^{{Kan}}$', size=25)
plt.xticks(size=20)
plt.xlabel(r'SM-Higgs mass $m_h$ [GeV]', size=25)
ax.grid()

plt.legend(fontsize=20)

plt.show()

mh = mhSM
SM = anyBSM('SM',scheme_name = 'OS')


#%% Plots para checar higgs corrections

###                         Plot 1 - THDM comparison with 1-loop as function of mh

mh_array = np.linspace(-50,50,101)+mhSM #Plot mhSM +-50 GeV
δΓ = np.copy(mh_array)
δΓ2 = np.copy(mh_array)

#Standard definitions in anyBSM
mA = 300 #Pseudoscalar mass
mHpm = 300 #Charged scalar mass
mH = 300  #Heavy Higgs mass
sBmA = 0.9 #Standard definition of sin(β-α) in anyBSM
M = 600 #Standard definition of M in anyBSM
BmA = np.arcsin(sBmA)
Beta = np.arctan(10)

###                         Plot 3 - 2HDM 1-loop comparison as function of M

M_array = np.linspace(0,500,100)
Γ = np.copy(M_array)
Γ0 = np.copy(M_array)
if sBmA == 1:
    scheme_str = 'OSalignment'
else:
    scheme_str = 'OS'

THDM2 = anyBSM('THDMI', scheme_name = scheme_str)
THDM2.setparameters({'MAh2': mA, 'MHm2': mHpm, 'Mh2': mH, 'M': M, 'SinBmA': sBmA, 'TanBeta': np.tan(Beta)}) #Define new mass in anyBSM


for i,x in enumerate(mh_array):
    THDM2.setparameters({'Mh1': x}) #Define new mass in anyBSM
    lamb = THDM2.lambdahhh()  #Recalculate lambda
    mh = x  #Redefine mh globally to recalculate Gamma
    kappa_kan = (-Gammahhh_oneloop_cos(np.arcsin(sBmA), Beta, M, mH, mA, mHpm))/(-Gammahhh_treelevel(0, 0))
    kappa_kan_full = (-Gammahhh_oneloop_cos_correc(np.arcsin(sBmA), Beta, M, mH, mA, mHpm))/(-Gammahhh_treelevel(0, 0))
    kappa_any = np.real(lamb["total"])/(-Gammahhh_treelevel(0, 0))
    δΓ[i] = (kappa_any-kappa_kan)/kappa_kan
    δΓ2[i] = (kappa_any-kappa_kan_full)/kappa_kan_full


fig, ax = plt.subplots(figsize=(10, 10))

plt.plot(mh_array, δΓ,label='Without Higgs loop correction')
plt.plot(mh_array, δΓ2,label='With Higgs loop correction')
plt.yticks(size=20)
plt.ylabel(r'$(\kappa_\lambda^{{Kan}}-\kappa_\lambda)/\kappa_\lambda^{{Kan}}$', size=25)
plt.xticks(size=20)
plt.xlabel(r'SM-Higgs mass $m_h$ [GeV]', size=25)
ax.grid()

textstr = '\n'.join((
    r'\textbf{2HDM-II}',
    r'$M_A=%d$' %(mA),
    r'$M_{H}=%d$' %(mH),
    r'$M_{H^\pm}=%d$' %(mHpm),
    r'$M=%d$' %(M),
    r'$\sin(\beta-\alpha)=%.2f$' %(sBmA),
    r'$\tan(\beta)=%.1f$' %(np.tan(Beta)),
    ))

ax.text(0.70, 0.5, textstr, transform=ax.transAxes, fontsize=18, verticalalignment='top')

plt.legend(fontsize=20)

plt.show()

mh = mhSM
SM = anyBSM('SM',scheme_name = 'OS')

#%%

###                         Plot 2 - SM comparison with 1-loop as function of mt

mt_array = np.linspace(-50,450,101)+mtSM
δΓ = np.copy(mt_array)

for i,x in enumerate(mt_array):
    SM.setparameters({'Mu3': x})
    lamb = SM.lambdahhh()
    mt = x
    δΓ[i] = (np.real(lamb["total"])-(-Gammahhh_oneloop_SM_like(0, 0)))/(-Gammahhh_treelevel(0, 0))


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
mHpm = 300 #Charged scalar mass
mH = 300  #Heavy Higgs mass
sBmA = 1 #Standard definition of sin(β-α) in anyBSM
M = 300 #Standard definition of M in anyBSM
BmA = np.arcsin(sBmA)
Beta = np.arctan(1)

###                         Plot 3 - 2HDM 1-loop comparison as function of M

M_array = np.linspace(0,500,100)
Γ = np.copy(M_array)
Γ0 = np.copy(M_array)
if sBmA == 1:
    scheme_str = 'OSalignment'
else:
    scheme_str = 'OS'
    

THDM2 = anyBSM('THDMI', scheme_name = scheme_str)
THDM2.setparameters({'MAh2': mA, 'MHm2': mHpm, 'Mh2': mH, 'M': M, 'SinBmA': sBmA, 'TanBeta': np.tan(Beta)}) #Define new mass in anyBSM

for i,x in enumerate(M_array):
    THDM2.setparameters({'M': x})
    lamb = THDM2.lambdahhh()
    Γ[i] = np.real(lamb["total"])
    Γ0[i] = np.real(lamb["treelevel"])

fig, ax = plt.subplots(figsize=(10, 10))

plt.plot(M_array, Γ/lambdahhhSM,c='k',label='AnyBSM')
plt.plot(M_array, Γ0/lambdahhhSM,c='b',label='AnyBSM (tree-level)')
plt.plot(M_array, -Gammahhh_treelevel(np.arcsin(sBmA)-np.pi/2, M_array)/lambdahhhSM,c='y',label='Kanemura (tree-level) (Approx. x)',ls='dashed')
plt.plot(M_array, -Gammahhh_oneloop(np.arcsin(sBmA)-np.pi/2, M_array,mH,mA,mHpm)/lambdahhhSM,c='orange',label='Kanemura (one-loop) (Approx. x)',ls='dashed')
plt.plot(M_array, -Gammahhh_treelevel_cos(BmA,Beta, M_array)/lambdahhhSM,c='y',label='Kanemura (tree-level)',ls='dashdot')
plt.plot(M_array, -Gammahhh_oneloop_cos_kan(BmA,Beta, M_array,mH,mA,mHpm)/lambdahhhSM,c='green',label='Kanemura (one-loop)')
plt.plot(M_array, -Gammahhh_oneloop_cos(BmA,Beta, M_array,mH,mA,mHpm)/lambdahhhSM,c='orange',ls='dashdot',label='Kanemura (Top) (one-loop)')
plt.plot(M_array, -Gammahhh_oneloop_cos_correc(BmA,Beta, M_array,mH,mA,mHpm)/lambdahhhSM,c='red',ls='dashed',label='Kanemura (Top+Higgs). (one-loop)')
plt.yticks(size=20)
plt.ylabel(r'$\kappa_\lambda$', size=25)
plt.xticks(size=20)
plt.xlabel(r'$M$ [GeV]', size=25)
ax.grid()

textstr = '\n'.join((
    r'\textbf{2HDM-II}',
    r'$M_A=%d$' %(mA),
    r'$M_{H}=%d$' %(mH),
    r'$M_{H^\pm}=%d$' %(mHpm),
#    r'$M=%d$' %(M),
    r'$\sin(\beta-\alpha)=%.2f$' %(sBmA),
    r'$\tan(\beta)=%.1f$' %(np.tan(Beta)),
    ))

ax.yaxis.get_offset_text().set_fontsize(20)
ax.text(0.05, 0.75, textstr, transform=ax.transAxes, fontsize=18, verticalalignment='top')


plt.legend(fontsize=20)

plt.show()

#%%

###                         Plot 4 - 2HDM 1-loop comparison as function of sβα

#Standard definitions in anyBSM
mA = 300 #Pseudoscalar mass
mHpm = 300 #Charged scalar mass
mH = 300  #Heavy Higgs mass
sBmA = 0.9 #Standard definition of sin(β-α) in anyBSM
M = 300 #Standard definition of M in anyBSM
BmA = np.arcsin(sBmA)
Beta = np.arctan(10)

sBmA_array = np.linspace(0,1,20,endpoint=False)
Γ = np.copy(sBmA_array)
Γ0 = np.copy(sBmA_array)
THDM2 = anyBSM('THDMII', scheme_name = 'OS')
THDM2.setparameters({'MAh2': mA, 'MHm2': mHpm, 'Mh2': mH, 'M': M, 'SinBmA': sBmA, 'TanBeta': np.tan(Beta)}) #Define new mass in anyBSM

for i,x in enumerate(sBmA_array):
    THDM2.setparameters({'SinBmA': x})
    lamb = THDM2.lambdahhh()
    Γ[i] = np.real(lamb["total"])
    Γ0[i] = np.real(lamb["treelevel"])

fig, ax = plt.subplots(figsize=(10, 10))

plt.plot(sBmA_array, Γ0/lambdahhhSM,c='b',label='AnyBSM (tree-level)')
plt.plot(sBmA_array, Γ/lambdahhhSM,c='k',label='AnyBSM')
#plt.plot(sBmA_array, -Gammahhh_treelevel(np.arcsin(sBmA_array)-np.pi/2, M)/lambdahhhSM,c='y',label='Kanemura (tree-level) Approx',ls='dashed')
#plt.plot(sBmA_array, -Gammahhh_oneloop(np.arcsin(sBmA_array)-np.pi/2, M,mH,mA,mHpm)/lambdahhhSM,c='orange',label='Kanemura (one-loop) Approx',ls='dashed')
plt.plot(sBmA_array, -Gammahhh_treelevel_cos(np.arcsin(sBmA_array),Beta, M)/lambdahhhSM,c='y',label='Kanemura (tree-level)',ls='dashdot')
plt.plot(sBmA_array, -Gammahhh_oneloop_cos_kan(np.arcsin(sBmA_array),Beta, M,mH,mA,mHpm)/lambdahhhSM,c='green',ls='dashed',label='Kanemura (one-loop)')
plt.plot(sBmA_array, -Gammahhh_oneloop_cos(np.arcsin(sBmA_array),Beta, M,mH,mA,mHpm)/lambdahhhSM,c='orange',label='Kanemura (Top) (one-loop)')
plt.plot(sBmA_array, -Gammahhh_oneloop_cos_correc(np.arcsin(sBmA_array),Beta, M,mH,mA,mHpm)/lambdahhhSM,c='red',ls='dashed',label='Kanemura (Top+Higgs). (one-loop)')
plt.yticks(size=20)
plt.ylabel(r'$\kappa_\lambda$', size=25)
plt.xticks(size=20)
plt.xlabel(r'$s_{\beta-\alpha}$', size=25)
ax.grid()

textstr = '\n'.join((
    r'\textbf{2HDM-II}',
    r'$M_A=%d$' %(mA),
    r'$M_{H}=%d$' %(mH),
    r'$M_{H^\pm}=%d$' %(mHpm),
    r'$M=%d$' %(M),
#    r'$\sin(\beta-\alpha)=%.2f$' %(sBmA),
    r'$\tan(\beta)=%.1f$' %(np.tan(Beta)),
    ))

ax.yaxis.get_offset_text().set_fontsize(20)
ax.text(0.05, 0.35, textstr, transform=ax.transAxes, fontsize=18, verticalalignment='top')

plt.legend(fontsize=20)

plt.show()

#%%%

###                         Plot 5 - 2HDM 1-loop comparison as function of tanb

#Standard definitions in anyBSM
mA = 300 #Pseudoscalar mass
mHpm = 300 #Charged scalar mass
mH = 300  #Heavy Higgs mass
sBmA = 0.9 #Standard definition of sin(β-α) in anyBSM
M = 800 #Standard definition of M in anyBSM
BmA = np.arcsin(sBmA)
Beta = np.arctan(2)


sBmA_array = np.linspace(0.8,40,40)
Γ = np.copy(sBmA_array)
Γ0 = np.copy(sBmA_array)
THDM2 = anyBSM('THDMII', scheme_name = 'OSalignment')
THDM2.setparameters({'MAh2': mA, 'MHm2': mHpm, 'Mh2': mH, 'M': M, 'SinBmA': sBmA, 'TanBeta': np.tan(Beta)}) #Define new mass in anyBSM

for i,x in enumerate(sBmA_array):
    THDM2.setparameters({'TanBeta': x})
    lamb = THDM2.lambdahhh()
    Γ[i] = np.real(lamb["total"])
    Γ0[i] = np.real(lamb["treelevel"])

fig, ax = plt.subplots(figsize=(10, 10))

plt.plot(sBmA_array, Γ0/lambdahhhSM,c='b',label='AnyBSM (tree-level)')
plt.plot(sBmA_array, Γ/lambdahhhSM,c='k',label='AnyBSM')
#plt.plot(sBmA_array, -Gammahhh_treelevel(np.arcsin(sBmA_array)-np.pi/2, M)/lambdahhhSM,c='y',label='Kanemura (tree-level) Approx',ls='dashed')
#plt.plot(sBmA_array, -Gammahhh_oneloop(np.arcsin(sBmA_array)-np.pi/2, M,mH,mA,mHpm)/lambdahhhSM,c='orange',label='Kanemura (one-loop) Approx',ls='dashed')
plt.plot(sBmA_array, -Gammahhh_treelevel_cos(np.arcsin(sBmA),np.arctan(sBmA_array), M)/lambdahhhSM,c='y',label='Kanemura (tree-level)',ls='dashdot')
plt.plot(sBmA_array, -Gammahhh_oneloop_cos_kan(np.arcsin(sBmA),np.arctan(sBmA_array), M,mH,mA,mHpm)/lambdahhhSM,c='green',ls='dashed',label='Kanemura (one-loop)')
plt.plot(sBmA_array, -Gammahhh_oneloop_cos(np.arcsin(sBmA),np.arctan(sBmA_array), M,mH,mA,mHpm)/lambdahhhSM,c='orange',label='Kanemura (Top) (one-loop)')
plt.plot(sBmA_array, -Gammahhh_oneloop_cos_correc(np.arcsin(sBmA),np.arctan(sBmA_array), M,mH,mA,mHpm)/lambdahhhSM,c='red',ls='dashed',label='Kanemura (Top+Higgs). (one-loop)')
plt.yticks(size=20)
plt.ylabel(r'$\kappa_\lambda$', size=25)
plt.xticks(size=20)
plt.xlabel(r'$\tan{\beta}$', size=25)
ax.grid()

textstr = '\n'.join((
    r'\textbf{2HDM-II}',
    r'$M_A=%d$' %(mA),
    r'$M_{H}=%d$' %(mH),
    r'$M_{H^\pm}=%d$' %(mHpm),
    r'$M=%d$' %(M),
    r'$\sin(\beta-\alpha)=%.2f$' %(sBmA),
#    r'$\tan(\beta)=%.1f$' %(np.tan(Beta)),
    ))

ax.yaxis.get_offset_text().set_fontsize(20)
ax.text(0.05, 0.6, textstr, transform=ax.transAxes, fontsize=18, verticalalignment='top')

plt.legend(fontsize=20,loc='lower left')

plt.show()

#%%%

###                         Plot 5 - 2HDM 1-loop comparison as function of tanb - Models

#Standard definitions in anyBSM
mA = 500 #Pseudoscalar mass
mHpm = 400 #Charged scalar mass
mH = 100  #Heavy Higgs mass
sBmA = 0.9 #Standard definition of sin(β-α) in anyBSM
M = 800 #Standard definition of M in anyBSM
BmA = np.arcsin(sBmA)
Beta = np.arctan(2)


sBmA_array = np.linspace(0.8,40,40)
Γ = np.copy(sBmA_array)
Γ2 = np.copy(sBmA_array)
Γ0 = np.copy(sBmA_array)
THDM1 = anyBSM('THDMI', scheme_name = 'OSalignment')
THDM2 = anyBSM('THDMII', scheme_name = 'OSalignment')
THDM1.setparameters({'MAh2': mA, 'MHm2': mHpm, 'Mh2': mH, 'M': M, 'SinBmA': sBmA, 'TanBeta': np.tan(Beta)}) #Define new mass in anyBSM
THDM2.setparameters({'MAh2': mA, 'MHm2': mHpm, 'Mh2': mH, 'M': M, 'SinBmA': sBmA, 'TanBeta': np.tan(Beta)}) #Define new mass in anyBSM

for i,x in enumerate(sBmA_array):
    THDM1.setparameters({'TanBeta': x})
    lamb1 = THDM1.lambdahhh()
    lamb2 = THDM2.lambdahhh()
    Γ[i] = np.real(lamb1["total"])
    Γ2[i] = np.real(lamb2["total"])
    Γ0[i] = np.real(lamb1["treelevel"])

fig, ax = plt.subplots(figsize=(10, 10))

plt.plot(sBmA_array, Γ0/lambdahhhSM,c='y',label='Tree-level')
plt.plot(sBmA_array, Γ/lambdahhhSM,c='k',label='Type I')
plt.plot(sBmA_array, Γ/lambdahhhSM,c='r',label='Type II')


plt.yticks(size=20)
plt.ylabel(r'$\kappa_\lambda$', size=25)
plt.xticks(size=20)
plt.xlabel(r'$\tan{\beta}$', size=25)
ax.grid()

textstr = '\n'.join((
    r'\textbf{2HDM-II}',
    r'$M_A=%d$' %(mA),
    r'$M_{H}=%d$' %(mH),
    r'$M_{H^\pm}=%d$' %(mHpm),
    r'$M=%d$' %(M),
    r'$\sin(\beta-\alpha)=%.2f$' %(sBmA),
#    r'$\tan(\beta)=%.1f$' %(np.tan(Beta)),
    ))

ax.yaxis.get_offset_text().set_fontsize(20)
ax.text(0.05, 0.25, textstr, transform=ax.transAxes, fontsize=18, verticalalignment='top')

plt.legend(fontsize=20,loc='lower right')

plt.show()

#%%                         Plots 2HDM

#Standard definitions in anyBSM
mA = 300 #Pseudoscalar mass
mHpm = 300 #Charged scalar mass
mH = 300  #Heavy Higgs mass
sBmA = 1 #Standard definition of sin(β-α) in anyBSM
M = 300 #Standard definition of M in anyBSM
BmA = np.arcsin(sBmA)
Beta = np.arctan(1)

###                         Plot 3 - 2HDM 1-loop comparison as function of M_A

M_array = np.linspace(0,500,100)
Γ = np.copy(M_array)
Γ0 = np.copy(M_array)
if sBmA == 1:
    scheme_str = 'OSalignment'
else:
    scheme_str = 'OS'
    

THDM2 = anyBSM('THDMI', scheme_name = scheme_str)
THDM2.setparameters({'MAh2': mA, 'MHm2': mHpm, 'Mh2': mH, 'M': M, 'SinBmA': sBmA, 'TanBeta': np.tan(Beta)}) #Define new mass in anyBSM

for i,x in enumerate(M_array):
    THDM2.setparameters({'M': x})
    lamb = THDM2.lambdahhh()
    Γ[i] = np.real(lamb["total"])
    Γ0[i] = np.real(lamb["treelevel"])

fig, ax = plt.subplots(figsize=(10, 10))

plt.plot(M_array, Γ/lambdahhhSM,c='k',label='AnyBSM')
plt.plot(M_array, Γ0/lambdahhhSM,c='b',label='AnyBSM (tree-level)')
plt.plot(M_array, -Gammahhh_treelevel(np.arcsin(sBmA)-np.pi/2, M_array)/lambdahhhSM,c='y',label='Kanemura (tree-level) (Approx. x)',ls='dashed')
plt.plot(M_array, -Gammahhh_oneloop(np.arcsin(sBmA)-np.pi/2, M_array,mH,mA,mHpm)/lambdahhhSM,c='orange',label='Kanemura (one-loop) (Approx. x)',ls='dashed')
plt.plot(M_array, -Gammahhh_treelevel_cos(BmA,Beta, M_array)/lambdahhhSM,c='y',label='Kanemura (tree-level)',ls='dashdot')
plt.plot(M_array, -Gammahhh_oneloop_cos_kan(BmA,Beta, M_array,mH,mA,mHpm)/lambdahhhSM,c='green',label='Kanemura (one-loop)')
plt.plot(M_array, -Gammahhh_oneloop_cos(BmA,Beta, M_array,mH,mA,mHpm)/lambdahhhSM,c='orange',ls='dashdot',label='Kanemura (Top) (one-loop)')
plt.plot(M_array, -Gammahhh_oneloop_cos_correc(BmA,Beta, M_array,mH,mA,mHpm)/lambdahhhSM,c='red',ls='dashed',label='Kanemura (Top+Higgs). (one-loop)')
plt.yticks(size=20)
plt.ylabel(r'$\kappa_\lambda$', size=25)
plt.xticks(size=20)
plt.xlabel(r'$M$ [GeV]', size=25)
ax.grid()

textstr = '\n'.join((
    r'\textbf{2HDM-II}',
    r'$M_A=%d$' %(mA),
    r'$M_{H}=%d$' %(mH),
    r'$M_{H^\pm}=%d$' %(mHpm),
#    r'$M=%d$' %(M),
    r'$\sin(\beta-\alpha)=%.2f$' %(sBmA),
    r'$\tan(\beta)=%.1f$' %(np.tan(Beta)),
    ))

ax.yaxis.get_offset_text().set_fontsize(20)
ax.text(0.05, 0.75, textstr, transform=ax.transAxes, fontsize=18, verticalalignment='top')


plt.legend(fontsize=20)

plt.show()

#%%                         Plots 2HDM

#Standard definitions in anyBSM
mA = 300 #Pseudoscalar mass
mHpm = 300 #Charged scalar mass
mH = 300  #Heavy Higgs mass
sBmA = 0.9 #Standard definition of sin(β-α) in anyBSM
M = 800 #Standard definition of M in anyBSM
BmA = np.arcsin(sBmA)
Beta = np.arctan(10)

###                         Plot 3 - 2HDM 1-loop comparison as function of M_A

M_array = np.linspace(200,800,100)
Γ = np.copy(M_array)
Γ0 = np.copy(M_array)
if sBmA == 1:
    scheme_str = 'OSalignment'
else:
    scheme_str = 'OS'
    

THDM2 = anyBSM('THDMII', scheme_name = scheme_str)
THDM2.setparameters({'MAh2': mA, 'MHm2': mHpm, 'Mh2': mH, 'M': M, 'SinBmA': sBmA, 'TanBeta': np.tan(Beta)}) #Define new mass in anyBSM

for i,x in enumerate(M_array):
    THDM2.setparameters({'MAh2': x})
    lamb = THDM2.lambdahhh()
    Γ[i] = np.real(lamb["total"])
    Γ0[i] = np.real(lamb["treelevel"])

fig, ax = plt.subplots(figsize=(10, 10))

plt.plot(M_array, Γ/lambdahhhSM,c='k',label='AnyBSM')
#plt.plot(M_array, Γ0/lambdahhhSM,c='b',label='AnyBSM (tree-level)')
#plt.plot(M_array, -Gammahhh_treelevel(np.arcsin(sBmA)-np.pi/2, M)/lambdahhhSM,c='y',label='Kanemura (tree-level) (Approx. x)',ls='dashed')
plt.plot(M_array, -Gammahhh_oneloop(np.arcsin(sBmA)-np.pi/2, M,mH,M_array,mHpm)/lambdahhhSM,c='orange',label='Kanemura (one-loop) (Approx. x)',ls='dashed')
#plt.plot(M_array, -Gammahhh_treelevel_cos(BmA,Beta, M)/lambdahhhSM,c='y',label='Kanemura (tree-level)',ls='dashdot')
plt.plot(M_array, -Gammahhh_oneloop_cos_kan(BmA,Beta, M,mH,M_array,mHpm)/lambdahhhSM,c='green',label='Kanemura (one-loop)')
plt.plot(M_array, -Gammahhh_oneloop_cos(BmA,Beta, M,mH,M_array,mHpm)/lambdahhhSM,c='orange',ls='dashdot',label='Kanemura (Top) (one-loop)')
plt.plot(M_array, -Gammahhh_oneloop_cos_correc(BmA,Beta, M,mH,M_array,mHpm)/lambdahhhSM,c='red',ls='dashed',label='Kanemura (Top+Higgs). (one-loop)')
plt.yticks(size=20)
plt.ylabel(r'$\kappa_\lambda$', size=25)
plt.xticks(size=20)
plt.xlabel(r'$M_A$ [GeV]', size=25)
ax.grid()

#plt.ylim(-10,10)

textstr = '\n'.join((
    r'\textbf{2HDM-II}',
#    r'$M_A=%d$' %(mA),
    r'$M_{H}=%d$' %(mH),
    r'$M_{H^\pm}=%d$' %(mHpm),
    r'$M=%d$' %(M),
    r'$\sin(\beta-\alpha)=%.2f$' %(sBmA),
    r'$\tan(\beta)=%.1f$' %(np.tan(Beta)),
    ))

ax.yaxis.get_offset_text().set_fontsize(20)
ax.text(0.05, 0.85, textstr, transform=ax.transAxes, fontsize=18, verticalalignment='top')


plt.legend(fontsize=20)

plt.show()

#%%                         Plots 2HDM - 331

#Standard definitions in anyBSM
mA = 300 #Pseudoscalar mass
mHpm = 300 #Charged scalar mass
mH = 300  #Heavy Higgs mass
sBmA = 0.9 #Standard definition of sin(β-α) in anyBSM
M = 300 #Standard definition of M in anyBSM
BmA = np.arcsin(sBmA)
Beta = np.arctan(1)

###                         Plot 3 - 2HDM 1-loop comparison as function of M_A

M_array = np.linspace(200,800,100)
Γ = np.copy(M_array)
Γ0 = np.copy(M_array)
if sBmA == 1:
    scheme_str = 'OSalignment'
else:
    scheme_str = 'OS'
    

THDM2 = anyBSM('THDMI', scheme_name = scheme_str)
THDM2.setparameters({'MAh2': mA, 'MHm2': mHpm, 'Mh2': mH, 'M': M, 'SinBmA': sBmA, 'TanBeta': np.tan(Beta)}) #Define new mass in anyBSM

for i,x in enumerate(M_array):
    THDM2.setparameters({'MAh2': x, 'M': x})
    lamb = THDM2.lambdahhh()
    Γ[i] = np.real(lamb["total"])
    Γ0[i] = np.real(lamb["treelevel"])

fig, ax = plt.subplots(figsize=(10, 10))

plt.plot(M_array, Γ/lambdahhhSM,c='k',label='AnyBSM')
#plt.plot(M_array, Γ0/lambdahhhSM,c='b',label='AnyBSM (tree-level)')
#plt.plot(M_array, -Gammahhh_treelevel(np.arcsin(sBmA)-np.pi/2, M)/lambdahhhSM,c='y',label='Kanemura (tree-level) (Approx. x)',ls='dashed')
plt.plot(M_array, -Gammahhh_oneloop(np.arcsin(sBmA)-np.pi/2, M_array,mH,M_array,mHpm)/lambdahhhSM,c='orange',label='Kanemura (one-loop) (Approx. x)',ls='dashed')
#plt.plot(M_array, -Gammahhh_treelevel_cos(BmA,Beta, M)/lambdahhhSM,c='y',label='Kanemura (tree-level)',ls='dashdot')
plt.plot(M_array, -Gammahhh_oneloop_cos_kan(BmA,Beta, M_array,mH,M_array,mHpm)/lambdahhhSM,c='green',label='Kanemura (one-loop)')
plt.plot(M_array, -Gammahhh_oneloop_cos(BmA,Beta, M_array,mH,M_array,mHpm)/lambdahhhSM,c='orange',ls='dashdot',label='Kanemura (Top) (one-loop)')
plt.plot(M_array, -Gammahhh_oneloop_cos_correc(BmA,Beta, M_array,mH,M_array,mHpm)/lambdahhhSM,c='red',ls='dashed',label='Kanemura (Top+Higgs). (one-loop)')
plt.yticks(size=20)
plt.ylabel(r'$\kappa_\lambda$', size=25)
plt.xticks(size=20)
plt.xlabel(r'$M_A$ [GeV]', size=25)
ax.grid()

#plt.ylim(-10,10)

textstr = '\n'.join((
    r'\textbf{2HDM-331 EFT}',
#    r'$M_A=%d$' %(mA),
    r'$M_{H}=%d$' %(mH),
    r'$M_{H^\pm}=%d$' %(mHpm),
#    r'$M=%d$' %(M),
    r'$\sin(\beta-\alpha)=%.2f$' %(sBmA),
    r'$\tan(\beta)=%.1f$' %(np.tan(Beta)),
    ))

ax.yaxis.get_offset_text().set_fontsize(20)
ax.text(0.05, 0.55, textstr, transform=ax.transAxes, fontsize=18, verticalalignment='top')


plt.legend(fontsize=20)

plt.show()
