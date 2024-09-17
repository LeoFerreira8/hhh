#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Oct 23 17:53:30 2023

@author: leonardoferreira
"""

from anyBSM import anyBSM
import numpy as np
import matplotlib.pyplot as plt
import scan_parameterspace as spr
import scan_parameterspace_funcs as fcs
import quartic_couplings as qtcp
from mpl_toolkits.axes_grid1.inset_locator import zoomed_inset_axes 
from mpl_toolkits.axes_grid1.inset_locator import mark_inset


SM = anyBSM('SM',scheme_name = 'OS')
THDM1 = anyBSM('THDMI', scheme_name = 'OS')
THDM2 = anyBSM('THDMII', scheme_name = 'OS')

THDM2.setparameters({'MAh2': 200, 'MHm2': 200, 'Mh2': 200, 'M': 200, 'SinBmA': 0.9, 'TanBeta': 40, 'Mh1': 125.6}) #Define new mass in anyBSM

#lamb = SM.lambdahhh()
#lamb1 = THDM1.lambdahhh()
#THDM2.set_evaluation_mode('analytical')
#lamb2 = THDM2.lambdahhh(draw=0.50,simplify=True)
#expre = THDM2.SolveDependencies(lamb2["total"])
#lamb2 = THDM2.lambdahhh()

plt.rc('text', usetex=True)
plt.rc('xtick',labelsize=40)
plt.rc('ytick',labelsize=40)
resol=400

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

def Gammahhh_oneloop_cos_correc_coup(ba,b,M,mH,mA,mHpm):
    R47 = (np.sin(ba)*(2*mA**2-mh**2)+((np.cos(ba)*(1-np.tan(b)**2))/(2*np.tan(b))+np.sin(ba))*(2*mh**2-(M**2-mA**2)))/(mA**2+mh**2-M**2)
    G = Gammahhh_treelevel_cos(ba, b, M)-(3*mh**2/v)*(-(np.cos(ba)/np.tan(b)+np.sin(ba))**3*3*mt**4/(3*np.pi**2*mh**2*v**2)+mh**4/(2*np.pi**2*mh**2*v**2)+(mH**4/(12*np.pi**2*mh**2*v**2))*(1-M**2/(mH**2))**3+R47**3*(mA**4/(12*np.pi**2*mh**2*v**2))*(1-M**2/(mA**2))**3+(mHpm**4/(6*np.pi**2*mh**2*v**2))*(1-M**2/(mHpm**2))**3)
    #G = Gammahhh_treelevel_cos(ba, b, M)-(3*mh**2/v)*(-3*mt**4/(3*np.pi**2*mh**2*v**2)+(mH**4/(12*np.pi**2*mh**2*v**2))*(1-M**2/(mH**2))**3+(mA**4/(12*np.pi**2*mh**2*v**2))*(1-M**2/(mA**2))**3+(mHpm**4/(6*np.pi**2*mh**2*v**2))*(1-M**2/(mHpm**2))**3)
    
    return G

def Gammahhh_oneloop_cos_correc_coup_full(ba,b,M,mH,mA,mHpm):
    a = b-ba
    chH = -(np.sin(ba)/np.sin(2*b)*np.sin(2*a)*(mh**2 + 2*mH**2) - M**2*(3*np.sin(2*a) + np.sin(2*b)) + mh**2*np.sin(ba))/mH**2/2
    chHp = -(np.sin(ba)*(mh**2 - 2*mHpm**2) - np.cos(a+b)/np.sin(2*b)*(2*mh**2 - 2*M**2) + mh**2)/2/mHpm**2
    chA = (-np.sin(b - a)*(mh**2 - 2*mA**2) + np.cos(a + b)/np.sin(2*b)*(2*mh**2 - 2*M**2) - mh**2*np.sin(b - a))/2/mA**2
    ccH = 1-M**2/(mH**2)+((2*mH**2+mh**2-3*M**2)/(2*mH**2))*(np.tan(b)-1/np.tan(b))*np.arccos(np.sin(ba))
    ccA = 1-M**2/(mA**2)+((M**2-mh**2)/(2*mA**2))*(np.tan(b)-1/np.tan(b))*np.arccos(np.sin(ba))
    ccHp = 1-M**2/(mHpm**2)+((M**2-mh**2)/(2*mHpm**2))*(np.tan(b)-1/np.tan(b))*np.arccos(np.sin(ba))
    #R47 = (np.sin(ba)*(2*mA**2-mh**2)+((np.cos(ba)*(1-np.tan(b)**2))/(2*np.tan(b))+np.sin(ba))*(2*mh**2-(M**2-mA**2)))/(mA**2+mh**2-M**2)
    G = Gammahhh_treelevel_cos(ba, b, M)-(3*mh**2/(v))*(-(np.cos(ba)/np.tan(b)+np.sin(ba))**3*3*mt**4/(3*np.pi**2*mh**2*v**2)+mh**4/(2*np.pi**2*mh**2*v**2)+(mH**4/(12*np.pi**2*mh**2*v**2))*(chH)**3+(mA**4/(12*np.pi**2*mh**2*v**2))*(chA)**3+(mHpm**4/(6*np.pi**2*mh**2*v**2))*(chHp)**3)
    G = Gammahhh_treelevel_cos(ba, b, M)-(3*mh**2/(v))*(-(np.cos(ba)/np.tan(b)+np.sin(ba))**3*3*mt**4/(3*np.pi**2*mh**2*v**2)+mh**4/(2*np.pi**2*mh**2*v**2)+(mH**4/(12*np.pi**2*mh**2*v**2))*(ccH)**3+(mA**4/(12*np.pi**2*mh**2*v**2))*(ccA)**3+(mHpm**4/(6*np.pi**2*mh**2*v**2))*(ccHp)**3)
    #G = Gammahhh_treelevel_cos(ba, b, M)*(1-(-(np.cos(ba)/np.tan(b)+np.sin(ba))**3*3*mt**4/(3*np.pi**2*mh**2*v**2)+mh**4/(2*np.pi**2*mh**2*v**2)+(mH**4/(12*np.pi**2*mh**2*v**2))*(chH)**3+(mA**4/(12*np.pi**2*mh**2*v**2))*(chA)**3+(mHpm**4/(6*np.pi**2*mh**2*v**2))*(chHp)**3))
    #G = Gammahhh_treelevel_cos(ba, b, M)-(3*mh**2/(v*np.sin(2*b)))*(-(np.cos(ba)/np.tan(b)+np.sin(ba))**3*3*mt**4/(3*np.pi**2*mh**2*v**2)+mh**4/(2*np.pi**2*mh**2*v**2)+(mH**4/(12*np.pi**2*mh**2*v**2))*(chH)**3+(mA**4/(12*np.pi**2*mh**2*v**2))*(chA)**3+(mHpm**4/(6*np.pi**2*mh**2*v**2))*(chHp)**3)
    
    return G

# #%% Plots para checar renormalização

# ###                         Plot 1 - SM comparison with 1-loop as function of mh

# mh_array = np.linspace(-50,100,101)+mhSM #Plot mhSM +-50 GeV
# δΓ = np.copy(mh_array)
# δΓ2 = np.copy(mh_array)

# for i,x in enumerate(mh_array):
#     SM.setparameters({'Mh': x}) #Define new mass in anyBSM
#     lamb = SM.lambdahhh()  #Recalculate lambda
#     mh = x  #Redefine mh globally to recalculate Gamma
#     kappa_kan = (-Gammahhh_oneloop_SM_like(0, 0))/(-Gammahhh_treelevel(0, 0))
#     kappa_kan_full = (-Gammahhh_oneloop_SM_like_full(0, 0))/(-Gammahhh_treelevel(0, 0))
#     kappa_any = np.real(lamb["total"])/(-Gammahhh_treelevel(0, 0))
#     δΓ[i] = (kappa_any-kappa_kan)/kappa_kan
#     δΓ2[i] = (kappa_any-kappa_kan_full)/kappa_kan_full


# fig, ax = plt.subplots(figsize=(10, 10))

# plt.plot(mh_array, δΓ,label='Without Higgs loop correction')
# plt.plot(mh_array, δΓ2,label='With Higgs loop correction')
# plt.yticks(size=20)
# plt.ylabel(r'$(\kappa_\lambda^{{Kan}}-\kappa_\lambda)/\kappa_\lambda^{{Kan}}$', size=25)
# plt.xticks(size=20)
# plt.xlabel(r'SM-Higgs mass $m_h$ [GeV]', size=25)
# ax.grid()

# plt.legend(fontsize=20)

# plt.show()

# mh = mhSM
# SM = anyBSM('SM',scheme_name = 'OS')


#%% Plots para checar higgs corrections

###                         Plot 1 - THDM comparison with 1-loop as function of mh

mh_array = np.linspace(-50,50,101)+mhSM #Plot mhSM +-50 GeV
δΓ = np.copy(mh_array)
δΓ2 = np.copy(mh_array)

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

THDM2 = anyBSM('THDMII', scheme_name = scheme_str)
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


fig, ax = plt.subplots(figsize=(10, 10),dpi=resol)

plt.plot(mh_array, δΓ,label='Without Higgs loop correction')
plt.plot(mh_array, δΓ2,label='With Higgs loop correction')
plt.yticks(size=25)
plt.ylabel(r'$(\bar{\kappa}_\lambda-\kappa_\lambda)/\bar{\kappa}_\lambda$', size=25)
plt.xticks(size=25)
plt.xlabel(r'SM-Higgs mass $m_h$ [GeV]', size=25)
plt.minorticks_on()
plt.tick_params(axis='both', which='both', right=True, top=True, direction='in')
ax.grid()

textstr = '\n'.join((
    r'\textbf{2HDM-II}',
    r'$m_A=%d$ GeV' %(mA),
    r'$m_{H}=%d$ GeV' %(mH),
    r'$m_{H^\pm}=%d$ GeV' %(mHpm),
    r'$M=%d$ GeV' %(M),
    r'$\sin(\beta-\alpha)=%.2f$' %(sBmA),
    r'$\tan(\beta)=%.1f$' %(np.tan(Beta)),
    ))

ax.text(0.70, 0.5, textstr, transform=ax.transAxes, fontsize=20, verticalalignment='top')

plt.legend(fontsize=25,loc='upper left')

plt.show()

mh = mhSM

# #%%

# ###                         Plot 2 - SM comparison with 1-loop as function of mt

# mt_array = np.linspace(-50,450,101)+mtSM
# δΓ = np.copy(mt_array)

# for i,x in enumerate(mt_array):
#     SM.setparameters({'Mu3': x})
#     lamb = SM.lambdahhh()
#     mt = x
#     δΓ[i] = (np.real(lamb["total"])-(-Gammahhh_oneloop_SM_like(0, 0)))/(-Gammahhh_treelevel(0, 0))


# fig, ax = plt.subplots(figsize=(10, 10))

# plt.plot(mt_array, δΓ)
# plt.yticks(size=20)
# plt.ylabel(r'$\delta \Gamma_{(1)} / \Gamma_0$', size=25)
# plt.xticks(size=20)
# plt.xlabel(r'Top-quark mass $m_t$ [GeV]', size=25)
# ax.grid()
    

# plt.show()

# mt = mtSM
# SM = anyBSM('SM',scheme_name = 'OS')

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

M_array = np.linspace(0,500,300)
Γ = np.copy(M_array)
Γ0 = np.copy(M_array)
Γ1 = np.copy(M_array)
Γ2 = np.copy(M_array)
Γ3 = np.copy(M_array)
if sBmA == 1:
    scheme_str = 'OSalignment'
else:
    scheme_str = 'OS'
    

THDM2 = anyBSM('THDMII', scheme_name = scheme_str)
THDM2.setparameters({'MAh2': mA, 'MHm2': mHpm, 'Mh2': mH, 'M': M, 'SinBmA': sBmA, 'TanBeta': np.tan(Beta)}) #Define new mass in anyBSM

for i,x in enumerate(M_array):
    THDM2.setparameters({'M': x})
    lamb = THDM2.lambdahhh()
    Γ[i] = np.real(lamb["total"])/lambdahhhSM
    Γ0[i] = np.real(lamb["treelevel"])/lambdahhhSM
    Γ1[i] = (-Gammahhh_oneloop_cos_kan(BmA,Beta, x,mH,mA,mHpm)/lambdahhhSM)
    Γ2[i] = (-Gammahhh_oneloop_cos(BmA,Beta, x,mH,mA,mHpm)/lambdahhhSM)
    Γ3[i] = (-Gammahhh_oneloop_cos_correc(BmA,Beta, x,mH,mA,mHpm)/lambdahhhSM)

#%%

fig, ax = plt.subplots(figsize=(10, 10),dpi=resol)

plt.plot(M_array, Γ,c='k',label='Full calculation')
#plt.plot(M_array, (Γ0/lambdahhhSM),c='b',label='AnyBSM (tree-level)')
#plt.plot(M_array, (-Gammahhh_treelevel(np.arcsin(sBmA)-np.pi/2, M_array)/lambdahhhSM),c='y',label='Kanemura (tree-level) (Approx. x)',ls='dashed')
#plt.plot(M_array, (-Gammahhh_oneloop(np.arcsin(sBmA)-np.pi/2, M_array,mH,mA,mHpm)/lambdahhhSM),c='orange',label='Kanemura (one-loop) (Approx. x)',ls='dashed')
#plt.plot(M_array, (-Gammahhh_treelevel_cos(BmA,Beta, M_array)/lambdahhhSM),c='y',label='Kanemura (tree-level)',ls='dashdot')
plt.plot(M_array, Γ1, c='green',label='Leading contribution')
plt.plot(M_array, Γ2,c='orange',ls='dashdot',label='Leading contribution (+Top)')
plt.plot(M_array, Γ3,c='red',ls='dashed',label='Leading contribution (+Top+Higgs)')
plt.yticks(size=25)
plt.ylabel(r'$\kappa_\lambda$', size=25)
plt.xticks(size=25)
plt.xlabel(r'$M$ [GeV]', size=25)
plt.minorticks_on()
plt.tick_params(axis='both', which='both', right=True, top=True, direction='in')
#plt.yscale('log')
#plt.ylim(0.9,1.1)

ax.grid()

textstr = '\n'.join((
    r'\textbf{2HDM-II}',
    r'$m_A=%d$ GeV' %(mA),
    r'$m_{H}=%d$ GeV' %(mH),
    r'$m_{H^\pm}=%d$ GeV' %(mHpm),
#    r'$M=%d$' %(M),
    r'$\sin(\beta-\alpha)=%.2f$' %(sBmA),
    r'$\tan(\beta)=%.1f$' %(np.tan(Beta)),
    ))

ax.yaxis.get_offset_text().set_fontsize(20)
ax.text(0.05, 0.69, textstr, transform=ax.transAxes, fontsize=20, verticalalignment='top')


plt.legend(loc='lower left',fontsize=18)

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
Beta = np.arctan(1)

sBmA_array = np.linspace(0,1,200)
Γ = np.copy(sBmA_array)
Γ0 = np.copy(sBmA_array)
THDM2 = anyBSM('THDMII', scheme_name = 'OS')
THDM2.setparameters({'MAh2': mA, 'MHm2': mHpm, 'Mh2': mH, 'M': M, 'SinBmA': sBmA, 'TanBeta': np.tan(Beta)}) #Define new mass in anyBSM

for i,x in enumerate(sBmA_array):
    THDM2.setparameters({'SinBmA': x})
    if x==1:
        THDM2.load_renormalization_scheme('OSalignment')
    lamb = THDM2.lambdahhh()
    Γ[i] = np.real(lamb["total"])
    Γ0[i] = np.real(lamb["treelevel"])
    
#%%

fig, ax = plt.subplots(figsize=(10, 10),dpi=resol)

plt.plot(sBmA_array, Γ0/lambdahhhSM,c='b',label='Full calculation (tree-level)')
plt.plot(sBmA_array, Γ/lambdahhhSM,c='k',label='Full calculation')
#plt.plot(sBmA_array, -Gammahhh_treelevel(np.arcsin(sBmA_array)-np.pi/2, M)/lambdahhhSM,c='y',label='Kanemura (tree-level) Approx',ls='dashed')
#plt.plot(sBmA_array, -Gammahhh_oneloop(np.arcsin(sBmA_array)-np.pi/2, M,mH,mA,mHpm)/lambdahhhSM,c='orange',label='Kanemura (one-loop) Approx',ls='dashed')
plt.plot(sBmA_array, -Gammahhh_treelevel_cos(np.arcsin(sBmA_array),Beta, M)/lambdahhhSM,c='y',label='Leading contribution (tree-level)',ls='dashdot')
plt.plot(sBmA_array, -Gammahhh_oneloop_cos_kan(np.arcsin(sBmA_array),Beta, M,mH,mA,mHpm)/lambdahhhSM,c='green',ls='dashed',label='Leading contribution')
plt.plot(sBmA_array, -Gammahhh_oneloop_cos(np.arcsin(sBmA_array),Beta, M,mH,mA,mHpm)/lambdahhhSM,c='orange',label='Leading contribution (+Top)')
plt.plot(sBmA_array, -Gammahhh_oneloop_cos_correc(np.arcsin(sBmA_array),Beta, M,mH,mA,mHpm)/lambdahhhSM,c='red',ls='dashed',label='Leading contribution (+Top+Higgs)')
plt.yticks(size=25)
plt.ylabel(r'$\kappa_\lambda$', size=25)
plt.xticks(size=25)
plt.xlabel(r'$s_{\beta-\alpha}$', size=25)
ax.grid()
plt.minorticks_on()
plt.tick_params(axis='both', which='both', right=True, top=True, direction='in')

textstr = '\n'.join((
    r'\textbf{2HDM-II}',
    r'$m_A=%d$ GeV' %(mA),
    r'$m_{H}=%d$ GeV' %(mH),
    r'$m_{H^\pm}=%d$ GeV' %(mHpm),
    r'$M=%d$ GeV' %(M),
#    r'$\sin(\beta-\alpha)=%.2f$' %(sBmA),
    r'$\tan(\beta)=%.1f$' %(np.tan(Beta)),
    ))

ax.yaxis.get_offset_text().set_fontsize(20)
ax.text(0.05, 0.25, textstr, transform=ax.transAxes, fontsize=20, verticalalignment='top')

# Make the zoom-in plot:
    
x1 = 0.99
x2 = 1.00

# select y-range for zoomed region
y1 = 0.85
y2 = 1.005

axins = zoomed_inset_axes(ax, 28, loc='lower right', bbox_to_anchor=(0.67, 0.18, 0.5, 0.5), bbox_transform=ax.figure.transFigure,axes_kwargs={'aspect': .08}) # zoom = 10
axins.plot(sBmA_array, Γ0/lambdahhhSM,c='b')
axins.plot(sBmA_array, Γ/lambdahhhSM,c='k')
#plt.plot(sBmA_array, -Gammahhh_treelevel(np.arcsin(sBmA_array)-np.pi/2, M)/lambdahhhSM,c='y',label='Kanemura (tree-level) Approx',ls='dashed')
#plt.plot(sBmA_array, -Gammahhh_oneloop(np.arcsin(sBmA_array)-np.pi/2, M,mH,mA,mHpm)/lambdahhhSM,c='orange',label='Kanemura (one-loop) Approx',ls='dashed')
axins.plot(sBmA_array, -Gammahhh_treelevel_cos(np.arcsin(sBmA_array),Beta, M)/lambdahhhSM,c='y',ls='dashdot')
axins.plot(sBmA_array, -Gammahhh_oneloop_cos_kan(np.arcsin(sBmA_array),Beta, M,mH,mA,mHpm)/lambdahhhSM,c='green',ls='dashed')
axins.plot(sBmA_array, -Gammahhh_oneloop_cos(np.arcsin(sBmA_array),Beta, M,mH,mA,mHpm)/lambdahhhSM,c='orange')
axins.plot(sBmA_array, -Gammahhh_oneloop_cos_correc(np.arcsin(sBmA_array),Beta, M,mH,mA,mHpm)/lambdahhhSM,c='red',ls='dashed')
axins.set_xlim(x1, x2)
axins.set_ylim(y1, y2)
plt.xticks(visible=True,fontsize=20)
plt.yticks(visible=True,fontsize=20)
axins.grid()
mark_inset(ax, axins,loc1=2,loc2=1, fc="none", ec="0.65",clip_on=True)
plt.draw()
axins.minorticks_on()
plt.tick_params(axis='both', which='both', right=True, top=True, direction='in')


ax.legend(fontsize=18,framealpha=1,loc='upper center')

plt.savefig('../temp/kappa_sino.png', dpi=600, bbox_inches='tight')

plt.show()

#%%%

###                         Plot 5 - 2HDM 1-loop comparison as function of tanb

#Standard definitions in anyBSM
mA = 300 #Pseudoscalar mass
mHpm = 300 #Charged scalar mass
mH = 300  #Heavy Higgs mass
sBmA = 0.98 #Standard definition of sin(β-α) in anyBSM
M = 300 #Standard definition of M in anyBSM
BmA = np.arcsin(sBmA)
Beta = np.arctan(2)


sBmA_array = np.linspace(0.8,40,200)
Γ = np.copy(sBmA_array)
Γ0 = np.copy(sBmA_array)
if sBmA == 1:
    THDM2 = anyBSM('THDMII', scheme_name = 'OSalignment')
else:
    THDM2 = anyBSM('THDMII', scheme_name = 'OS')
THDM2.setparameters({'MAh2': mA, 'MHm2': mHpm, 'Mh2': mH, 'M': M, 'SinBmA': sBmA, 'TanBeta': np.tan(Beta)}) #Define new mass in anyBSM

β=sBmA_array
α=-BmA+β
l5=fcs.lamb5(mA, fcs.m122M(M, np.sin(β), np.cos(β)), np.sin(β), np.cos(β), v, 0, 0)

C_xs = np.array([qtcp.C_93(α, β, mH, l5),qtcp.C_94(α, β, mH, mA, l5),qtcp.C_102(α, β, mH, mA, l5),qtcp.C_123(α, β, mH, l5),qtcp.C_140(α, β, mH, l5)])
cnd = spr.perturbative_unitarity_const(C_xs)

THDM2.progress=False
THDM2.warnSSSS=False
a0=[THDM2.eigSSSS(parameters={'TanBeta':x}) for x in sBmA_array]

for i,x in enumerate(sBmA_array):
    THDM2.setparameters({'TanBeta': x})
    if (a0[i]>0.5):
        Γ[i] = np.nan
        Γ0[i] = np.nan
    else:
        lamb = THDM2.lambdahhh()
        Γ[i] = np.real(lamb["total"])
        Γ0[i] = np.real(lamb["treelevel"])
    # lamb = THDM2.lambdahhh()
    # Γ[i] = np.real(lamb["total"])
    # Γ0[i] = np.real(lamb["treelevel"])

#%%

y1=-Gammahhh_treelevel_cos(np.arcsin(sBmA),np.arctan(sBmA_array), M)/lambdahhhSM
y2=-Gammahhh_oneloop_cos_kan(np.arcsin(sBmA),np.arctan(sBmA_array), M,mH,mA,mHpm)/lambdahhhSM
y3=-Gammahhh_oneloop_cos(np.arcsin(sBmA),np.arctan(sBmA_array), M,mH,mA,mHpm)/lambdahhhSM
y4=-Gammahhh_oneloop_cos_correc(np.arcsin(sBmA),np.arctan(sBmA_array), M,mH,mA,mHpm)/lambdahhhSM
y5=-Gammahhh_oneloop_cos_correc_coup_full(np.arcsin(sBmA),np.arctan(sBmA_array), M,mH,mA,mHpm)/lambdahhhSM

y1,y2,y3,y4,y5=[np.where(np.array(a0)<0.5,y,np.nan) for y in [y1,y2,y3,y4,y5]]

fig, ax = plt.subplots(figsize=(10, 10),dpi=resol)

plt.vlines(sBmA_array[np.where(np.isnan(Γ))[0][0]-1],0.35,2.5,colors='k',alpha=0.2)
plt.plot(sBmA_array, Γ0/lambdahhhSM,c='b',label='Full calculation (tree-level)')
plt.plot(sBmA_array, Γ/lambdahhhSM,c='k',label='Full calculation')
plt.plot(sBmA_array, y1,c='y',label='Leading contribution (tree-level)',ls='dashdot')
plt.plot(sBmA_array, y2,c='green',ls='dashed',label='Leading contribution')
plt.plot(sBmA_array, y3,c='orange',label='Leading contribution (+Top)')
plt.plot(sBmA_array, y4,c='red',ls='dashed',label='Leading contribution (+Top+Higgs)')
plt.plot(sBmA_array, y5,c='magenta',ls='solid',label='Leading contribution coup(+Top+Higgs)')
plt.yticks(size=25)
plt.ylabel(r'$\kappa_\lambda$', size=25)
plt.xticks(size=25)
plt.xlabel(r'$\tan{\beta}$', size=25)
ax.grid()

textstr = '\n'.join((
    r'\textbf{2HDM-II}',
    r'$m_A=%d$ GeV' %(mA),
    r'$m_{H}=%d$ GeV' %(mH),
    r'$m_{H^\pm}=%d$ GeV' %(mHpm),
    r'$M=%d$ GeV' %(M),
    r'$\sin(\beta-\alpha)=%.2f$' %(sBmA),
#    r'$\tan(\beta)=%.1f$' %(np.tan(Beta)),
    ))

ax.yaxis.get_offset_text().set_fontsize(20)
ax.text(0.05, 0.91, textstr, transform=ax.transAxes, fontsize=20, verticalalignment='top')
ax.text(0.72, 0.55, r'\textbf{Perturbative unitarity}'+'\n'+r'\textbf{constraint}', transform=ax.transAxes, fontsize=18, verticalalignment='top',horizontalalignment='center',alpha=0.8)

plt.ylim(0.35,2.5)
plt.xlim(0,40)
ax.fill_between(sBmA_array[np.where(np.isnan(Γ))[0][0]-1:],0.35,2.5,hatch="/",color='k',alpha=0.2)
plt.minorticks_on()
plt.tick_params(axis='both', which='both', right=True, top=True, direction='in')

plt.legend(fontsize=18,loc='lower right',framealpha=0.99)

# Make the zoom-in plot:
    
x1 = 0
x2 = 5

# select y-range for zoomed region
y1 = 0.35
y2 = 0.7

axins = zoomed_inset_axes(ax, 2, loc='upper right', bbox_to_anchor=(0.35, 0.38, 0.5, 0.5), bbox_transform=ax.figure.transFigure) # zoom = 2
axins.plot(sBmA_array, Γ0/lambdahhhSM,c='b',label='AnyBSM (tree-level)')
axins.plot(sBmA_array, Γ/lambdahhhSM,c='k',label='AnyBSM')
#plt.plot(sBmA_array, -Gammahhh_treelevel(np.arcsin(sBmA_array)-np.pi/2, M)/lambdahhhSM,c='y',label='Kanemura (tree-level) Approx',ls='dashed')
#plt.plot(sBmA_array, -Gammahhh_oneloop(np.arcsin(sBmA_array)-np.pi/2, M,mH,mA,mHpm)/lambdahhhSM,c='orange',label='Kanemura (one-loop) Approx',ls='dashed')
axins.plot(sBmA_array, -Gammahhh_treelevel_cos(np.arcsin(sBmA),np.arctan(sBmA_array), M)/lambdahhhSM,c='y',label='Kanemura (tree-level)',ls='dashdot')
axins.plot(sBmA_array, -Gammahhh_oneloop_cos_kan(np.arcsin(sBmA),np.arctan(sBmA_array), M,mH,mA,mHpm)/lambdahhhSM,c='green',ls='dashed',label='Kanemura (one-loop)')
axins.plot(sBmA_array, -Gammahhh_oneloop_cos(np.arcsin(sBmA),np.arctan(sBmA_array), M,mH,mA,mHpm)/lambdahhhSM,c='orange',label='Kanemura (Top) (one-loop)')
axins.plot(sBmA_array, -Gammahhh_oneloop_cos_correc(np.arcsin(sBmA),np.arctan(sBmA_array), M,mH,mA,mHpm)/lambdahhhSM,c='red',ls='dashed',label='Kanemura (Top+Higgs). (one-loop)')
axins.set_xlim(x1, x2)
axins.set_ylim(y1, y2)
plt.minorticks_on()
plt.tick_params(axis='both', which='both', right=True, top=True, direction='in')
plt.xticks(visible=True,fontsize=20)
plt.yticks(visible=True,fontsize=20)
axins.grid()
mark_inset(ax, axins, loc1=3, loc2=2, fc="none", ec="0.7")
plt.draw()

plt.savefig('../temp/kappa_tanb.png', dpi=600,bbox_inches='tight')

plt.show()

# #%%%

# ###                         Plot 5 - 2HDM 1-loop comparison as function of tanb - Models

# #Standard definitions in anyBSM
# mA = 500 #Pseudoscalar mass
# mHpm = 400 #Charged scalar mass
# mH = 100  #Heavy Higgs mass
# sBmA = 0.9 #Standard definition of sin(β-α) in anyBSM
# M = 800 #Standard definition of M in anyBSM
# BmA = np.arcsin(sBmA)
# Beta = np.arctan(2)


# sBmA_array = np.linspace(0.8,40,40)
# Γ = np.copy(sBmA_array)
# Γ2 = np.copy(sBmA_array)
# Γ0 = np.copy(sBmA_array)
# THDM1 = anyBSM('THDMI', scheme_name = 'OSalignment')
# THDM2 = anyBSM('THDMII', scheme_name = 'OSalignment')
# THDM1.setparameters({'MAh2': mA, 'MHm2': mHpm, 'Mh2': mH, 'M': M, 'SinBmA': sBmA, 'TanBeta': np.tan(Beta)}) #Define new mass in anyBSM
# THDM2.setparameters({'MAh2': mA, 'MHm2': mHpm, 'Mh2': mH, 'M': M, 'SinBmA': sBmA, 'TanBeta': np.tan(Beta)}) #Define new mass in anyBSM

# for i,x in enumerate(sBmA_array):
#     THDM1.setparameters({'TanBeta': x})
#     lamb1 = THDM1.lambdahhh()
#     lamb2 = THDM2.lambdahhh()
#     Γ[i] = np.real(lamb1["total"])
#     Γ2[i] = np.real(lamb2["total"])
#     Γ0[i] = np.real(lamb1["treelevel"])

# fig, ax = plt.subplots(figsize=(10, 10))

# plt.plot(sBmA_array, Γ0/lambdahhhSM,c='y',label='Tree-level')
# plt.plot(sBmA_array, Γ/lambdahhhSM,c='k',label='Type I')
# plt.plot(sBmA_array, Γ/lambdahhhSM,c='r',label='Type II')


# plt.yticks(size=20)
# plt.ylabel(r'$\kappa_\lambda$', size=25)
# plt.xticks(size=20)
# plt.xlabel(r'$\tan{\beta}$', size=25)
# ax.grid()

# textstr = '\n'.join((
#     r'\textbf{2HDM-II}',
#     r'$M_A=%d$' %(mA),
#     r'$M_{H}=%d$' %(mH),
#     r'$M_{H^\pm}=%d$' %(mHpm),
#     r'$M=%d$' %(M),
#     r'$\sin(\beta-\alpha)=%.2f$' %(sBmA),
# #    r'$\tan(\beta)=%.1f$' %(np.tan(Beta)),
#     ))

# ax.yaxis.get_offset_text().set_fontsize(20)
# ax.text(0.05, 0.25, textstr, transform=ax.transAxes, fontsize=18, verticalalignment='top')

# plt.legend(fontsize=20,loc='lower right')

# plt.show()

# #%%                         Plots 2HDM

# #Standard definitions in anyBSM
# mA = 300 #Pseudoscalar mass
# mHpm = 300 #Charged scalar mass
# mH = 600  #Heavy Higgs mass
# sBmA = 1 #Standard definition of sin(β-α) in anyBSM
# M = 600 #Standard definition of M in anyBSM
# BmA = np.arcsin(sBmA)
# Beta = np.arctan(2)

# ###                         Plot 3 - 2HDM 1-loop comparison as function of M_A

# M_array = np.linspace(0,500,100)
# Γ = np.copy(M_array)
# Γ0 = np.copy(M_array)
# if sBmA == 1:
#     scheme_str = 'OSalignment'
# else:
#     scheme_str = 'OS'
    

# THDM2 = anyBSM('THDMI', scheme_name = scheme_str)
# THDM2.setparameters({'MAh2': mA, 'MHm2': mHpm, 'Mh2': mH, 'M': M, 'SinBmA': sBmA, 'TanBeta': np.tan(Beta)}) #Define new mass in anyBSM

# for i,x in enumerate(M_array):
#     THDM2.setparameters({'M': x})
#     lamb = THDM2.lambdahhh()
#     Γ[i] = np.real(lamb["total"])
#     Γ0[i] = np.real(lamb["treelevel"])

# fig, ax = plt.subplots(figsize=(10, 10))

# plt.plot(M_array, Γ/lambdahhhSM,c='k',label='AnyBSM')
# plt.plot(M_array, Γ0/lambdahhhSM,c='b',label='AnyBSM (tree-level)')
# plt.plot(M_array, -Gammahhh_treelevel(np.arcsin(sBmA)-np.pi/2, M_array)/lambdahhhSM,c='y',label='Kanemura (tree-level) (Approx. x)',ls='dashed')
# plt.plot(M_array, -Gammahhh_oneloop(np.arcsin(sBmA)-np.pi/2, M_array,mH,mA,mHpm)/lambdahhhSM,c='orange',label='Kanemura (one-loop) (Approx. x)',ls='dashed')
# plt.plot(M_array, -Gammahhh_treelevel_cos(BmA,Beta, M_array)/lambdahhhSM,c='y',label='Kanemura (tree-level)',ls='dashdot')
# plt.plot(M_array, -Gammahhh_oneloop_cos_kan(BmA,Beta, M_array,mH,mA,mHpm)/lambdahhhSM,c='green',label='Kanemura (one-loop)')
# plt.plot(M_array, -Gammahhh_oneloop_cos(BmA,Beta, M_array,mH,mA,mHpm)/lambdahhhSM,c='orange',ls='dashdot',label='Kanemura (Top) (one-loop)')
# plt.plot(M_array, -Gammahhh_oneloop_cos_correc(BmA,Beta, M_array,mH,mA,mHpm)/lambdahhhSM,c='red',ls='dashed',label='Kanemura (Top+Higgs). (one-loop)')
# plt.yticks(size=20)
# plt.ylabel(r'$\kappa_\lambda$', size=25)
# plt.xticks(size=20)
# plt.xlabel(r'$M$ [GeV]', size=25)
# ax.grid()

# textstr = '\n'.join((
#     r'\textbf{2HDM-II}',
#     r'$M_A=%d$' %(mA),
#     r'$M_{H}=%d$' %(mH),
#     r'$M_{H^\pm}=%d$' %(mHpm),
# #    r'$M=%d$' %(M),
#     r'$\sin(\beta-\alpha)=%.2f$' %(sBmA),
#     r'$\tan(\beta)=%.1f$' %(np.tan(Beta)),
#     ))

# ax.yaxis.get_offset_text().set_fontsize(20)
# ax.text(0.05, 0.75, textstr, transform=ax.transAxes, fontsize=18, verticalalignment='top')


# plt.legend(fontsize=20)

# plt.show()

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
    
#%%

fig, ax = plt.subplots(figsize=(10, 10),dpi=resol)

plt.plot(M_array, Γ/lambdahhhSM,c='k',label='Full calculation')
#plt.plot(M_array,-Gammahhh_treelevel_cos(BmA, Beta, M)/lambdahhhSM+(M_array**4/(12*np.pi**2*mh**2*v**2)))#-3*M**2*M_array**2/(12*np.pi**2*mh**2*v**2))#*(1-M**2/(M_array**2))**3)
#plt.plot(M_array, Γ0/lambdahhhSM,c='b',label='AnyBSM (tree-level)')
#plt.plot(M_array, -Gammahhh_treelevel(np.arcsin(sBmA)-np.pi/2, M)/lambdahhhSM,c='y',label='Kanemura (tree-level) (Approx. x)',ls='dashed')
#plt.plot(M_array, -Gammahhh_oneloop(np.arcsin(sBmA)-np.pi/2, M,mH,M_array,mHpm)/lambdahhhSM,c='orange',label='Kanemura (one-loop) (Approx. x)',ls='dashed')
#plt.plot(M_array, -Gammahhh_treelevel_cos(BmA,Beta, M)/lambdahhhSM,c='y',label='Kanemura (tree-level)',ls='dashdot')
plt.plot(M_array, -Gammahhh_oneloop_cos_kan(BmA,Beta, M,mH,M_array,mHpm)/lambdahhhSM,c='green',label='Leading contribution')
plt.plot(M_array, -Gammahhh_oneloop_cos(BmA,Beta, M,mH,M_array,mHpm)/lambdahhhSM,c='orange',ls='dashdot',label='Leading contribution (+Top)')
plt.plot(M_array, -Gammahhh_oneloop_cos_correc(BmA,Beta, M,mH,M_array,mHpm)/lambdahhhSM,c='red',ls='dashed',label='Leading contribution (+Top+Higgs)')
plt.yticks(size=25)
#plt.yscale('log')
#plt.xscale('log')
plt.ylabel(r'$\kappa_\lambda$', size=25)
plt.xticks(size=25)
plt.xlabel(r'$M_A$ [GeV]', size=25)
ax.grid()
plt.minorticks_on()
plt.tick_params(axis='both', which='both', right=True, top=True, direction='in')

#plt.ylim(-10,10)

textstr = '\n'.join((
    r'\textbf{2HDM-II}',
#    r'$M_A=%d$' %(mA),
    r'$m_{H}=%d$ GeV' %(mH),
    r'$m_{H^\pm}=%d$ GeV' %(mHpm),
    r'$M=%d$ GeV' %(M),
    r'$\sin(\beta-\alpha)=%.2f$' %(sBmA),
    r'$\tan(\beta)=%.1f$' %(np.tan(Beta)),
    ))

ax.yaxis.get_offset_text().set_fontsize(20)
ax.text(0.05, 0.92, textstr, transform=ax.transAxes, fontsize=20, verticalalignment='top')


plt.legend(fontsize=18,loc='center left')

plt.show()

#%%                         Plots 2HDM - 331

#Standard definitions in anyBSM
mA = 300 #Pseudoscalar mass
mHpm = 300 #Charged scalar mass
mH = 300  #Heavy Higgs mass
sBmA = 1 #Standard definition of sin(β-α) in anyBSM
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
    

THDM2 = anyBSM('THDMII', scheme_name = scheme_str)
THDM2.setparameters({'MAh2': mA, 'MHm2': mHpm, 'Mh2': mH, 'M': M, 'SinBmA': sBmA, 'TanBeta': np.tan(Beta)}) #Define new mass in anyBSM

for i,x in enumerate(M_array):
    THDM2.setparameters({'MAh2': x, 'M': x})
    lamb = THDM2.lambdahhh()
    Γ[i] = np.real(lamb["total"])
    Γ0[i] = np.real(lamb["treelevel"])
    
#%%

fig, ax = plt.subplots(figsize=(10, 10),dpi=resol)

plt.plot(M_array, Γ/lambdahhhSM,c='k',label='Full calculation')
#plt.plot(M_array, Γ0/lambdahhhSM,c='b',label='AnyBSM (tree-level)')
#plt.plot(M_array, -Gammahhh_treelevel(np.arcsin(sBmA)-np.pi/2, M)/lambdahhhSM,c='y',label='Kanemura (tree-level) (Approx. x)',ls='dashed')
#plt.plot(M_array, -Gammahhh_oneloop(np.arcsin(sBmA)-np.pi/2, M_array,mH,M_array,mHpm)/lambdahhhSM,c='orange',label='Kanemura (one-loop) (Approx. x)',ls='dashed')
#plt.plot(M_array, -Gammahhh_treelevel_cos(BmA,Beta, M)/lambdahhhSM,c='y',label='Kanemura (tree-level)',ls='dashdot')
plt.plot(M_array, -Gammahhh_oneloop_cos_kan(BmA,Beta, M_array,mH,M_array,mHpm)/lambdahhhSM,c='green',label='Leading contribution')
plt.plot(M_array, -Gammahhh_oneloop_cos(BmA,Beta, M_array,mH,M_array,mHpm)/lambdahhhSM,c='orange',ls='dashdot',label='Leading contribution (+Top)')
plt.plot(M_array, -Gammahhh_oneloop_cos_correc(BmA,Beta, M_array,mH,M_array,mHpm)/lambdahhhSM,c='red',ls='dashed',label='Leading contribution (+Top+Higgs)')
plt.yticks(size=25)
plt.ylabel(r'$\kappa_\lambda$', size=25)
plt.xticks(size=25)
plt.xlabel(r'$M_A$ [GeV]', size=25)
ax.grid()
plt.minorticks_on()
plt.tick_params(axis='both', which='both', right=True, top=True, direction='in')

#plt.ylim(-10,10)

textstr = '\n'.join((
    r'\textbf{2HDM-331 EFT}',
#    r'$M_A=%d$' %(mA),
    r'$m_{H}=%d$ GeV' %(mH),
    r'$m_{H^\pm}=%d$ GeV' %(mHpm),
#    r'$M=%d$' %(M),
    r'$\sin(\beta-\alpha)=%.2f$' %(sBmA),
    r'$\tan(\beta)=%.1f$' %(np.tan(Beta)),
    ))

ax.yaxis.get_offset_text().set_fontsize(20)
ax.text(0.05, 0.6, textstr, transform=ax.transAxes, fontsize=20, verticalalignment='top')


plt.legend(fontsize=18)

plt.show()


# #%%                         Plots 2HDM - Weiglein

# #Standard definitions in anyBSM
# mA = 300 #Pseudoscalar mass
# mHpm = 300 #Charged scalar mass
# mH = 600  #Heavy Higgs mass
# sBmA = 1 #Standard definition of sin(β-α) in anyBSM
# M = 600 #Standard definition of M in anyBSM
# BmA = np.arcsin(sBmA)
# Beta = np.arctan(2)

# ###                         Plot 3 - 2HDM 1-loop comparison as function of M_A

# M_array = np.linspace(600,1050,100)
# Γ = np.copy(M_array)
# Γ0 = np.copy(M_array)
# if sBmA == 1:
#     scheme_str = 'OSalignment'
# else:
#     scheme_str = 'OS'
    

# THDM2 = anyBSM('THDMI', scheme_name = scheme_str)
# THDM2.setparameters({'MAh2': mA, 'MHm2': mHpm, 'Mh2': mH, 'M': M, 'SinBmA': sBmA, 'TanBeta': np.tan(Beta)}) #Define new mass in anyBSM

# for i,x in enumerate(M_array):
#     THDM2.setparameters({'MAh2': x, 'MHm2': x})
#     lamb = THDM2.lambdahhh()
#     Γ[i] = np.real(lamb["total"])
#     Γ0[i] = np.real(lamb["treelevel"])

# fig, ax = plt.subplots(figsize=(10, 10))

# plt.plot(M_array, Γ/lambdahhhSM,c='k',label='AnyBSM')
# #plt.plot(M_array, Γ0/lambdahhhSM,c='b',label='AnyBSM (tree-level)')
# #plt.plot(M_array, -Gammahhh_treelevel(np.arcsin(sBmA)-np.pi/2, M)/lambdahhhSM,c='y',label='Kanemura (tree-level) (Approx. x)',ls='dashed')
# #plt.plot(M_array, -Gammahhh_oneloop(np.arcsin(sBmA)-np.pi/2, M_array,mH,M_array,mHpm)/lambdahhhSM,c='orange',label='Kanemura (one-loop) (Approx. x)',ls='dashed')
# #plt.plot(M_array, -Gammahhh_treelevel_cos(BmA,Beta, M)/lambdahhhSM,c='y',label='Kanemura (tree-level)',ls='dashdot')
# #plt.plot(M_array, -Gammahhh_oneloop_cos_kan(BmA,Beta, M_array,mH,M_array,mHpm)/lambdahhhSM,c='green',label='Kanemura (one-loop)')
# #plt.plot(M_array, -Gammahhh_oneloop_cos(BmA,Beta, M_array,mH,M_array,mHpm)/lambdahhhSM,c='orange',ls='dashdot',label='Kanemura (Top) (one-loop)')
# #plt.plot(M_array, -Gammahhh_oneloop_cos_correc(BmA,Beta, M_array,mH,M_array,mHpm)/lambdahhhSM,c='red',ls='dashed',label='Kanemura (Top+Higgs). (one-loop)')
# plt.yticks(size=20)
# plt.ylabel(r'$\kappa_\lambda$', size=25)
# plt.xticks(size=20)
# plt.xlabel(r'$M_A$ [GeV]', size=25)
# ax.grid()

# #plt.ylim(-10,10)

# textstr = '\n'.join((
#     r'\textbf{2HDM-I}',
# #    r'$M_A=%d$' %(mA),
#     r'$M_{H}=%d$' %(mH),
# #    r'$M_{H^\pm}=%d$' %(mHpm),
#     r'$M=%d$' %(M),
#     r'$\sin(\beta-\alpha)=%.2f$' %(sBmA),
#     r'$\tan(\beta)=%.1f$' %(np.tan(Beta)),
#     ))

# ax.yaxis.get_offset_text().set_fontsize(20)
# ax.text(0.05, 0.45, textstr, transform=ax.transAxes, fontsize=18, verticalalignment='top')


# plt.legend(fontsize=20)

# plt.show()
