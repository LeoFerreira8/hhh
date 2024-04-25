#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Apr 22 11:09:52 2024

@author: leonardoferreira
"""

import pandas as pd
import numpy as np

alignment = False
small_l5 = True
l5_size_max = 1e-10
non_alignment_max = 0.2

### ---------      SM param
MW = 80.379
MZ = 91.187
mtSM = 172.5
mt = mtSM
alphaQED = 137.035999679

e = 2*np.sqrt(np.pi/(alphaQED))
vSM = 2*MW*np.sqrt(1-MW**2/MZ**2)/e
v = vSM
mhSM = 125.1 #GeV
mh = mhSM
### ----------

### ----------      intervals
mA_min = 125
mA_max = 1000

mH_min = 125
mH_max = 1000

mHpm_min = 125
mHpm_max = 1000

cosa_min = -1
cosa_max = 1

tanb_min = 1.5
tanb_max = 5

M_min = 1e2
M_max = 1e3
### ----------

def Gammahhh_treelevel(x,M):
    G = -(3*mh**2/v)*(1+(3/2)*(1-4*M**2/(3*mh**2))*x**2)
    
    return G

def Gammahhh_oneloop(x,M,mH,mA,mHpm):
    G = -(3*mh**2/v)*(1+(3/2)*(1-4*M**2/(3*mh**2))*x**2-3*mt**4/(3*np.pi**2*mh**2*v**2)+(mH**4/(12*np.pi**2*mh**2*v**2))*(1-M**2/(mH**2))**3+(mA**4/(12*np.pi**2*mh**2*v**2))*(1-M**2/(mA**2))**3+(mHpm**4/(6*np.pi**2*mh**2*v**2))*(1-M**2/(mHpm**2))**3)
    
    return G

def Gammahhh_treelevel_cos(ba,b,M):
    G = -(3*mh**2/(2*v*np.sin(2*b)))*(np.cos(-3*ba+2*b)+3*np.cos(-ba+2*b)-4*np.cos(ba)**2*np.cos(-ba+2*b)*(M**2/mh**2))
    
    return G

def Gammahhh_oneloop_cos(ba,b,M,mH,mA,mHpm):
    G = Gammahhh_treelevel_cos(ba, b, M)-(3*mh**2/v)*(-(np.cos(ba)/np.tan(b)+np.sin(ba))**3*3*mt**4/(3*np.pi**2*mh**2*v**2)+(mH**4/(12*np.pi**2*mh**2*v**2))*(1-M**2/(mH**2))**3+(mA**4/(12*np.pi**2*mh**2*v**2))*(1-M**2/(mA**2))**3+(mHpm**4/(6*np.pi**2*mh**2*v**2))*(1-M**2/(mHpm**2))**3)
    #G = Gammahhh_treelevel_cos(ba, b, M)-(3*mh**2/v)*(-3*mt**4/(3*np.pi**2*mh**2*v**2)+(mH**4/(12*np.pi**2*mh**2*v**2))*(1-M**2/(mH**2))**3+(mA**4/(12*np.pi**2*mh**2*v**2))*(1-M**2/(mA**2))**3+(mHpm**4/(6*np.pi**2*mh**2*v**2))*(1-M**2/(mHpm**2))**3)
    
    return G

def find_random_points(N):
    '''Generate N random points in the parameter space. mA, mH, mHpm varying from 125 to 1000, cos(alpha) from -1 to 1, tan(beta) from 1.5 to 5 and M from 1e3 to 1e7.
        Returns: mA, mH, mHpm, cosa, tanb, M.
    '''    
    
    if not alignment:    
        DataFrame = pd.DataFrame(columns = ['mA', 'mH', 'mHpm', 'cosa', 'tanb', 'M'])
        
        numb = np.random.rand(6,N)
        
        DataFrame['mA'] = mA_min+numb[0]*(mA_max-mA_min)
        DataFrame['mH'] = mH_min+numb[1]*(mH_max-mH_min)
        DataFrame['mHpm'] = mHpm_min+numb[2]*(mHpm_max-mHpm_min)
        DataFrame['tanb'] = tanb_min+numb[4]*(tanb_max-tanb_min)
        DataFrame['cosa'] = np.cos(beta(DataFrame['tanb'])-(np.pi/2)*(1+np.random.choice([-1,1],size=N)*numb[3]*non_alignment_max))
        #DataFrame['M'] = 10**(3+numb[5]*(7-3))
        if small_l5:
            DataFrame['M'] = np.sqrt(m122(DataFrame['mA'], np.sin(beta(DataFrame['tanb'])), np.cos(beta(DataFrame['tanb'])),vSM,numb[5]*l5_size_max,0,0)/(np.sin(beta(DataFrame['tanb']))*np.cos(beta(DataFrame['tanb']))))
        else:
            DataFrame['M'] = M_min+numb[5]*(M_max-M_min)
    else:
        DataFrame = pd.DataFrame(columns = ['mA', 'mH', 'mHpm', 'cosa', 'tanb', 'M'])
        
        numb = np.random.rand(6,N)
        
        DataFrame['mA'] = mA_min+numb[0]*(mA_max-mA_min)
        DataFrame['mH'] = mH_min+numb[1]*(mH_max-mH_min)
        DataFrame['mHpm'] = mHpm_min+numb[2]*(mHpm_max-mHpm_min)
        DataFrame['tanb'] = tanb_min+numb[4]*(tanb_max-tanb_min)
        DataFrame['cosa'] = np.cos(beta(DataFrame['tanb'])-np.pi/2)
        #DataFrame['M'] = 10**(3+numb[5]*(7-3))
        if small_l5:
            DataFrame['M'] = np.sqrt(m122(DataFrame['mA'], np.sin(beta(DataFrame['tanb'])), np.cos(beta(DataFrame['tanb'])),vSM,numb[5]*l5_size_max,0,0)/(np.sin(beta(DataFrame['tanb']))*np.cos(beta(DataFrame['tanb']))))
        else:
            DataFrame['M'] = M_min+numb[5]*(M_max-M_min)
        
    return DataFrame
    
def lamb1(mH,mh,m122,ca,sa,sb,cb,v,lambda6, lambda7):
    return (mH**2*ca**2+mh**2*sa**2-m122*sb/cb)/(v**2*cb**2)-(3/2)*(lambda6*sb/cb)+(1/2)*lambda7*(sb/cb)**3
    
def lamb2(mH,mh,m122,ca,sa,sb,cb,v,lambda6, lambda7):
    return (mH**2*sa**2+mh**2*ca**2-m122*cb/sb)/(v**2*sb**2)-(3/2)*(lambda7*cb/sb)+(1/2)*lambda6*(cb/sb)**3
    
def lamb3(mH,mh,mHpm,m122,ca,sa,sb,cb,v,lambda6, lambda7):
    return ((mH**2-mh**2)*sa*ca+2*mHpm**2*sb*cb-m122)/(v**2*sb*cb)-(1/2)*(lambda6*cb/sb)-(1/2)*lambda7*(sb/cb)
    
def lamb4(mA,mHpm,m122,sb,cb,v,lambda6, lambda7):
    return ((mA**2-2*mHpm**2)*sb*cb+m122)/(v**2*sb*cb)-(1/2)*(lambda6*cb/sb)-(1/2)*lambda7*(sb/cb)

def lamb5(mA,m122,sb,cb,v,lambda6, lambda7):
    return ((m122-mA**2*sb*cb))/(v**2*sb*cb)-(1/2)*(lambda6*cb/sb)-(1/2)*lambda7*(sb/cb)
    
def m122(mA,sb,cb,v,lambda5,lambda6, lambda7):
    return mA**2*sb*cb+(lambda5+(1/2)*lambda6*cb/sb+(1/2)*lambda7*sb/cb)*(v**2*sb*cb)

def m122M(M,sb,cb):
    return M**2*sb*cb
    
def sinbma(alpha,beta):
    return np.sin(alpha-beta)
    
def beta(tanb):
    return np.arctan(tanb)
    
def alpha(cosa):
    return -np.arccos(cosa)

