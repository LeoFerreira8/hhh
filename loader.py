"""
Created on Sun May 26 00:13:53 2024

@author: leo
"""

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import scan_parameterspace_funcs as fcs
import scan_parameterspace as spr
import quartic_couplings as qtcp

#%%%                                Model parameters to load
THDM_type = ''
small_l5 = False
alignment = True
non_alignment_max = 0.2
load_CDF = True

######## loading

if THDM_type=='':
    latex_model = '2HDM - Type I'
    THDM_type=''
elif THDM_type=='II' and not small_l5:
    latex_model = '2HDM - Type II'
elif THDM_type=='II' and small_l5:
    latex_model = '2HDM - 331 EFT'
    
if alignment:
    latex_alignment = r'($\beta-\alpha=\frac{\pi}{2}$)'
else:
    latex_alignment = r'($\beta-\alpha\in[\frac{\pi}{2}\pm%.2f]$)' %(non_alignment_max,)

if alignment:
    strga = 'A'
else:
    strg = str(non_alignment_max)
    strga = 'NA'+strg
    
if small_l5:
    strgl5 = '331lk'
else:
    strgl5 = ''

set_dir = 'data_'+'THDM'+THDM_type+strgl5+'-'+strga+'/'

TableTot = pd.read_csv('./'+set_dir+'/THDM'+THDM_type+strgl5+'-'+strga+'-Theo_PDG.csv')
TableTot_STU = pd.read_csv('./'+set_dir+'/THDM'+THDM_type+strgl5+'-'+strga+'-STU_PDG.csv')
TableTot_STU_Collid = pd.read_csv('./'+set_dir+'/THDM'+THDM_type+strgl5+'-'+strga+'-Collid_PDG.csv')
TableTot_STU_Collid_BSG = pd.read_csv('./'+set_dir+'/THDM'+THDM_type+strgl5+'-'+strga+'-BSG_PDG.csv')
TableTot_STU_Collid_BSG_unit = pd.read_csv('./'+set_dir+'/THDM'+THDM_type+strgl5+'-'+strga+'-PU_PDG.csv')

Dataset_teo = {'name': r'T', 'data': TableTot}
Dataset_stu = {'name': r'T+EW', 'data': TableTot_STU}
Dataset_col = {'name': r'T+EW+C', 'data': TableTot_STU_Collid}
Dataset_bsg = {'name': r'T+EW+C+bs$\gamma$', 'data': TableTot_STU_Collid_BSG}
Dataset_unit = {'name': r'T+EW+C+bs$\gamma$+PU', 'data': TableTot_STU_Collid_BSG_unit}

####### loading CDF data

if load_CDF:

    set_dir = 'data_'+'THDM'+THDM_type+strgl5+'-'+strga+'_CDF_'+'/'
    
    TableTot_CDF = pd.read_csv('./'+set_dir+'/THDM'+THDM_type+strgl5+'-'+strga+'-Theo.csv')
    TableTot_STU_CDF = pd.read_csv('./'+set_dir+'/THDM'+THDM_type+strgl5+'-'+strga+'-STU.csv')
    TableTot_STU_Collid_CDF = pd.read_csv('./'+set_dir+'/THDM'+THDM_type+strgl5+'-'+strga+'-Collid.csv')
    TableTot_STU_Collid_BSG_CDF = pd.read_csv('./'+set_dir+'/THDM'+THDM_type+strgl5+'-'+strga+'-BSG.csv')
    TableTot_STU_Collid_BSG_unit_CDF = pd.read_csv('./'+set_dir+'/THDM'+THDM_type+strgl5+'-'+strga+'-PU.csv')
    
    Dataset_teo_cdf = {'name': r'T', 'data': TableTot_CDF}
    Dataset_stu_cdf = {'name': r'T+EW (CDF)', 'data': TableTot_STU_CDF}
    Dataset_col_cdf = {'name': r'T+EW (CDF)+C', 'data': TableTot_STU_Collid_CDF}
    Dataset_bsg_cdf = {'name': r'T+EW (CDF)+C+bs$\gamma$', 'data': TableTot_STU_Collid_BSG_CDF}
    Dataset_unit_cdf = {'name': r'T+EW (CDF)+C+bs$\gamma$+PU', 'data': TableTot_STU_Collid_BSG_unit_CDF}
    
    set_dir = 'data_'+'THDM'+THDM_type+strgl5+'-'+strga+'_Comb_'+'/'
    
    TableTot_Comb = pd.read_csv('./'+set_dir+'/THDM'+THDM_type+strgl5+'-'+strga+'-Theo.csv')
    TableTot_STU_Comb = pd.read_csv('./'+set_dir+'/THDM'+THDM_type+strgl5+'-'+strga+'-STU.csv')
    TableTot_STU_Collid_Comb = pd.read_csv('./'+set_dir+'/THDM'+THDM_type+strgl5+'-'+strga+'-Collid.csv')
    TableTot_STU_Collid_BSG_Comb = pd.read_csv('./'+set_dir+'/THDM'+THDM_type+strgl5+'-'+strga+'-BSG.csv')
    TableTot_STU_Collid_BSG_unit_Comb = pd.read_csv('./'+set_dir+'/THDM'+THDM_type+strgl5+'-'+strga+'-PU.csv')
    
    Dataset_teo_comb = {'name': r'T', 'data': TableTot_Comb}
    Dataset_stu_comb = {'name': r'T+EW (Comb)', 'data': TableTot_STU_Comb}
    Dataset_col_comb = {'name': r'T+EW (Comb)+C', 'data': TableTot_STU_Collid_Comb}
    Dataset_bsg_comb = {'name': r'T+EW (Comb)+C+bs$\gamma$', 'data': TableTot_STU_Collid_BSG_Comb}
    Dataset_unit_comb = {'name': r'T+EW (Comb)+C+bs$\gamma$+PU', 'data': TableTot_STU_Collid_BSG_unit_Comb}

#%%                                 Plots

plt.rc('text', usetex=True)
plt.rc('xtick',labelsize=40)
plt.rc('ytick',labelsize=40)

def str_to_tex(strg):
    latex_parameters = [r'$m_A$ [GeV]',r'$m_H$ [GeV]',r'$m_{H^\pm}$ [GeV]',r'$\cos{\alpha}$',r'$\tan{\beta}$',r'$M$ [GeV]',r'$m_{12}^2$ [GeV$^2$]', r'$\lambda_1$', r'$\lambda_2$', r'$\lambda_3$', r'$\lambda_4$', r'$\lambda_5$',r'$\kappa_\lambda$',r'$\kappa_\lambda^{(0)}$',r'$\tilde{\kappa}_\lambda$',r'$\bar{\kappa}_\lambda$',r'$\sin{(\beta-\alpha)}$',r'$C_{93}$',r'$C_{94}$',r'$C_{102}$',r'$C_{123}$',r'$C_{140}$']
    
    if strg=='mA':
        return latex_parameters[0]
    elif strg=='mH':
        return latex_parameters[1]
    elif strg=='mHpm':
        return latex_parameters[2]
    elif strg=='cosa':
        return latex_parameters[3]
    elif strg=='tanb':
        return latex_parameters[4]
    elif strg=='M':
        return latex_parameters[5]
    elif strg=='m122':
        return latex_parameters[6]
    elif strg=='l1':
        return latex_parameters[7]
    elif strg=='l2':
        return latex_parameters[8]
    elif strg=='l3':
        return latex_parameters[9]
    elif strg=='l4':
        return latex_parameters[10]
    elif strg=='l5':
        return latex_parameters[11]
    elif strg=='kappa':
        return latex_parameters[12]
    elif strg=='kappa-tree':
        return latex_parameters[13]
    elif strg=='kappa-kan-x':
        return latex_parameters[14]
    elif strg=='kappa-kan':
        return latex_parameters[15]
    elif strg=='sino':
        return latex_parameters[16]
    elif strg=='c93':
        return latex_parameters[17]
    elif strg=='c94':
        return latex_parameters[18]
    elif strg=='c102':
        return latex_parameters[19]
    elif strg=='c123':
        return latex_parameters[20]
    elif strg=='c140':
        return latex_parameters[21]
    else:
        raise ValueError("Invalid parameter.")
    
def plotter(param1,param2,*dataset):

    fig, ax = plt.subplots(figsize=(10, 10))
    
    plt.title(latex_model +' '+ latex_alignment,size=20)
    
    for tbl in dataset:
        if param1=='sino':
            sino = np.sin(fcs.beta(tbl['data']['tanb'])-fcs.alpha(tbl['data']['cosa']))
            plt.scatter(sino,tbl['data'][param2],label=tbl['name'])
        elif param2=='sino':
            sino = np.sin(fcs.beta(tbl['data']['tanb'])-fcs.alpha(tbl['data']['cosa']))
            plt.scatter(tbl['data'][param1],sino,label=tbl['name'])
        else:
            plt.scatter(tbl['data'][param1],tbl['data'][param2],label=tbl['name'])
    plt.xlabel(str_to_tex(param1), size=25)
    plt.xticks(size=20)
    #plt.xlim(125,900)
    plt.ylabel(str_to_tex(param2), size=25)
    plt.yticks(size=20)
    #plt.ylim(125,900)
    ax.grid()
    
    plt.legend(fontsize=20)
    
    plt.show()

    return 0

def plotter_3(param1,param2,param3,*dataset):

    fig, ax = plt.subplots(figsize=(10, 10))
    
    plt.title(latex_model +' '+ latex_alignment,size=20)
    
    for tbl in dataset:
        sino = np.sin(fcs.beta(tbl['data']['tanb'])-fcs.alpha(tbl['data']['cosa']))
        proxer=tbl['data'].copy()
        proxer.insert(0,"sino", sino)
        proxer = proxer.sort_values(by=[param3],ascending=True)
        plt.scatter(proxer[param1],proxer[param2],c=proxer[param3])
        
    cbar=plt.colorbar(label=str_to_tex(param3))
    plt.xlabel(str_to_tex(param1), size=25)
    plt.xticks(size=20)
    #plt.xlim(125,900)
    plt.ylabel(str_to_tex(param2), size=25)
    plt.yticks(size=20)
    #plt.ylim(125,900)
    ax.grid()
    cbar.ax.tick_params(labelsize=20)  # Change size of tick labels
    cbar.ax.yaxis.label.set_size(20)   # Change size of colorbar label
    
    plt.show()

    return 0

def plotter_comp(param1,param2,param3,param4,*dataset):    

    fig, ax = plt.subplots(figsize=(10, 10))
    
    plt.title(latex_model +' '+ latex_alignment,size=20)
    
    for tbl in dataset:
        ratio = tbl['data'][param3]/tbl['data'][param4]
        if param1=='sino':
            sino = np.sin(fcs.beta(tbl['data']['tanb'])-fcs.alpha(tbl['data']['cosa']))
            plt.scatter(sino,tbl['data'][param2],c=ratio)
        elif param2=='sino':
            sino = np.sin(fcs.beta(tbl['data']['tanb'])-fcs.alpha(tbl['data']['cosa']))
            plt.scatter(tbl['data'][param1],sino,c=ratio)
        else:
            plt.scatter(tbl['data'][param1],tbl['data'][param2],c=ratio)
    cbar=plt.colorbar(label=str_to_tex(param3)+'/'+str_to_tex(param4))
    plt.xlabel(str_to_tex(param1), size=25)
    plt.xticks(size=20)
    #plt.xlim(125,900)
    plt.ylabel(str_to_tex(param2), size=25)
    plt.yticks(size=20)
    #plt.ylim(125,900)
    ax.grid()
    cbar.ax.tick_params(labelsize=20)  # Change size of tick labels
    cbar.ax.yaxis.label.set_size(20)   # Change size of colorbar label
    
    plt.show()

    return 0

def plotter_diff(param1,param2,param3,param4,*dataset):

    fig, ax = plt.subplots(figsize=(10, 10))
    
    plt.title(latex_model +' '+ latex_alignment,size=20)
    
    for tbl in dataset:
        plt.scatter(tbl['data'][param1]-tbl['data'][param2],tbl['data'][param3]-tbl['data'][param4],label=tbl['name'])
        
    plt.xlabel(str_to_tex(param1).replace(' [GeV]','')+'-'+str_to_tex(param2), size=25)
    plt.xticks(size=20)
    #plt.xlim(-450,200)
    plt.ylabel(str_to_tex(param3).replace(' [GeV]','')+'-'+str_to_tex(param4), size=25)
    plt.yticks(size=20)
    #plt.ylim(-350,350)
    ax.grid()
    
    plt.legend(loc='upper left',fontsize=16)
    
    plt.show()

    return 0

def plotter_diff_color(param1,param2,param3,param4,param5,*dataset):

    fig, ax = plt.subplots(figsize=(10, 10))
    
    plt.title(latex_model +' '+ latex_alignment,size=20)
    
    for tbl in dataset:
        sino = np.sin(fcs.beta(tbl['data']['tanb'])-fcs.alpha(tbl['data']['cosa']))
        proxer=tbl['data'].copy()
        proxer.insert(0,"sino", sino)
        proxer = proxer.sort_values(by=[param5],ascending=True)
        plt.scatter(proxer[param1]-proxer[param2],proxer[param3]-proxer[param4],c=proxer[param5])
    
    cbar=plt.colorbar(label=str_to_tex(param5))
    plt.xlabel(str_to_tex(param1).replace(' [GeV]','')+'-'+str_to_tex(param2), size=25)
    plt.xticks(size=20)
    #plt.xlim(-450,200)
    plt.ylabel(str_to_tex(param3).replace(' [GeV]','')+'-'+str_to_tex(param4), size=25)
    plt.yticks(size=20)
    #plt.ylim(-320,350)
    ax.grid()
    cbar.ax.tick_params(labelsize=20)  # Change size of tick labels
    cbar.ax.yaxis.label.set_size(20)   # Change size of colorbar label
    
    plt.savefig('../temp/Colormap_'+param5+THDM_type+strgl5+'-'+strga+'.png', dpi=600)
    
    plt.show()

    return 0

def plotter_diff_comp(param1,param2,param3,param4,param5,param6,*dataset):    

    fig, ax = plt.subplots(figsize=(10, 10))
    
    plt.title(latex_model +' '+ latex_alignment,size=20)
    
    for tbl in dataset:
        ratio = tbl['data'][param5]/tbl['data'][param6]
        proxer=tbl['data'].copy()
        proxer.insert(0,"ratio", ratio)
        proxer = proxer.sort_values(by=['ratio'],ascending=True)
        plt.scatter(proxer[param1]-proxer[param2],proxer[param3]-proxer[param4],c=proxer['ratio'])
    cbar=plt.colorbar(label=str_to_tex(param5)+'/'+str_to_tex(param6))
    plt.xlabel(str_to_tex(param1).replace(' [GeV]','')+'-'+str_to_tex(param2), size=25)
    plt.xticks(size=20)
    plt.xlim(-450,200)
    plt.ylabel(str_to_tex(param3).replace(' [GeV]','')+'-'+str_to_tex(param4), size=25)
    plt.yticks(size=20)
    plt.ylim(-320,350)
    ax.grid()
    cbar.ax.tick_params(labelsize=20)  # Change size of tick labels
    cbar.ax.yaxis.label.set_size(20)   # Change size of colorbar label
    
    plt.show()

    return 0

#%%%  Examples

# plotter_diff_color('mH','mHpm','mA','mHpm','kappa',Dataset_bsg)
# plotter_diff_color('mH','mHpm','mA','mHpm','tanb',Dataset_bsg)
# plotter_diff('mH','mHpm','mA','mHpm',Dataset_teo,Dataset_stu,Dataset_col,Dataset_bsg)
# plotter_diff_comp('mH','mHpm','mA','mHpm','kappa','kappa-kan',Dataset_bsg)
# plotter('mH','mHpm',Dataset_teo,Dataset_stu,Dataset_col,Dataset_bsg)
# plotter('mA','mHpm',Dataset_teo,Dataset_stu,Dataset_col,Dataset_bsg)
# plotter('mHpm','tanb',Dataset_teo,Dataset_stu,Dataset_col,Dataset_bsg)
# plotter('mHpm','M',Dataset_teo,Dataset_stu,Dataset_col,Dataset_bsg)
# plotter_3('mHpm','M','kappa',Dataset_bsg)
# plotter('M','tanb',Dataset_teo,Dataset_stu,Dataset_col,Dataset_bsg)
# plotter('M','kappa',Dataset_teo,Dataset_stu,Dataset_col,Dataset_bsg)
# plotter('tanb','kappa',Dataset_teo,Dataset_stu,Dataset_col,Dataset_bsg)
# plotter('mA','kappa',Dataset_teo,Dataset_stu,Dataset_col,Dataset_bsg)
# plotter('mH','kappa',Dataset_teo,Dataset_stu,Dataset_col,Dataset_bsg)
# plotter('mHpm','kappa',Dataset_teo,Dataset_stu,Dataset_col,Dataset_bsg)
# plotter('l1','kappa',Dataset_teo,Dataset_stu,Dataset_col,Dataset_bsg)
# plotter('l2','kappa',Dataset_teo,Dataset_stu,Dataset_col,Dataset_bsg)
# plotter('l3','kappa',Dataset_teo,Dataset_stu,Dataset_col,Dataset_bsg)
# plotter('l4','kappa',Dataset_teo,Dataset_stu,Dataset_col,Dataset_bsg)
# plotter('l5','kappa',Dataset_teo,Dataset_stu,Dataset_col,Dataset_bsg)

#%%%    

#TableTot_STU_Collid_BSG_unit = pd.concat([TableTot_STU_Collid_BSG,pd.DataFrame(np.abs(np.array(spr.calculate_quartics(TableTot_STU_Collid_BSG))).T,columns=['c93','c94','c102','c123','c140'])],axis=1)

# for cs in ['c93','c94','c102','c123','c140']:
#     cnd = spr.perturbative_unitarity_const(TableTot_STU_Collid_BSG_unit[cs])
#     TableTot_STU_Collid_BSG_unit = TableTot_STU_Collid_BSG_unit.drop(TableTot_STU_Collid_BSG_unit[cnd].index)

#TableTot_STU_Collid_BSG_unit = pd.concat([TableTot_STU_Collid_BSG_unit,pd.DataFrame(np.array(spr.calculate_eigenvalues(TableTot_STU_Collid_BSG_unit)).T,columns=['a0'])],axis=1)

#cnd = spr.perturbative_unitarity_const_a0(TableTot_STU_Collid_BSG_unit['a0'])
#TableTot_STU_Collid_BSG_unit = TableTot_STU_Collid_BSG_unit.drop(TableTot_STU_Collid_BSG_unit[cnd].index)
    

#%%

# plotter_diff_color('mH','mHpm','mA','mHpm','kappa',Dataset_unit)
# plotter_diff_color('mH','mHpm','mA','mHpm','tanb',Dataset_unit)
# plotter_diff('mH','mHpm','mA','mHpm',Dataset_teo,Dataset_stu,Dataset_col,Dataset_bsg,Dataset_unit)
# plotter_diff_comp('mH','mHpm','mA','mHpm','kappa','kappa-kan',Dataset_unit)
# plotter('mH','mHpm',Dataset_teo,Dataset_stu,Dataset_col,Dataset_bsg,Dataset_unit)
# plotter('mA','mHpm',Dataset_teo,Dataset_stu,Dataset_col,Dataset_bsg,Dataset_unit)
# plotter('mHpm','tanb',Dataset_teo,Dataset_stu,Dataset_col,Dataset_bsg,Dataset_unit)
# plotter('mHpm','M',Dataset_teo,Dataset_stu,Dataset_col,Dataset_bsg,Dataset_unit)
# plotter('M','tanb',Dataset_teo,Dataset_stu,Dataset_col,Dataset_bsg,Dataset_unit)
# plotter('M','kappa',Dataset_teo,Dataset_stu,Dataset_col,Dataset_bsg,Dataset_unit)
# plotter('tanb','kappa',Dataset_teo,Dataset_stu,Dataset_col,Dataset_bsg,Dataset_unit)
# plotter('mA','kappa',Dataset_teo,Dataset_stu,Dataset_col,Dataset_bsg,Dataset_unit)
# plotter('mH','kappa',Dataset_teo,Dataset_stu,Dataset_col,Dataset_bsg,Dataset_unit)
# plotter('mHpm','kappa',Dataset_teo,Dataset_stu,Dataset_col,Dataset_bsg,Dataset_unit)
# plotter('l1','kappa',Dataset_teo,Dataset_stu,Dataset_col,Dataset_bsg,Dataset_unit)
# plotter('l2','kappa',Dataset_teo,Dataset_stu,Dataset_col,Dataset_bsg,Dataset_unit)
# plotter('l3','kappa',Dataset_teo,Dataset_stu,Dataset_col,Dataset_bsg,Dataset_unit)
# plotter('l4','kappa',Dataset_teo,Dataset_stu,Dataset_col,Dataset_bsg,Dataset_unit)
# plotter('l5','kappa',Dataset_teo,Dataset_stu,Dataset_col,Dataset_bsg,Dataset_unit)

#%%%
# plotter('c93','kappa',Dataset_unit)
# plotter('c94','kappa',Dataset_unit)
# plotter('c102','kappa',Dataset_unit)
# plotter('c123','kappa',Dataset_unit)
# plotter('c140','kappa',Dataset_unit)

#%%

plotter_diff_color('mHpm','M','mH','M','kappa',Dataset_unit)
plotter_diff_color('mHpm','M','mA','M','kappa',Dataset_unit)
plotter_diff_color('mHpm','M','mH','M','tanb',Dataset_unit)
plotter_diff_color('mHpm','M','mA','M','tanb',Dataset_unit)

#%%

if load_CDF:

    plotter_diff('mH','mHpm','mA','mHpm',Dataset_stu,Dataset_stu_comb,Dataset_stu_cdf)
    plotter_diff('mH','mHpm','mA','mHpm',Dataset_stu_cdf)
    
#%%                 Recalculating CDF and Comb PU

# TableTot_STU_Collid_BSG_unit_CDF= pd.concat([TableTot_STU_Collid_BSG_CDF,pd.DataFrame(np.abs(np.array(spr.calculate_quartics(TableTot_STU_Collid_BSG_CDF))).T,columns=['c93','c94','c102','c123','c140'])],axis=1)

# for cs in ['c93','c94','c102','c123','c140']:
#     cnd = spr.perturbative_unitarity_const(TableTot_STU_Collid_BSG_unit[cs])
#     TableTot_STU_Collid_BSG_unit = TableTot_STU_Collid_BSG_unit.drop(TableTot_STU_Collid_BSG_unit[cnd].index)

# TableTot_STU_Collid_BSG_unit_CDF = pd.concat([TableTot_STU_Collid_BSG_unit_CDF,pd.DataFrame(np.array(spr.calculate_eigenvalues(TableTot_STU_Collid_BSG_unit_CDF)).T,columns=['a0'])],axis=1)

# cnd = spr.perturbative_unitarity_const_a0(TableTot_STU_Collid_BSG_unit_CDF['a0'])
# TableTot_STU_Collid_BSG_unit_CDF = TableTot_STU_Collid_BSG_unit_CDF.drop(TableTot_STU_Collid_BSG_unit_CDF[cnd].index)

#%%

# TableTot_STU_Collid_BSG_unit_Comb= pd.concat([TableTot_STU_Collid_BSG_Comb,pd.DataFrame(np.abs(np.array(spr.calculate_quartics(TableTot_STU_Collid_BSG_Comb))).T,columns=['c93','c94','c102','c123','c140'])],axis=1)

# for cs in ['c93','c94','c102','c123','c140']:
#     cnd = spr.perturbative_unitarity_const(TableTot_STU_Collid_BSG_unit[cs])
#     TableTot_STU_Collid_BSG_unit = TableTot_STU_Collid_BSG_unit.drop(TableTot_STU_Collid_BSG_unit[cnd].index)

# TableTot_STU_Collid_BSG_unit_Comb = pd.concat([TableTot_STU_Collid_BSG_unit_Comb,pd.DataFrame(np.array(spr.calculate_eigenvalues(TableTot_STU_Collid_BSG_unit_Comb)).T,columns=['a0'])],axis=1)

# cnd = spr.perturbative_unitarity_const_a0(TableTot_STU_Collid_BSG_unit_Comb['a0'])
# TableTot_STU_Collid_BSG_unit_Comb = TableTot_STU_Collid_BSG_unit_Comb.drop(TableTot_STU_Collid_BSG_unit_Comb[cnd].index)

    