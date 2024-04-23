from anyBSM import anyBSM

import numpy as np
import matplotlib.pyplot as plt

sm = anyBSM(
    'SM', # load the singlet extended SM
    scheme_name = "OS", # load the OS scheme (tadpoles renormalised MSbar)
    progress = True, # dont show any progress bars (default = True)
    caching = 2 # enable the cache (default = 2)
)

smM = anyBSM(
    'SM', # load the singlet extended SM
    scheme_name = "MS", # load the MS scheme (tadpoles renormalised MSbar)
    progress = True, # dont show any progress bars (default = True)
    caching = 2 # enable the cache (default = 2)
)


#sm.load_renormalization_scheme("MS")
#print(sm.lambdahhh())

#sm.set_evaluation_mode('analytical')
#print(sm.lambdahhh())

'''
MH= np.linspace(start=120,stop=600, num=10)
lamh = [sm.lambdahhh(parameters={'Mh': M})['total'].real/smM.lambdahhh(parameters={'Mh': M})['total'].real for M in MH]
lamTree = [sm.lambdahhh(parameters={'Mh': M})['treelevel'].real/smM.lambdahhh(parameters={'Mh': M})['treelevel'].real for M in MH]
lamOne = [sm.lambdahhh(parameters={'Mh': M})['genuine'].real/smM.lambdahhh(parameters={'Mh': M})['genuine'].real for M in MH]
lamWFR = [sm.lambdahhh(parameters={'Mh': M})['wfr'].real/smM.lambdahhh(parameters={'Mh': M})['wfr'].real for M in MH]
lamTad = [sm.lambdahhh(parameters={'Mh': M})['tads'].real/smM.lambdahhh(parameters={'Mh': M})['tads'].real for M in MH] 
lamMass = [smM.lambdahhh(parameters={'Mh': M})['massren'].real for M in MH]
lamVev = [smM.lambdahhh(parameters={'Mh': M})['vevren'].real for M in MH]

fig, ax = plt.subplots()
ax.plot(MH,lamh, label = "Total")
ax.plot(MH,lamTree, label = "Tree")
ax.plot(MH,lamOne, label = "One")
ax.plot(MH,lamWFR, label = "WFR")
ax.plot(MH,lamTad, label = "Tadpole")
ax.plot(MH,lamMass, label = "Mass")
ax.plot(MH,lamVev, label = "VeV")

ax.set_xlabel('$M_h$ [GeV]')
ax.set_ylabel(r'${\lambda_{hhh}^{\mathrm{SM}}}$')
ax.margins(0,0)

plt.legend()
plt.show()
'''

print(sm.lambdahhh())
print(smM.lambdahhh())
print(sm.MW())
print(smM.MW())

############################ top #
'''
MH= np.linspace(start=120.,stop=220, num=10)
lamh = [sm.lambdahhh(parameters={'Mu3': M})['total'].real for M in MH]

lamTree = [sm.lambdahhh(parameters={'Mu3': M})['treelevel'].real for M in MH]
lamOne = [smM.lambdahhh(parameters={'Mu3': M})['genuine'].real for M in MH]
lamWFR = [sm.lambdahhh(parameters={'Mu3': M})['wfr'].real for M in MH]
lamTad = [sm.lambdahhh(parameters={'Mu3': M})['tads'].real for M in MH] 
lamMass = [smM.lambdahhh(parameters={'Mu3': M})['massren'].real for M in MH]
lamVev = [smM.lambdahhh(parameters={'Mu3': M})['vevren'].real for M in MH]

[smM.lambdahhh(parameters={'Mu3': M})['total'].real for M in MH]


print(sm.lambdahhh())
print(smM.lambdahhh())
print(sm.MW())
print(smM.MW())


mW = 80.379
MHH = 125.1
Mt = 172.5
Mx = MHH

#lamThTree=[3*MHH**2/(6.60454*mW*np.sqrt(1 - 0.000120264*mW**2)) for M in MH]
#lamThOne=[-3*MHH**2/(6.60454*sm.MW(parameters={'Mu3': M})*np.sqrt(1 - 0.000120264*sm.MW(parameters={'Mu3': M})**2))*M**4/(np.pi**2*MHH**2*(6.60454*sm.MW(parameters={'Mu3': M})*np.sqrt(1 - 0.000120264*sm.MW(parameters={'Mu3': M})**2))**2) for M in MH]
lamThTotal=[3*Mx**2/(6.60454*sm.MW(parameters={'Mh':125.1, 'Mu3': 172.5})*np.sqrt(1 - 0.000120264*sm.MW(parameters={'Mh':125.1,'Mu3': 172.5})**2))-3*Mx**2/(6.60454*sm.MW(parameters={'Mh':MHH,'Mu3': Mt})*np.sqrt(1 - 0.000120264*sm.MW(parameters={'Mh':MHH,'Mu3': Mt})**2))*M**4/(np.pi**2*Mx**2*(6.60454*sm.MW(parameters={'Mh':MHH, 'Mu3': Mt})*np.sqrt(1 - 0.000120264*sm.MW(parameters={'Mh':MHH,'Mu3': Mt})**2))**2) for M in MH]

print(MH)
print(lamThTotal)
print(lamh)
print([m/n for m, n in zip(lamh, lamThTotal)])
print((6.60454*sm.MW(parameters={'Mh':MHH, 'Mu3': Mt})*np.sqrt(1 - 0.000120264*sm.MW(parameters={'Mh':MHH,'Mu3': Mt})**2)))
print((6.60454*sm.MW(parameters={'Mh':125.1, 'Mu3': 172.5})*np.sqrt(1 - 0.000120264*sm.MW(parameters={'Mh':125.1,'Mu3': 172.5})**2)))



fig, ax = plt.subplots()
ax.plot(MH,lamh, label = "Total")
#ax.plot(MH,lamTree, label = "Tree")
#ax.plot(MH,lamOne, label = "One")
#ax.plot(MH,lamWFR, label = "WFR")
#ax.plot(MH,lamTad, label = "Tadpole")
#ax.plot(MH,lamMass, label = "Mass")
#ax.plot(MH,lamVev, label = "VeV")

#ax.plot(MH,lamThTree, label = "Tree Th")
#ax.plot(MH,lamThOne, label = "One Th")
ax.plot(MH,lamThTotal, label = "Total Th")

ax.set_xlabel('$M_t$ [GeV]')
ax.set_ylabel(r'${\lambda_{hhh}^{\mathrm{SM}}}$')
ax.margins(0,0)

plt.legend()
plt.show()

ratio = [m/n for m, n in zip(lamh, lamThTotal)]

fig, ax = plt.subplots()
ax.plot(MH,ratio, label = "Total/Total Th")

ax.set_xlabel('$M_t$ [GeV]')
ax.set_ylabel(r'${\lambda_{hhh}^{\mathrm{SM}}}$')
ax.margins(0,0)

plt.legend()
plt.show()
'''
############################ higgs #
'''
print(sm.lambdahhh(parameters={'Mu3': Mt}))
print(smM.lambdahhh(parameters={'Mu3': Mt}))
print(sm.MW(parameters={'Mu3': Mt}))
print(smM.MW(parameters={'Mu3': Mt}))
'''
MH= np.linspace(start=120.,stop=600, num=10)
lamh = [sm.lambdahhh(parameters={'Mh': M})['total'].real for M in MH]

lamTree = [sm.lambdahhh(parameters={'Mh': M})['treelevel'].real for M in MH]
lamOne = [sm.lambdahhh(parameters={'Mh': M})['genuine'].real for M in MH]
lamWFR = [sm.lambdahhh(parameters={'Mh': M})['wfr'].real for M in MH]
lamTad = [sm.lambdahhh(parameters={'Mh': M})['tads'].real for M in MH] 
lamMass = [sm.lambdahhh(parameters={'Mh': M})['massren'].real for M in MH]
lamVev = [sm.lambdahhh(parameters={'Mh': M})['vevren'].real for M in MH]

[smM.lambdahhh(parameters={'Mh': M})['total'].real for M in MH]


print(sm.lambdahhh())
print(smM.lambdahhh())
print(sm.MW())
print(smM.MW())


mW = 80.379
MHH = 125.1
Mt = 172.5
Mx = MHH

vSM = (6.60454*mW*np.sqrt(1 - 0.000120264*mW**2))

lamThTree=[3*Mx**2/(6.60454*mW*np.sqrt(1 - 0.000120264*mW**2)) for Mx in MH]
lamThOne=[-(3*Mx**2/vSM)*Mt**4/(np.pi**2*Mx**2*vSM**2) + (3*Mx**2/vSM)*(7/2)*Mt**2/(16*np.pi**2*vSM**2) + (3*Mx**2/vSM)*6*Mx**4/(12*np.pi**2*Mx**2*vSM**2) for Mx in MH]
lamThTotal=[3*Mx**2/vSM-(3*Mx**2/vSM)*Mt**4/(np.pi**2*Mx**2*vSM**2) + (3*Mx**2/vSM)*(7/2)*Mt**2/(16*np.pi**2*vSM**2) + (3*Mx**2/vSM)*6*Mx**4/(12*np.pi**2*Mx**2*vSM**2) for Mx in MH]

print(MH)
print(lamThTotal)
print(lamh)
print([m/n for m, n in zip(lamh, lamThTotal)])
print((6.60454*sm.MW()*np.sqrt(1 - 0.000120264*sm.MW()**2)))
print((6.60454*sm.MW(parameters={'Mh':125.1, 'Mu3': 172.5})*np.sqrt(1 - 0.000120264*sm.MW(parameters={'Mh':125.1,'Mu3': 172.5})**2)))

lamRestOS= [m1+m2+m3+m4+m5 for m1,m2,m3,m4,m5 in zip(lamOne, lamWFR, lamTad, lamMass, lamVev)]
lamRestMS= [m1+m2+m3 for m1,m2,m3 in zip(lamOne, lamWFR, lamTad)]


fig, ax = plt.subplots()
ax.plot(MH,lamh, label = "Total")
ax.plot(MH,lamTree, label = "Tree")
#ax.plot(MH,lamOne, label = "One")
#ax.plot(MH,lamWFR, label = "WFR")
#ax.plot(MH,lamTad, label = "Tadpole")
#ax.plot(MH,lamMass, label = "Mass")
#ax.plot(MH,lamVev, label = "VeV")
ax.plot(MH,lamRestOS, label = "Rest OS")
ax.plot(MH,lamRestMS, label = "Rest MS")



ax.plot(MH,lamThTree, label = "Tree Th")
ax.plot(MH,lamThOne, label = "One Th")
ax.plot(MH,lamThTotal, label = "Total Th")

ax.set_xlabel('$M_h$ [GeV]')
ax.set_ylabel(r'${\lambda_{hhh}^{\mathrm{SM}}}$')
ax.margins(0,0)

plt.legend()
plt.show()

ratio = [m/n for m, n in zip(lamh, lamThTotal)]

fig, ax = plt.subplots()
ax.plot(MH,ratio, label = "Total/Total Th")

ax.set_xlabel('$M_h$ [GeV]')
ax.set_ylabel(r'${\lambda_{hhh}^{\mathrm{SM}}}$')
ax.margins(0,0)

plt.legend()
plt.show()
