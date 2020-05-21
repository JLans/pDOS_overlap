#!python2
# -*- coding: utf-8 -*-

"""
Created on Wed Jan 11 14:33:54 2017

@author: lansford
"""
from __future__ import division
from vasp_dos.coordination import Coordination
from ase.io import read
import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import trapz
from matplotlib import rcParams
from vasp_dos import VASP_DOS
import os
rcParams['lines.markersize'] = 10
NumPt = [147,19,44,79,6]
GCNList = []
EdList = []
dfillingList = []
EsList = []
sfillingList = []
EpList = []
pfillingList = []
EspList = []
spfillingList = []
SurfaceAtom = []
for i in NumPt:
    folder = os.path.expanduser('C:/Users/lansf/Documents/Data/DOS/spin/Pt'+str(i) + '/')
    nanoparticle = read(folder+'CONTCAR')
    # read and return densityofstates object
    dos = VASP_DOS(folder + 'DOSCAR')
    # Get energy grid
    energies = dos.get_energies()
    fermi = dos.e_fermi
    #why cutting off at energies+fermi
    #might not need to add back eferim
    idx = (np.abs(energies+fermi)).argmin()
    #idx = len(energies)
    CN = Coordination(nanoparticle,cutoff=1.25)
    CN.get_coordination_numbers()
    for Atomindex in range(i):
        if CN.cn[Atomindex] < 12:
            SurfaceAtom.append('surface')
        else:
            SurfaceAtom.append('bulk')
        GCNList.append(CN.get_gcn([Atomindex]))
        d_atom = dos.get_atom_dos(Atomindex)
        d_atom = dos.get_site_dos(Atomindex,['s','p','d'], sum_density=True)
        s = d_atom[0]
        p = d_atom[1]
        dband = d_atom[3]
        Energy = trapz(dband[0:idx]*energies[0:idx],energies[0:idx])
        d_filling = trapz(dband[0:idx],energies[0:idx])
        Ed = Energy/d_filling
        EdList.append(Ed)
        dfillingList.append(d_filling)
        Energy = trapz(s[0:idx]*energies[0:idx],energies[0:idx])
        s_filling = trapz(s[0:idx],energies[0:idx])
        Es = Energy/s_filling
        EsList.append(Es)
        sfillingList.append(s_filling)
        Energy = trapz(p[0:idx]*energies[0:idx],energies[0:idx])
        p_filling = trapz(p[0:idx],energies[0:idx])
        Ep = Energy/p_filling
        EpList.append(Ep)
        pfillingList.append(p_filling)
        
        Energy = trapz((s[0:idx]+p[0:idx])*energies[0:idx],energies[0:idx])
        sp_filling = trapz((s[0:idx]+p[0:idx]),energies[0:idx])
        Esp = Energy/sp_filling
        EspList.append(Esp)
        spfillingList.append(sp_filling)
        
GCNList = []
EdListfilled = []
dfillingList = []
EsListfilled = []
sfillingList = []
EpListfilled = []
pfillingList = []
SurfaceAtom = []
for i in NumPt:
    folder = os.path.expanduser('C:/Users/lansf/Documents/Data/DOS/spin/Pt'+str(i) + '/')
    nanoparticle = read(folder+'CONTCAR')
    # read and return densityofstates object
    # Get energy grid
    energies = dos.energies
    fermi = dos.e_fermi
    idx = (np.abs(energies)).argmin()
    CN = Coordination(nanoparticle,cutoff=1.25)
    CN.get_coordination_numbers()
    for Atomindex in range(i):
        if CN.cn[Atomindex] < 12:
            SurfaceAtom.append('surface')
        else:
            SurfaceAtom.append('bulk')
        GCNList.append(CN.get_gcn([Atomindex]))
        d_atom = dos.get_atom_dos(Atomindex)
        d_atom = dos.get_site_dos(Atomindex,['s','p','d'], sum_density=True)
        s = d_atom[0]
        p = d_atom[1]
        dband = d_atom[3]
        Energy = trapz(dband[0:idx]*energies[0:idx],energies[0:idx])
        d_filling = trapz(dband[0:idx],energies[0:idx])
        Ed = Energy/d_filling
        EdListfilled.append(Ed)
        dfillingList.append(d_filling)
        Energy = trapz(s[0:idx]*energies[0:idx],energies[0:idx])
        s_filling = trapz(s[0:idx],energies[0:idx])
        Es = Energy/s_filling
        EsListfilled.append(Es)
        sfillingList.append(s_filling)
        Energy = trapz(p[0:idx]*energies[0:idx],energies[0:idx])
        p_filling = trapz(p[0:idx],energies[0:idx])
        Ep = Energy/p_filling
        EpListfilled.append(Ep)
        pfillingList.append(p_filling)

GCNList = np.array(GCNList)
EsList = np.array(EsList)
EpList = np.array(EpList)
EspList = np.array(EspList)
EdList = np.array(EdList)
EsListfilled = np.array(EsListfilled)
EpListfilled = np.array(EpListfilled)
EdListfilled = np.array(EdListfilled)
sfillingList = np.array(sfillingList)
pfillingList = np.array(pfillingList)
spfillingList = np.array(spfillingList)
dfillingList = np.array(dfillingList)
SurfaceAtom = np.array(SurfaceAtom)

plt.figure(0,figsize=(7,5))
Esfit = np.polyfit(GCNList, EsList, 1)
Epfit = np.polyfit(GCNList, EpList, 1)
Edfit = np.polyfit(GCNList, EdList, 1)
plt.plot(GCNList,EsList,'bo',GCNList,EpList,'go',GCNList,EdList,'ro')
plt.plot(np.sort(GCNList), np.poly1d(Esfit)(np.sort(GCNList)),'--b')
plt.plot(np.sort(GCNList), np.poly1d(Epfit)(np.sort(GCNList)),'--g')
plt.plot(np.sort(GCNList), np.poly1d(Edfit)(np.sort(GCNList)),'--r')
plt.legend([r'${\epsilon}_{s}$=%.2fGCN + %.2f eV' %(Esfit[0],Esfit[1])
,r'${\epsilon}_{p}$=%.2fGCN + %.2f eV' %(Epfit[0],Epfit[1])
,r'${\epsilon}_{d}$=%.2fGCN + %.2f eV' %(Edfit[0],Edfit[1])]
,loc=3,prop={'size':14},frameon=False)
plt.xlabel('Generalized Coordination Number (GCN)')
plt.ylabel('Band Center (${\epsilon}$) [eV]')
plt.gcf().subplots_adjust(bottom=0.15)
plt.gcf().subplots_adjust(top=0.99)
plt.gcf().subplots_adjust(left=0.12)
plt.gcf().subplots_adjust(right=0.99)
plt.xlim([1,12.5])
plt.ylim([-6,1.5])

plt.figure(1,figsize=(7,5))
EsfitS = np.polyfit(GCNList[SurfaceAtom=='surface'], EsList[SurfaceAtom=='surface'], 1)
EpfitS = np.polyfit(GCNList[SurfaceAtom=='surface'], EpList[SurfaceAtom=='surface'], 1)
EdfitS = np.polyfit(GCNList[SurfaceAtom=='surface'], EdList[SurfaceAtom=='surface'], 1)
plt.plot(GCNList[SurfaceAtom=='surface'],EsList[SurfaceAtom=='surface'],'bo',GCNList[SurfaceAtom=='surface']
,EpList[SurfaceAtom=='surface'],'go',GCNList[SurfaceAtom=='surface'],EdList[SurfaceAtom=='surface'],'ro')
plt.plot(np.sort(GCNList[SurfaceAtom=='surface']), np.poly1d(EsfitS)(np.sort(GCNList[SurfaceAtom=='surface'])),'--b')
plt.plot(np.sort(GCNList[SurfaceAtom=='surface']), np.poly1d(EpfitS)(np.sort(GCNList[SurfaceAtom=='surface'])),'--g')
plt.plot(np.sort(GCNList[SurfaceAtom=='surface']), np.poly1d(EdfitS)(np.sort(GCNList[SurfaceAtom=='surface'])),'--r')
plt.legend([r'${\epsilon}_{s}$=%.2fGCN + %.2f eV' %(EsfitS[0],EsfitS[1])
,r'${\epsilon}_{p}$=%.2fGCN + %.2f eV' %(EpfitS[0],EpfitS[1])
,r'${\epsilon}_{d}$=%.2fGCN + %.2f eV' %(EdfitS[0],EdfitS[1])]
,loc=3,prop={'size':14},frameon=False)
plt.xlabel('Generalized Coordination Number (GCN)')
plt.ylabel('Band Center (${\epsilon}$) [eV]')
plt.gcf().subplots_adjust(bottom=0.15)
plt.gcf().subplots_adjust(top=0.99)
plt.gcf().subplots_adjust(left=0.12)
plt.gcf().subplots_adjust(right=0.99)
plt.xlim([1,7])
plt.ylim([-5,1.5])

plt.figure(5,figsize=(7,5))
EsfitS = np.polyfit(GCNList[SurfaceAtom=='surface'], EsList[SurfaceAtom=='surface'], 1)
EspfitS = np.polyfit(GCNList[SurfaceAtom=='surface'], EspList[SurfaceAtom=='surface'], 1)
EdfitS = np.polyfit(GCNList[SurfaceAtom=='surface'], EdList[SurfaceAtom=='surface'], 1)
plt.plot(GCNList[SurfaceAtom=='surface'],EsList[SurfaceAtom=='surface'],'bo',GCNList[SurfaceAtom=='surface']
,EspList[SurfaceAtom=='surface'],'go',GCNList[SurfaceAtom=='surface'],EdList[SurfaceAtom=='surface'],'ro')
plt.plot(np.sort(GCNList[SurfaceAtom=='surface']), np.poly1d(EsfitS)(np.sort(GCNList[SurfaceAtom=='surface'])),'--b')
plt.plot(np.sort(GCNList[SurfaceAtom=='surface']), np.poly1d(EspfitS)(np.sort(GCNList[SurfaceAtom=='surface'])),'--g')
plt.plot(np.sort(GCNList[SurfaceAtom=='surface']), np.poly1d(EdfitS)(np.sort(GCNList[SurfaceAtom=='surface'])),'--r')
plt.legend([r'${\epsilon}_{s}$=%.2fGCN + %.2f eV' %(EsfitS[0],EsfitS[1])
,r'${\epsilon}_{sp}$=%.2fGCN + %.2f eV' %(EspfitS[0],EspfitS[1])
,r'${\epsilon}_{d}$=%.2fGCN + %.2f eV' %(EdfitS[0],EdfitS[1])]
,loc=3,prop={'size':14},frameon=False)
plt.xlabel('Generalized Coordination Number (GCN)')
plt.ylabel('Band Center (${\epsilon}$) [eV]')
plt.gcf().subplots_adjust(bottom=0.15)
plt.gcf().subplots_adjust(top=0.99)
plt.gcf().subplots_adjust(left=0.12)
plt.gcf().subplots_adjust(right=0.99)
plt.xlim([1,7])
plt.ylim([-5,1.5])

plt.figure(2,figsize=(9,5))
rcParams['lines.markersize'] = 8
plt.plot(EsList[SurfaceAtom=='surface'],EsListfilled[SurfaceAtom=='surface'],'bo',zorder=9)
plt.plot(EsList[SurfaceAtom=='bulk'],EsListfilled[SurfaceAtom=='bulk'],'bs',zorder=10)
plt.plot(EpList[SurfaceAtom=='surface'],EpListfilled[SurfaceAtom=='surface'],'go',zorder=3)
plt.plot(EpList[SurfaceAtom=='bulk'],EpListfilled[SurfaceAtom=='bulk'],'gs',zorder=4)
plt.plot(EdList[SurfaceAtom=='surface'],EdListfilled[SurfaceAtom=='surface'],'ro',zorder=5)
plt.plot(EdList[SurfaceAtom=='bulk'],EdListfilled[SurfaceAtom=='bulk'],'rs',zorder=6)
plt.legend(['${\epsilon}_{s}$/surface','${\epsilon}_{s}$/bulk'
            ,'${\epsilon}_{p}$/surface','${\epsilon}_{p}$/bulk'
            ,'${\epsilon}_{d}$/surface','${\epsilon}_{d}$/bulk']
,loc=2,title='Band/Location',prop={'size':14},frameon=False)
plt.gcf().subplots_adjust(bottom=0.15)
plt.gcf().subplots_adjust(top=0.99)
plt.gcf().subplots_adjust(left=0.10)
plt.gcf().subplots_adjust(right=0.99)
plt.xlabel('Total Band Center (${\epsilon}$) [eV]')
plt.ylabel('Occupied Band Center (${\epsilon}^{*}$) [eV]')

plt.figure(3,figsize=(7,5))
EsfitS = np.polyfit(GCNList[SurfaceAtom=='surface'], EsListfilled[SurfaceAtom=='surface'], 1)
EpfitS = np.polyfit(GCNList[SurfaceAtom=='surface'], EpListfilled[SurfaceAtom=='surface'], 1)
EdfitS = np.polyfit(GCNList[SurfaceAtom=='surface'], EdListfilled[SurfaceAtom=='surface'], 1)
plt.plot(GCNList[SurfaceAtom=='surface'],EsListfilled[SurfaceAtom=='surface'],'bo',GCNList[SurfaceAtom=='surface']
,EpListfilled[SurfaceAtom=='surface'],'go',GCNList[SurfaceAtom=='surface'],EdListfilled[SurfaceAtom=='surface'],'ro')
plt.plot(np.sort(GCNList[SurfaceAtom=='surface']), np.poly1d(EsfitS)(np.sort(GCNList[SurfaceAtom=='surface'])),'--b')
plt.plot(np.sort(GCNList[SurfaceAtom=='surface']), np.poly1d(EpfitS)(np.sort(GCNList[SurfaceAtom=='surface'])),'--g')
plt.plot(np.sort(GCNList[SurfaceAtom=='surface']), np.poly1d(EdfitS)(np.sort(GCNList[SurfaceAtom=='surface'])),'--r')
plt.legend([r'${\epsilon}_{s}$=%.2fGCN + %.2f eV' %(EsfitS[0],EsfitS[1])
,r'${\epsilon}_{p}$=%.2fGCN + %.2f eV' %(EpfitS[0],EpfitS[1])
,r'${\epsilon}_{d}$=%.2fGCN + %.2f eV' %(EdfitS[0],EdfitS[1])],loc=3)
"""Getting predicted slope for surface"""
f = np.mean(EsListfilled[SurfaceAtom=='surface']/EsList[SurfaceAtom=='surface'])
Theta = np.mean(sfillingList)
Ecohesive = -5.85
Spredicted = Ecohesive/(2*12*f*Theta)
print(Spredicted)