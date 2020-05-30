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
from vasp_dos import VASP_DOS
from vasp_dos import get_all_VASP_files

GCNList = []
atom_type = []
band_list = []
occupied_band_list = []
DOSCAR_files, CONTCAR_files = get_all_VASP_files(\
        r'C:\Users\lansf\Documents\Data\PROBE_PDOS\vasp_dos_files\Pt_nano')


for DOSCAR, CONTCAR in zip(DOSCAR_files, CONTCAR_files):
    nanoparticle = read(CONTCAR)
    CN = Coordination(nanoparticle,cutoff=1.25)
    CN.get_coordination_numbers()
    # read and return densityofstates object
    PDOS = VASP_DOS(DOSCAR)   
    for atom_index in range(len(nanoparticle)):
        GCNList.append(CN.get_gcn([atom_index]))
        if CN.cn[atom_index] < 12:
            atom_type.append('surface')
        else:
            atom_type.append('bulk')
        band_center = PDOS.get_band_center(atom_index, ['s','p','d']\
                                , sum_density=True) - PDOS.e_fermi
        occupied_band_center = PDOS.get_band_center(atom_index, ['s','p','d']\
                                , sum_density=True, max_energy=PDOS.e_fermi) - PDOS.e_fermi
            
        band_list.append(band_center)
        occupied_band_list.append(occupied_band_center)

GCNList = np.array(GCNList)
atom_type = np.array(atom_type)
band_list = np.array(band_list).T
occupied_band_list = np.array(occupied_band_list).T

#plotting scaling of band center with GCN
plt.figure(0,figsize=(7,5))
colors = ['b','g','r']
Efit = []
for count, color in enumerate(colors):
    Efit.append(np.polyfit(GCNList,band_list[count], 1))
    plt.plot(np.sort(GCNList), np.poly1d(Efit[count])(np.sort(GCNList)), color + '--')
for count, color in enumerate(colors):
    plt.plot(GCNList, band_list[count], color + 'o')
plt.legend([r'${\epsilon}_{s}$=%.2fGCN + %.2f eV' %(Efit[0][0],Efit[0][1])
,r'${\epsilon}_{p}$=%.2fGCN + %.2f eV' %(Efit[1][0],Efit[1][1])
,r'${\epsilon}_{d}$=%.2fGCN + %.2f eV' %(Efit[2][0],Efit[2][1])]
,loc=3,prop={'size':14},frameon=False)
plt.xlabel('Generalized coordination number (GCN)')
plt.ylabel('Band center (${\epsilon}$) [eV]')
plt.show()

#plotting scaling of occupied band center with GCN
plt.figure(1,figsize=(7,5))
colors = ['b','g','r']
Efit = []
for count, color in enumerate(colors):
    Efit.append(np.polyfit(GCNList, occupied_band_list[count], 1))
    plt.plot(np.sort(GCNList), np.poly1d(Efit[count])(np.sort(GCNList)), color + '--')
for count, color in enumerate(colors):
    plt.plot(GCNList, occupied_band_list[count], color + 'o')
plt.legend([r'${\epsilon}_{s}$=%.2fGCN + %.2f eV' %(Efit[0][0],Efit[0][1])
,r'${\epsilon}_{p}$=%.2fGCN + %.2f eV' %(Efit[1][0],Efit[1][1])
,r'${\epsilon}_{d}$=%.2fGCN + %.2f eV' %(Efit[2][0],Efit[2][1])]
,loc=3,prop={'size':14},frameon=False)
plt.xlabel('Generalized coordination number (GCN)')
plt.ylabel('Occupied band center (${\epsilon}^{*}$) [eV]')
plt.show()
#plt.xlim([1,12.5])
#plt.ylim([-6,1.5])

#plotting scaling of band center with GCN for surface sites
plt.figure(2,figsize=(7,5))
colors = ['b','g','r']
Efit = []
for count, color in enumerate(colors):
    Efit.append(np.polyfit(GCNList[atom_type=='surface']\
                ,band_list[count][atom_type=='surface'], 1))
    plt.plot(np.sort(GCNList[atom_type=='surface'])\
             , np.poly1d(Efit[count])\
             (np.sort(GCNList[atom_type=='surface'])), color + '--')
for count, color in enumerate(colors):
    plt.plot(GCNList[atom_type=='surface'], band_list[count][atom_type=='surface'], color + 'o')
plt.legend([r'${\epsilon}_{s}$=%.2fGCN + %.2f eV' %(Efit[0][0],Efit[0][1])
,r'${\epsilon}_{p}$=%.2fGCN + %.2f eV' %(Efit[1][0],Efit[1][1])
,r'${\epsilon}_{d}$=%.2fGCN + %.2f eV' %(Efit[2][0],Efit[2][1])]
,loc=3,prop={'size':14},frameon=False)
plt.xlabel('Generalized coordination number (GCN)')
plt.ylabel('Band center (${\epsilon}$) [eV]')
plt.show()

#plotting scaling of occupied band center with GCN for surface sites
plt.figure(3,figsize=(7,5))
colors = ['b','g','r']
Efit = []
for count, color in enumerate(colors):
    Efit.append(np.polyfit(GCNList[atom_type=='surface']\
                ,occupied_band_list[count][atom_type=='surface'], 1))
    plt.plot(np.sort(GCNList[atom_type=='surface'])\
             , np.poly1d(Efit[count])\
             (np.sort(GCNList[atom_type=='surface'])), color + '--')
for count, color in enumerate(colors):
    plt.plot(GCNList[atom_type=='surface'], occupied_band_list[count][atom_type=='surface'], color + 'o')
plt.legend([r'${\epsilon}_{s}$=%.2fGCN + %.2f eV' %(Efit[0][0],Efit[0][1])
,r'${\epsilon}_{p}$=%.2fGCN + %.2f eV' %(Efit[1][0],Efit[1][1])
,r'${\epsilon}_{d}$=%.2fGCN + %.2f eV' %(Efit[2][0],Efit[2][1])]
,loc=3,prop={'size':14},frameon=False)
plt.xlabel('Generalized coordination number (GCN)')
plt.ylabel('Occupied band center (${\epsilon}^{*}$) [eV]')
plt.show()


#calculating f factor
plt.figure(4,figsize=(7,5))
colors = ['b','g','r']
for count, color in enumerate(colors):
    if count in [0,2]:
        plt.plot(GCNList[atom_type=='surface'], occupied_band_list[count][atom_type=='surface']/band_list[count][atom_type=='surface'], color + 'o')
for count, color in enumerate(colors):
    if count in [0,2]:
        plt.plot(GCNList[atom_type=='bulk'], occupied_band_list[count][atom_type=='bulk']/band_list[count][atom_type=='bulk'], color + 's')
plt.legend(['s/surface', 'd/surface', 's/bulk', 'd/bulk']
,loc=2,title='Band/Location',prop={'size':14},frameon=False)
plt.gcf().subplots_adjust(bottom=0.15)
plt.gcf().subplots_adjust(top=0.99)
plt.gcf().subplots_adjust(left=0.10)
plt.gcf().subplots_adjust(right=0.99)
plt.xlabel('Generalized coordination number (GCN)')
plt.ylabel(r'${\epsilon}^{*}$ $\times$ $({\epsilon})^{-1}$')

#plotting occupied vs total band center
plt.figure(5,figsize=(7,5))
colors = ['b','g','r']
for count, color in enumerate(colors):
    if count in [0,2]:
        plt.plot(band_list[count][atom_type=='surface'], occupied_band_list[count][atom_type=='surface'], color + 'o')
for count, color in enumerate(colors):
    if count in [0,2]:
        plt.plot(band_list[count][atom_type=='bulk'], occupied_band_list[count][atom_type=='bulk'], color + 's')
plt.legend(['s/surface', 'd/surface', 's/bulk', 'd/bulk']
,loc=2,title='Band/Location',prop={'size':14},frameon=False)
plt.gcf().subplots_adjust(bottom=0.15)
plt.gcf().subplots_adjust(top=0.99)
plt.gcf().subplots_adjust(left=0.10)
plt.gcf().subplots_adjust(right=0.99)
plt.xlabel('Band center (${\epsilon}$) [eV]')
plt.ylabel('Occupied band center (${\epsilon}^{*}$) [eV]')
