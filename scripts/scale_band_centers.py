#!python2
# -*- coding: utf-8 -*-

"""
Created on Wed Jan 11 14:33:54 2017

@author: lansford
"""
from __future__ import division
from pdos_overlap.coordination import get_geometric_data
import numpy as np
import matplotlib.pyplot as plt
from pdos_overlap.vasp_dos import VASP_DOS
from pdos_overlap.vasp_dos import get_all_VASP_files
from pdos_overlap.plotting_tools import set_figure_settings
set_figure_settings('paper')

GCNList = []
atom_type = []
band_list = []
band_width_list = []
occupied_band_list = []
unoccupied_band_list = []
filling_list = []
second_moment_list = []
bond_energy_list = []
DOSCAR_files, CONTCAR_files = get_all_VASP_files(\
        r'C:\Users\lansf\Documents\Data\PROBE_PDOS\vasp_dos_files\Pt_nano')


for DOSCAR, CONTCAR in zip(DOSCAR_files, CONTCAR_files):
    indices, GCNs, atom_types = get_geometric_data(CONTCAR)
    GCNList += GCNs.tolist()
    atom_type += atom_types.tolist()
    # read and return densityofstates object
    PDOS = VASP_DOS(DOSCAR)   
    for atom_index in indices:
        band_center = PDOS.get_band_center(atom_index, ['s','p','d']\
                                , sum_density=True) - PDOS.e_fermi
        occupied_band_center = PDOS.get_band_center(atom_index, ['s','p','d']\
                                , sum_density=True, max_energy=PDOS.e_fermi) - PDOS.e_fermi
        unoccupied_band_center = PDOS.get_band_center(atom_index, ['s','p','d']\
                                , sum_density=True, min_energy=PDOS.e_fermi) - PDOS.e_fermi
        band_width = PDOS.get_center_width(PDOS.e_fermi, atom_index, ['s','p','d']\
                                , sum_density=True)
        second_moment = PDOS.get_second_moment(atom_index, ['s','p','d']\
                                , sum_density=True)
        
        bond_energy = PDOS.get_bond_energy(atom_index, ['s','p','d']\
                                , sum_density=True)
            
        filling = PDOS.get_filling(atom_index, ['s','p','d']\
                                , sum_density=True, max_energy=PDOS.e_fermi)
            
        band_list.append(band_center)
        band_width_list.append(band_width)
        occupied_band_list.append(occupied_band_center)
        unoccupied_band_list.append(unoccupied_band_center)
        filling_list.append(filling)
        second_moment_list.append(second_moment)
        bond_energy_list.append(bond_energy)

GCNList = np.array(GCNList)
atom_type = np.array(atom_type)
band_list = np.array(band_list).T
band_width_list = np.array(band_width_list).T
occupied_band_list = np.array(occupied_band_list).T
unoccupied_band_list = np.array(unoccupied_band_list).T
filling_list = np.array(filling_list).T
second_moment_list = np.array(second_moment_list).T
bond_energy_list = np.array(bond_energy_list).T
band_list
#plotting scaling of band center with GCN
plt.figure(figsize=(3.7,3.2))
colors = ['b','g','r']
Efit = []
for count, color in enumerate(colors):
    Efit.append(np.polyfit(GCNList,band_list[count], 1))
    plt.plot(np.sort(GCNList), np.poly1d(Efit[count])(np.sort(GCNList)), color + '--')
for count, color in enumerate(colors):
    plt.plot(GCNList, band_list[count], color + 'o')
plt.legend([r'${\epsilon}_{s}$=%.2fGCN - %.2f eV' %(Efit[0][0],abs(Efit[0][1]))
,r'${\epsilon}_{p}$=%.2fGCN + %.2f eV' %(Efit[1][0],Efit[1][1])
,r'${\epsilon}_{d}$=%.2fGCN - %.2f eV' %(Efit[2][0],abs(Efit[2][1]))]
,loc='best',frameon=False)
#plt.xlabel('Generalized coordination number (GCN)')
plt.ylabel('Band center (${\epsilon}$ - ${\epsilon}_{fermi}$) [eV]')
plt.show()

#plotting scaling of band center with GCN for surface sites
plt.figure(figsize=(3.7,3.2))
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
plt.legend([r'${\epsilon}_{s}$=%.2fGCN - %.2f eV' %(Efit[0][0],abs(Efit[0][1]))
,r'${\epsilon}_{p}$=%.2fGCN + %.2f eV' %(Efit[1][0],Efit[1][1])
,r'${\epsilon}_{d}$=%.2fGCN - %.2f eV' %(Efit[2][0],abs(Efit[2][1]))]
,loc='best',frameon=False)
#plt.xlabel('Generalized coordination number (GCN)')
plt.ylabel('Band center (${\epsilon}$ - ${\epsilon}_{fermi}$) [eV]')
plt.show()

#plotting scaling of band center with GCN for surface sites
plt.figure(figsize=(3.5,3.2))
colors = ['b','g','r']
Efit = []
for count, color in enumerate(colors):
    Efit.append(np.polyfit(GCNList[atom_type=='surface']\
                ,filling_list[count][atom_type=='surface'], 1))
    plt.plot(np.sort(GCNList[atom_type=='surface'])\
             , np.poly1d(Efit[count])\
             (np.sort(GCNList[atom_type=='surface'])), color + '--')
for count, color in enumerate(colors):
    plt.plot(GCNList[atom_type=='surface'], filling_list[count][atom_type=='surface'], color + 'o')
#plt.legend([r'${filling}_{s}$=%.2fGCN + %.2f states' %(Efit[0][0],Efit[0][1])
#,r'${filling}_{p}$=%.2fGCN + %.2f states' %(Efit[1][0],Efit[1][1])
#,r'${filling}_{d}$=%.2fGCN + %.2f states' %(Efit[2][0],Efit[2][1])]
#,loc='best',frameon=False)
plt.xlabel('Generalized coordination number (GCN)')
plt.ylabel('Filling [states]')
plt.show()

#plotting scaling of band center with GCN for surface sites
plt.figure(figsize=(3.5,3.2))
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
#plt.legend([r'${\epsilon}_{s}^{*}$=%.2fGCN + %.2f eV' %(Efit[0][0],Efit[0][1])
#,r'${\epsilon}_{p}^{*}$=%.2fGCN + %.2f eV' %(Efit[1][0],Efit[1][1])
#,r'${\epsilon}_{d}^{*}$=%.2fGCN + %.2f eV' %(Efit[2][0],Efit[2][1])]
#,loc='best',frameon=False)
#plt.xlabel('Generalized coordination number (GCN)')
plt.ylabel('Occupied band center (${\epsilon}^{*}$ - ${\epsilon}_{fermi}$) [eV]')
plt.show()

#plotting scaling of band center with GCN for surface sites
plt.figure(figsize=(3.5,3.2))
colors = ['b','g','r']
Efit = []
for count, color in enumerate(colors):
    Efit.append(np.polyfit(GCNList[atom_type=='surface']\
                ,unoccupied_band_list[count][atom_type=='surface'], 1))
    plt.plot(np.sort(GCNList[atom_type=='surface'])\
             , np.poly1d(Efit[count])\
             (np.sort(GCNList[atom_type=='surface'])), color + '--')
for count, color in enumerate(colors):
    plt.plot(GCNList[atom_type=='surface'], unoccupied_band_list[count][atom_type=='surface'], color + 'o')
#plt.legend([r'${\epsilon}_{s}^{un}$=%.2fGCN + %.2f eV' %(Efit[0][0],Efit[0][1])
#,r'${\epsilon}_{p}^{un}$=%.2fGCN + %.2f eV' %(Efit[1][0],Efit[1][1])
#,r'${\epsilon}_{d}^{un}$=%.2fGCN + %.2f eV' %(Efit[2][0],Efit[2][1])]
#,loc='best',frameon=False)
#plt.xlabel('Generalized coordination number (GCN)')
plt.ylabel('Unoccupied band center (${\epsilon}^{un}$ - ${\epsilon}_{fermi}$) [eV]')
plt.show()

#plotting scaling of occupied band center with GCN for surface sites
plt.figure(figsize=(3.5,3.2))
colors = ['b','g','r']
Efit = []
for count, color in enumerate(colors):
    Efit.append(np.polyfit(GCNList[atom_type=='surface']\
                ,second_moment_list[count][atom_type=='surface'], 1))
    plt.plot(np.sort(GCNList[atom_type=='surface'])\
             , np.poly1d(Efit[count])\
             (np.sort(GCNList[atom_type=='surface'])), color + '--')
for count, color in enumerate(colors):
    plt.plot(GCNList[atom_type=='surface'], second_moment_list[count][atom_type=='surface'], color + 'o')
#plt.legend([r'${\sigma^{2}}_{s}$=%.2fGCN + %.2f eV$^{2}$' %(Efit[0][0],Efit[0][1])
#,r'${\sigma^{2}}_{p}$=%.2fGCN + %.2f eV$^{2}$' %(Efit[1][0],Efit[1][1])
#,r'${\sigma^{2}}_{d}$=%.2fGCN + %.2f eV$^{2}$' %(Efit[2][0],Efit[2][1])]
#,loc='best',frameon=False)
#plt.xlabel('Generalized coordination number (GCN)')
plt.ylabel('Second moment [eV$^{2}$]')
plt.show()

#bond energy
plt.figure(figsize=(3.5,3.2))
colors = ['b','g','r']
Efit = []
for count, color in enumerate(colors):
    Efit.append(np.polyfit(GCNList[atom_type=='surface']\
                ,bond_energy_list[count][atom_type=='surface'], 1))
    plt.plot(np.sort(GCNList[atom_type=='surface'])\
             , np.poly1d(Efit[count])\
             (np.sort(GCNList[atom_type=='surface'])), color + '--')
for count, color in enumerate(colors):
    plt.plot(GCNList[atom_type=='surface'], bond_energy_list[count][atom_type=='surface'], color + 'o')
#plt.legend([r'${be}_{s}$=%.2fGCN + %.2f eV' %(Efit[0][0],Efit[0][1])
#,r'${be}_{p}$=%.2fGCN + %.2f eV' %(Efit[1][0],Efit[1][1])
#,r'${be}_{d}$=%.2fGCN + %.2f eV' %(Efit[2][0],Efit[2][1])]
#,loc='best',frameon=False)
#plt.xlabel('Generalized coordination number (GCN)')
plt.ylabel('Bond energy [eV]')
plt.show()