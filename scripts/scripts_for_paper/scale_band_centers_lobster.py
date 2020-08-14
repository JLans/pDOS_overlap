#!python2
# -*- coding: utf-8 -*-

"""
Created on Wed Jan 11 14:33:54 2017

@author: lansford
"""
from __future__ import division
import os
from pdos_overlap.coordination import get_geometric_data
import numpy as np
import matplotlib.pyplot as plt
from pdos_overlap.vasp_dos import VASP_DOS
from pdos_overlap.vasp_dos import get_all_VASP_files
from pdos_overlap.plotting_tools import set_figure_settings
from scipy.stats import linregress
set_figure_settings('paper')
Downloads_folder = os.path.join(os.path.expanduser("~"),'Downloads')
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
        r'C:\Users\lansf\Documents\Data\PROBE_PDOS\lobster_files_(N+1)bands\nanoparticles_noW')
DOSCAR_files = [DOSCAR + '.lobster' for DOSCAR in DOSCAR_files]

for DOSCAR, CONTCAR in zip(DOSCAR_files, CONTCAR_files):
    indices, GCNs, atom_types = get_geometric_data(CONTCAR)
    GCNList += GCNs.tolist()
    atom_type += atom_types.tolist()
    # read and return densityofstates object
    PDOS = VASP_DOS(DOSCAR)   
    for atom_index in indices:
        band_center = PDOS.get_band_center(atom_index, ['s','d']\
                                , sum_density=True) - PDOS.e_fermi
        occupied_band_center = PDOS.get_band_center(atom_index, ['s','d']\
                                , sum_density=True, max_energy=PDOS.e_fermi) - PDOS.e_fermi
        unoccupied_band_center = PDOS.get_band_center(atom_index, ['s','d']\
                                , sum_density=True, min_energy=PDOS.e_fermi) - PDOS.e_fermi
        band_width = PDOS.get_center_width(PDOS.e_fermi, atom_index, ['s','d']\
                                , sum_density=True)
        second_moment = PDOS.get_second_moment(atom_index, ['s','d']\
                                , sum_density=True)
        
        bond_energy = PDOS.get_bond_energy(atom_index, ['s','d']\
                                , sum_density=True)
            
        filling = PDOS.get_filling(atom_index, ['s','d']\
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
#plotting scaling of band center with GCN for surface sites
colors = ['b', 'r']
plt.figure(figsize=(3.5,3.2))
Efit = []
for count, color in enumerate(colors):
    Efit.append(np.polyfit(GCNList[atom_type=='surface']\
                ,filling_list[count][atom_type=='surface'], 1))
    plt.plot(np.sort(GCNList[atom_type=='surface'])\
             , np.poly1d(Efit[count])\
             (np.sort(GCNList[atom_type=='surface'])), color + '--')
for count, color in enumerate(colors):
    plt.plot(GCNList[atom_type=='surface'], filling_list[count][atom_type=='surface'], color + 'o')
plt.legend([r'${filling}_{s}$=%.2fGCN + %.2f states' %(Efit[0][0],Efit[0][1])
,r'${filling}_{d}$=%.2fGCN + %.2f states' %(Efit[1][0],Efit[1][1])]
,loc='best',frameon=False)
plt.xlabel('Generalized coordination number (GCN)')
plt.ylabel('Filling [states]')
plt.show()
#plotting scaling of band center with GCN
fig = plt.figure(figsize=(7.2,5),dpi=400)
axes = fig.subplots(nrows=2, ncols=2)
#plotting function
Efit = []
for count, color in enumerate(colors):
    slope, intercept, r_value, p_value, std_err = linregress(GCNList, band_list[count])
    Efit.append([slope, intercept])
    print('band center R^2 value and std_err')
    print(r_value**2)
    print(std_err)
    axes[0,0].plot(np.sort(GCNList), np.poly1d(Efit[count])(np.sort(GCNList)), color + '--')
for count, color in enumerate(colors):
    axes[0,0].plot(GCNList, band_list[count], color + 'o')
axes[0,0].legend([r'${\epsilon}_{s}$=%.2fGCN + %.2f eV' %(Efit[0][0],Efit[0][1])
,r'${\epsilon}_{d}$=%.2fGCN + %.2f eV' %(Efit[1][0],Efit[1][1])]
,loc=3,frameon=False)
#plt.xlabel('Generalized coordination number (GCN)')
axes[0,0].set_ylabel('Band center [eV]')
axes[0,0].text(0.01,0.92,'(a)',transform=axes[0,0].transAxes)
axes[0,0].set_ylim([-8, 2])
#plotting scaling of band center with GCN for surface sites
Efit = []
for count, color in enumerate(colors):
    slope, intercept, r_value, p_value, std_err = linregress(
        GCNList[atom_type=='surface'], band_list[count][atom_type=='surface'])
    Efit.append([slope, intercept])
    print('surface band center R^2 value and std_err')
    print(r_value**2)
    print(std_err)
    axes[0,1].plot(np.sort(GCNList[atom_type=='surface'])\
             , np.poly1d(Efit[count])\
             (np.sort(GCNList[atom_type=='surface'])), color + '--')
for count, color in enumerate(colors):
    axes[0,1].plot(GCNList[atom_type=='surface'], band_list[count][atom_type=='surface'], color + 'o')

axes[0,1].legend([r'${\epsilon}_{s}$=%.2fGCN + %.2f eV' %(Efit[0][0],Efit[0][1])
,r'${\epsilon}_{d}$=%.2fGCN + %.2f eV' %(Efit[1][0],Efit[1][1])]
,loc=3,frameon=False)
#plt.xlabel('Generalized coordination number (GCN)')
axes[0,1].set_ylabel('Band center [eV]')
axes[0,1].text(0.01,0.92,'(b)',transform=axes[0,1].transAxes)
axes[0,1].set_ylim([-4, 1])

Efit = []
for count, color in enumerate(colors):
    slope, intercept, r_value, p_value, std_err = linregress(
        GCNList[atom_type=='surface'], occupied_band_list[count][atom_type=='surface'])
    Efit.append([slope, intercept])
    print('occupied band center R^2 value and std_err')
    print(r_value**2)
    print(std_err)
    axes[1,0].plot(np.sort(GCNList[atom_type=='surface'])\
             , np.poly1d(Efit[count])\
             (np.sort(GCNList[atom_type=='surface'])), color + '--')
for count, color in enumerate(colors):
    axes[1,0].plot(GCNList[atom_type=='surface'], occupied_band_list[count][atom_type=='surface'], color + 'o')
axes[1,0].legend([r'${\epsilon}_{s}^{*}$=%.2fGCN + %.2f eV' %(Efit[0][0],Efit[0][1])
,r'${\epsilon}_{d}^{*}$=%.2fGCN + %.2f eV' %(Efit[1][0],Efit[1][1])]
,loc=3,frameon=False)
#plt.xlabel('Generalized coordination number (GCN)')
axes[1,0].set_ylabel('Occupied band center [eV]')
axes[1,0].text(0.01,0.92,'(c)',transform=axes[1,0].transAxes)
axes[1,0].set_ylim([-6.5, -1])
#plotting scaling of band center with GCN for surface sites
Efit = []
for count, color in enumerate(colors):
    slope, intercept, r_value, p_value, std_err = linregress(
        GCNList[atom_type=='surface'], unoccupied_band_list[count][atom_type=='surface'])
    Efit.append([slope, intercept])
    print('unoccupied band center R^2 value and std_err')
    print(r_value**2)
    print(std_err)
    axes[1,1].plot(np.sort(GCNList[atom_type=='surface'])\
             , np.poly1d(Efit[count])\
             (np.sort(GCNList[atom_type=='surface'])), color + '--')
for count, color in enumerate(colors):
    axes[1,1].plot(GCNList[atom_type=='surface'], unoccupied_band_list[count][atom_type=='surface'], color + 'o')
axes[1,1].legend([r'${\epsilon}_{s}^{un}$=%.2fGCN + %.2f eV' %(Efit[0][0],Efit[0][1])
,r'${\epsilon}_{d}^{un}$=%.2fGCN + %.2f eV' %(Efit[1][0],Efit[1][1])]
,loc=4,frameon=False)
axes[1,1].set_ylabel('Unoccupied band center [eV]')
axes[1,1].text(0.01,0.92,'(d)',transform=axes[1,1].transAxes)
axes[1,1].set_ylim([-1, 3.5])
fig.text(0.5, 0.01, 'Generalized Coordination Number (GCN)', ha='center')
figure_path = os.path.join(Downloads_folder,'band_center_lobster.jpg')
fig.set_tight_layout({'pad':2,'w_pad':1,'h_pad':0.25})
plt.savefig(figure_path, format='jpg')
plt.close()

#plotting scaling of occupied band center with GCN for surface sites
fig = plt.figure(figsize=(3.5,3.5),dpi=400)
axes = fig.subplots(nrows=2, ncols=1)
Efit = []
for count, color in enumerate(colors):
    slope, intercept, r_value, p_value, std_err = linregress(
        GCNList[atom_type=='surface'], (second_moment_list[count][atom_type=='surface'])**0.5)
    Efit.append([slope, intercept])
    print('second moment center R^2 value and std_err')
    print(r_value**2)
    print(std_err)
    axes[0].plot(np.sort(GCNList[atom_type=='surface'])\
             , np.poly1d(Efit[count])\
             (np.sort(GCNList[atom_type=='surface'])), color + '--')
for count, color in enumerate(colors):
    axes[0].plot(GCNList[atom_type=='surface'], (second_moment_list[count][atom_type=='surface'])**0.5, color + 'o')
axes[0].legend([r'${\sigma}_{s}$=%.2fGCN + %.2f eV' %(Efit[0][0],Efit[0][1])
,r'${\sigma}_{d}$=%.2fGCN + %.2f eV' %(Efit[1][0],Efit[1][1])]
,loc=4,frameon=False)
#plt.xlabel('Generalized coordination number (GCN)')
axes[0].set_ylabel('Band width [eV]')
axes[0].text(0.01,0.9,'(a)',transform=axes[0].transAxes)
axes[0].set_ylim([0, 5.05])
axes[0].set_xticks([])

#bond energy
Efit = []
for count, color in enumerate(colors):
    slope, intercept, r_value, p_value, std_err = linregress(
        GCNList[atom_type=='surface'], bond_energy_list[count][atom_type=='surface'])
    Efit.append([slope, intercept])
    print('bond energy R^2 value and std_err')
    print(r_value**2)
    print(std_err)
    axes[1].plot(np.sort(GCNList[atom_type=='surface'])\
             , np.poly1d(Efit[count])\
             (np.sort(GCNList[atom_type=='surface'])), color + '--')
for count, color in enumerate(colors):
    axes[1].plot(GCNList[atom_type=='surface'], bond_energy_list[count][atom_type=='surface'], color + 'o')
axes[1].legend([r'${be}_{s}$=%.2fGCN + %.2f eV' %(Efit[0][0],Efit[0][1])
,r'${be}_{d}$=%.2fGCN + %.2f eV' %(Efit[1][0],Efit[1][1])]
,loc=3,frameon=False)
#plt.xlabel('Generalized coordination number (GCN)')
axes[1].set_ylabel('Bond energy [eV]')
axes[1].text(0.01,0.9,'(b)',transform=axes[1].transAxes)
axes[1].set_ylim([-10, -2])
fig.text(0.55, 0.01, 'Generalized Coordination Number (GCN)', ha='center')
figure_path = os.path.join(Downloads_folder,'band_width_lobster.jpg')
fig.set_tight_layout({'pad':1.5,'w_pad':1,'h_pad':0.25})
plt.savefig(figure_path, format='jpg')
plt.close()