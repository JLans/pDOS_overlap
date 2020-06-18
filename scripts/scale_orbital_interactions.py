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
from pdos_overlap import set_figure_settings
from pdos_overlap import get_adsorbate_indices
from pdos_overlap import PDOS_OVERLAP


gas = 'CO'
surface = 'Pt111'
set_figure_settings('paper')
np.set_printoptions(linewidth=100)
example_path = r'C:\Users\lansf\Documents\Data\PROBE_PDOS\vasp_dos_files'
GAS_DOSCAR = os.path.join(example_path, gas + '/DOSCAR')
GAS_CONTCAR = os.path.join(example_path, gas + '/CONTCAR')
ADSORBATE_DOSCAR = os.path.join(example_path, gas + '+' + surface + '/DOSCAR')
ADSORBATE_CONTCAR = os.path.join(example_path, gas + '+' + surface + '/CONTCAR')
BULK_DOSCAR = os.path.join(example_path,'Pt_nano/Pt147/DOSCAR')
# VASP_DOS objects for both the gas (vacuum) and the adsorbate+surface system
GAS_PDOS = VASP_DOS(GAS_DOSCAR)
REFERENCE_PDOS = VASP_DOS(ADSORBATE_DOSCAR)
BULK_PDOS = VASP_DOS(BULK_DOSCAR)

# Get adsorbate and site indices and initialize PDOS_OVERLAP object
adsorbate_indices, site_indices = get_adsorbate_indices(GAS_CONTCAR\
                                                        , ADSORBATE_CONTCAR)
#Initialize Coordination object. Repeat is necessary so it doesn't count itself
CO_overlap = PDOS_OVERLAP(GAS_PDOS, REFERENCE_PDOS, adsorbate_indices\
                          , site_indices, min_occupation=0.9\
                          , upshift=0.5, energy_weight=3)
CO_overlap.optimize_energy_shift(bound=[-0.5,1.5], reset=True)

GCNList = []
four_sigma_list = []
one_pi_list = []
five_sigma_list = []

DOSCAR_files, CONTCAR_files = get_all_VASP_files(\
        r'C:\Users\lansf\Documents\Data\PROBE_PDOS\vasp_dos_files\Pt_nano')


for nano_DOSCAR, nano_CONTCAR in zip(DOSCAR_files, CONTCAR_files):
    nano_indices, GCNs, atom_types = get_geometric_data(nano_CONTCAR)
    GCNList += GCNs[atom_types[...] == 'surface'].tolist()
    # read and return density of states object
    nano_PDOS = VASP_DOS(nano_DOSCAR)   
    for atom_index in nano_indices[atom_types[...] == 'surface']:
        four_sigma = CO_overlap.calculate_orbital_interaction(1\
                    , nano_PDOS, atom_index, ['s','pz','dz2']\
                    , BULK_PDOS, bulk_atom=43\
                    , method='orbital_bond_energy', use_orbital_proximity=False)
        one_pi = CO_overlap.calculate_orbital_interaction(2\
                    , nano_PDOS, atom_index, ['dyz','dxz']\
                    , BULK_PDOS, bulk_atom=43\
                    , method='orbital_bond_energy', use_orbital_proximity=False)
        five_sigma = CO_overlap.calculate_orbital_interaction(3\
                    , nano_PDOS, atom_index, ['s','pz','dz2']\
                    , BULK_PDOS, bulk_atom=43\
                    , method='orbital_bond_energy', use_orbital_proximity=False)
        
            
        four_sigma_list.append(four_sigma)
        one_pi_list.append(one_pi)
        five_sigma_list.append(five_sigma)

GCNList = np.array(GCNList)
four_sigma_list = np.array(four_sigma_list).T
one_pi_list = np.array(one_pi_list).T
five_sigma_list = np.array(five_sigma_list).T

#plotting scaling of band center with GCN
plt.figure(0,figsize=(7,5))
colors = ['b','g','r']
orbitalfit = []
for count, color in enumerate(colors):
    orbitalfit.append(np.polyfit(GCNList,four_sigma_list[count], 1))
    plt.plot(np.sort(GCNList), np.poly1d(orbitalfit[count])(np.sort(GCNList)), color + '--')
for count, color in enumerate(colors):
    plt.plot(GCNList, four_sigma_list[count], color + 'o')
plt.legend([r'${interaction}_{s}$=%.2fGCN + %.2f eV' %(orbitalfit[0][0], orbitalfit[0][1])
,r'${interaction}_{pz}$=%.2fGCN + %.2f eV' %(orbitalfit[1][0], orbitalfit[1][1])
,r'${interaction}_{dz2}$=%.2fGCN + %.2f eV' %(orbitalfit[2][0], orbitalfit[2][1])]
,loc='best',frameon=False)
plt.xlabel('Generalized coordination number (GCN)')
plt.ylabel('Interaction energy [states]')
plt.show()

#plotting scaling of band center with GCN
plt.figure(1,figsize=(7,5))
colors = ['b','g','r']
orbitalfit = []
for count, color in enumerate(colors):
    orbitalfit.append(np.polyfit(GCNList,five_sigma_list[count], 1))
    plt.plot(np.sort(GCNList), np.poly1d(orbitalfit[count])(np.sort(GCNList)), color + '--')
for count, color in enumerate(colors):
    plt.plot(GCNList, five_sigma_list[count], color + 'o')
plt.legend([r'${interaction}_{s}$=%.2fGCN + %.2f eV' %(orbitalfit[0][0], orbitalfit[0][1])
,r'${interaction}_{pz}$=%.2fGCN + %.2f eV' %(orbitalfit[1][0], orbitalfit[1][1])
,r'${interaction}_{dz2}$=%.2fGCN + %.2f eV' %(orbitalfit[2][0], orbitalfit[2][1])]
,loc='best',frameon=False)
plt.xlabel('Generalized coordination number (GCN)')
plt.ylabel('Interaction energy [eV]')
plt.show()

#plotting scaling of band center with GCN
plt.figure(2,figsize=(7,5))
colors = ['r','r']
orbitalfit = []
for count, color in enumerate(colors):
    orbitalfit.append(np.polyfit(GCNList,one_pi_list[count], 1))
    plt.plot(np.sort(GCNList), np.poly1d(orbitalfit[count])(np.sort(GCNList)), color + '--')
for count, color in enumerate(colors):
    plt.plot(GCNList, one_pi_list[count], color + 'o')
plt.legend([r'${interaction}_{dyz}$=%.2fGCN + %.2f eV' %(orbitalfit[0][0], orbitalfit[0][1])
,r'${interaction}_{dxz}$=%.2fGCN + %.2f eV' %(orbitalfit[1][0], orbitalfit[1][1])]
,loc='best',frameon=False)
plt.xlabel('Generalized coordination number (GCN)')
plt.ylabel('Interaction energy [eV]')
plt.show()

