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


gas = 'NO'
adsorbate='NO_fcc'
surface = 'Pt111'
set_figure_settings('paper')
np.set_printoptions(linewidth=100)
example_path = r'C:\Users\lansf\Documents\Data\PROBE_PDOS\vasp_dos_files'
lobster_path = r'C:\Users\lansf\Documents\Data\PROBE_PDOS\lobster_files'
GAS_DOSCAR = os.path.join(lobster_path, gas + '/DOSCAR.lobster')
GAS_CONTCAR = os.path.join(lobster_path, gas + '/CONTCAR')
ADSORBATE_DOSCAR = os.path.join(lobster_path, 'gas+Pt_G.03_noW/'+surface + '+'\
                          + adsorbate + '/DOSCAR.lobster')
ADSORBATE_CONTCAR = os.path.join(lobster_path, 'gas+Ptnano/'+surface + '+'\
                          + adsorbate + '/CONTCAR')
BULK_DOSCAR = os.path.join(example_path,'Pt_nano/Pt147/DOSCAR')
# VASP_DOS objects for both the gas (vacuum) and the adsorbate+surface system
GAS_PDOS = VASP_DOS(GAS_DOSCAR)
REFERENCE_PDOS = VASP_DOS(ADSORBATE_DOSCAR)
BULK_PDOS = VASP_DOS(BULK_DOSCAR)

# Get adsorbate and site indices and initialize PDOS_OVERLAP object
adsorbate_indices, site_indices = get_adsorbate_indices(GAS_CONTCAR\
                                                        , ADSORBATE_CONTCAR)
#Initialize Coordination object. Repeat is necessary so it doesn't count itself
NO_overlap = PDOS_OVERLAP(GAS_PDOS, REFERENCE_PDOS, adsorbate_indices\
                          , site_indices, min_occupation=1\
                          , upshift=0.5, energy_weight=3)
#NO_overlap.optimize_energy_shift(bound=[-0.5,1.5], reset=True)

GCNList = []
orbital_1_list = []
orbital_2_list = []
orbital_3_list = []
orbital_4_list = []

DOSCAR_files, CONTCAR_files = get_all_VASP_files(\
        r'C:\Users\lansf\Documents\Data\PROBE_PDOS\vasp_dos_files\Pt_nano')


for nano_DOSCAR, nano_CONTCAR in zip(DOSCAR_files, CONTCAR_files):
    nano_indices, GCNs, atom_types = get_geometric_data(nano_CONTCAR)
    GCNList += GCNs[atom_types[...] == 'surface'].tolist()
    # read and return density of states object
    nano_PDOS = VASP_DOS(nano_DOSCAR)   
    for atom_index in nano_indices[atom_types[...] == 'surface']:
        orbital_1 = NO_overlap.get_orbital_interaction(1\
                    , nano_PDOS, atom_index, ['s','d']\
                    , BULK_PDOS, bulk_atom=43\
                    , method='band_width', use_orbital_proximity=False\
                    , index_type='adsorbate'
                    , sum_interaction=True)
        orbital_2 = NO_overlap.get_orbital_interaction(2\
                    , nano_PDOS, atom_index, ['s','d']\
                    , BULK_PDOS, bulk_atom=43\
                    , method='band_width', use_orbital_proximity=False\
                    , index_type='adsorbate'
                    , sum_interaction=True)
        orbital_3 = NO_overlap.get_orbital_interaction(3\
                    , nano_PDOS, atom_index, ['s','d']\
                    , BULK_PDOS, bulk_atom=43\
                    , method='band_width', use_orbital_proximity=False\
                    , index_type='adsorbate'
                    , sum_interaction=True)
        orbital_4 = NO_overlap.get_orbital_interaction(4\
                    , nano_PDOS, atom_index, ['s','d']\
                    , BULK_PDOS, bulk_atom=43\
                    , method='band_width', use_orbital_proximity=False\
                    , index_type='adsorbate'
                    , sum_interaction=True)
        
            
        orbital_1_list.append(orbital_1)
        orbital_2_list.append(orbital_2)
        orbital_3_list.append(orbital_3)
        orbital_4_list.append(orbital_4)

GCNList = np.array(GCNList)
orbital_1_list = np.array(orbital_1_list).T
orbital_2_list = np.array(orbital_2_list).T
orbital_3_list = np.array(orbital_3_list).T
orbital_4_list = np.array(orbital_4_list).T

#plotting scaling of band center with GCN
for count, orbital_list in enumerate([orbital_1_list, orbital_2_list\
                                     , orbital_3_list]):
    plt.figure(figsize=(7,5))
    colors = ['b','r']
    orbitalfit = []
    for count, color in enumerate(colors):
        orbitalfit.append(np.polyfit(GCNList,orbital_list[count], 1))
        plt.plot(np.sort(GCNList), np.poly1d(orbitalfit[count])(np.sort(GCNList)), color + '--')
    for count, color in enumerate(colors):
        plt.plot(GCNList, orbital_list[count], color + 'o')
    plt.legend([r'${interaction}_{s}$=%.2fGCN + %.2f eV' %(orbitalfit[0][0], orbitalfit[0][1])
    ,r'${interaction}_{d}$=%.2fGCN + %.2f eV' %(orbitalfit[1][0], orbitalfit[1][1])]
    ,loc='best',frameon=False)
    plt.xlabel('Generalized coordination number (GCN)')
    plt.ylabel('Interaction energy [eV]')
    plt.show()

