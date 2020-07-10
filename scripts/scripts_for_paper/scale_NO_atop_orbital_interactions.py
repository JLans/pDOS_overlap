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
Downloads_folder = os.path.join(os.path.expanduser("~"),'Downloads')

gas = 'NO'
adsorbate='NO_atop'
surface = 'Pt111'
set_figure_settings('paper')
np.set_printoptions(linewidth=100)
example_path = r'C:\Users\lansf\Documents\Data\PROBE_PDOS\vasp_dos_files'
lobster_path = r'C:\Users\lansf\Documents\Data\PROBE_PDOS\lobster_files'
GAS_DOSCAR = os.path.join(lobster_path, gas + '/DOSCAR.lobster')
GAS_CONTCAR = os.path.join(lobster_path, gas + '/CONTCAR')
ADSORBATE_DOSCAR = os.path.join(lobster_path, 'surfaces_noW/'+surface + '+'\
                          + adsorbate + '/DOSCAR.lobster')
ADSORBATE_CONTCAR = os.path.join(lobster_path, 'surfaces_noW/'+surface + '+'\
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
Overlap = PDOS_OVERLAP(GAS_PDOS, REFERENCE_PDOS, adsorbate_indices\
                          , site_indices, min_occupation=0.9\
                          , upshift=0.5, energy_weight=4)
Overlap.optimize_energy_shift(bound=[-10, 10], reset=True)

GCNList = []
orbital_indices = [0, 1, 2, 3]
orbital_array = [[] for i in range(len(orbital_indices))]
DOSCAR_files, CONTCAR_files = get_all_VASP_files(\
        r'C:\Users\lansf\Documents\Data\PROBE_PDOS\vasp_dos_files\Pt_nano')

for nano_DOSCAR, nano_CONTCAR in zip(DOSCAR_files, CONTCAR_files):
    nano_indices, GCNs, atom_types = get_geometric_data(nano_CONTCAR)
    GCNList += GCNs[atom_types[...] == 'surface'].tolist()
    # read and return density of states object
    nano_PDOS = VASP_DOS(nano_DOSCAR)   
    for atom_index in nano_indices[atom_types[...] == 'surface']:
        for i, index in enumerate(orbital_indices):
            orbital_array[i].append(Overlap.get_orbital_interaction(index\
                    , nano_PDOS, atom_index, ['s','d']\
                    , BULK_PDOS, bulk_atom=43\
                    , method='band_width', use_orbital_proximity=False\
                    , index_type='adsorbate'
                    , sum_interaction=True))
            


GCNList = np.array(GCNList)
for i in range(len(orbital_indices)):
    orbital_array[i] = np.array(orbital_array[i]).T
#plotting scaling of band center with GCN
plt.figure(figsize=(7,5))
colors = ['b','r']
abc = ['(a)','(b)','(c)','(d)']
fig = plt.figure(figsize=(7.2,5),dpi=400)
axes = fig.subplots(nrows=2, ncols=2)
axis_list = [axes[0, 0], axes[0, 1], axes[1, 0], axes[1, 1]]
for count, orbital_list in enumerate(orbital_array):
    orbitalfit = []
    for count2, color in enumerate(colors):
        orbitalfit.append(np.polyfit(GCNList,orbital_list[count2], 1))
        axis_list[count].plot(np.sort(GCNList), np.poly1d(orbitalfit[count2])(np.sort(GCNList)), color + '--')
    for count2, color in enumerate(colors):
        axis_list[count].plot(GCNList, orbital_list[count2], color + 'o')
    axis_list[count].legend([r'${H}_{s, %d}$ = %.2fGCN + %.2f eV' %(orbital_indices[count], orbitalfit[0][0], orbitalfit[0][1])
,r'${H}_{d, %d}$ = %.2fGCN + %.2f eV' %(orbital_indices[count], orbitalfit[1][0], orbitalfit[1][1])]
    ,loc='best',frameon=False, handlelength=1)
    axis_list[count].text(-0.05,1.05,abc[count],transform=axis_list[count].transAxes)
fig.text(0.001, 0.5, 'Interaction energy [eV]', va='center', rotation='vertical')
fig.text(0.5, 0.01, 'Generalized coordination number (GCN)', ha='center')
figure_path = os.path.join(Downloads_folder,'NOatop_interaction.jpg')
fig.set_tight_layout({'pad':2,'w_pad':0.5,'h_pad':0.5})
plt.savefig(figure_path, format='jpg')
plt.close()
