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

gas = 'CO'
adsorbate='CO'
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
CO_overlap = PDOS_OVERLAP(GAS_PDOS, REFERENCE_PDOS, adsorbate_indices\
                          , site_indices, min_occupation=1.5\
                          , upshift=0.5, energy_weight=4)
CO_overlap.optimize_energy_shift(bound=[-10, 10], reset=True)

GCNList = []
three_sigma_list = []
four_sigma_list = []
one_pi_list = []
five_sigma_list = []

DOSCAR_files, CONTCAR_files = get_all_VASP_files(\
        r'C:\Users\lansf\Documents\Data\PROBE_PDOS\vasp_dos_files\Pt_nano')
DOSCAR_files, CONTCAR_files = get_all_VASP_files(\
        r'C:\Users\lansf\Documents\Data\PROBE_PDOS\lobster_files\nanoparticles_noW')
DOSCAR_files = [file_name + '.lobster' for file_name in DOSCAR_files]


for nano_DOSCAR, nano_CONTCAR in zip(DOSCAR_files, CONTCAR_files):
    nano_indices, GCNs, atom_types = get_geometric_data(nano_CONTCAR)
    GCNList += GCNs[atom_types[...] == 'surface'].tolist()
    # read and return density of states object
    nano_PDOS = VASP_DOS(nano_DOSCAR,add_p2s=False)   
    for atom_index in nano_indices[atom_types[...] == 'surface']:
        three_sigma = CO_overlap.get_orbital_interaction(0\
                    , nano_PDOS, atom_index, ['s','dz2']\
                    , BULK_PDOS, bulk_atom=43\
                    , method='band_width', use_orbital_proximity=False
                    , index_type='adsorbate')      
        four_sigma = CO_overlap.get_orbital_interaction(2\
                    , nano_PDOS, atom_index, ['s','dz2']\
                    , BULK_PDOS, bulk_atom=43\
                    , method='band_width', use_orbital_proximity=False
                    , index_type='adsorbate')
        one_pi = CO_overlap.get_orbital_interaction(3\
                    , nano_PDOS, atom_index, ['dyz','dxz']\
                    , BULK_PDOS, bulk_atom=43\
                    , method='band_width', use_orbital_proximity=False
                    , index_type='adsorbate')
        five_sigma = CO_overlap.get_orbital_interaction(1\
                    , nano_PDOS, atom_index, ['s','dz2']\
                    , BULK_PDOS, bulk_atom=43\
                    , method='band_width', use_orbital_proximity=False
                    , index_type='adsorbate')
        
        three_sigma_list.append(three_sigma)
        four_sigma_list.append(four_sigma)
        one_pi_list.append(one_pi)
        five_sigma_list.append(five_sigma)

GCNList = np.array(GCNList)
three_sigma_list = np.array(three_sigma_list).T
four_sigma_list = np.array(four_sigma_list).T
one_pi_list = np.array(one_pi_list).T
five_sigma_list = np.array(five_sigma_list).T

fig = plt.figure(figsize=(7.2,5),dpi=400)
axes = fig.subplots(nrows=2, ncols=2)
colors = ['b','r']
orbitalfit = []
for count, color in enumerate(colors):
    orbitalfit.append(np.polyfit(GCNList,three_sigma_list[count], 1))
    axes[0,0].plot(np.sort(GCNList), np.poly1d(orbitalfit[count])(np.sort(GCNList)), color + '--')
for count, color in enumerate(colors):
    axes[0,0].plot(GCNList, three_sigma_list[count], color + 'o')
axes[0,0].legend([r'${H}_{s,3\sigma}$ = %.2fGCN + %.2f eV' %(orbitalfit[0][0], orbitalfit[0][1])
,r'${H}_{dz^{2},3\sigma}$ = %.2fGCN + %.2f eV' %(orbitalfit[1][0], orbitalfit[1][1])]
,loc='best',frameon=False, handlelength=1)
axes[0,0].text(-0.05,1.05,'(a)',transform=axes[0,0].transAxes)

#plotting scaling of band center with GCN
orbitalfit = []
for count, color in enumerate(colors):
    orbitalfit.append(np.polyfit(GCNList,five_sigma_list[count], 1))
    axes[0,1].plot(np.sort(GCNList), np.poly1d(orbitalfit[count])(np.sort(GCNList)), color + '--')
for count, color in enumerate(colors):
    axes[0,1].plot(GCNList, five_sigma_list[count], color + 'o')
axes[0,1].legend([r'${H}_{s,5\sigma}$ = %.2fGCN + %.2f eV' %(orbitalfit[0][0], orbitalfit[0][1])
,r'${H}_{dz^{2},5\sigma}$ = %.2fGCN + %.2f eV' %(orbitalfit[1][0], orbitalfit[1][1])]
,loc='best',frameon=False, handlelength=1)
axes[0,1].text(-0.05,1.05,'(b)',transform=axes[0,1].transAxes)

orbitalfit = []
for count, color in enumerate(colors):
    orbitalfit.append(np.polyfit(GCNList,four_sigma_list[count], 1))
    axes[1,0].plot(np.sort(GCNList), np.poly1d(orbitalfit[count])(np.sort(GCNList)), color + '--')
for count, color in enumerate(colors):
    axes[1,0].plot(GCNList, four_sigma_list[count], color + 'o')
axes[1,0].legend([r'${H}_{s,4\sigma}$ = %.2fGCN + %.2f eV' %(orbitalfit[0][0], orbitalfit[0][1])
,r'${H}_{dz^{2},4\sigma}$ = %.2fGCN + %.2f eV' %(orbitalfit[1][0], orbitalfit[1][1])]
,loc='best',frameon=False, handlelength=1)
axes[1,0].text(-0.05,1.05,'(c)',transform=axes[1,0].transAxes)

#plotting scaling of band center with GCN
orbitalfit = []
for count, color in enumerate(colors):
    orbitalfit.append(np.polyfit(GCNList,one_pi_list[count], 1))
    axes[1,1].plot(np.sort(GCNList), np.poly1d(orbitalfit[count])(np.sort(GCNList)), color + '--')
for count, color in enumerate(colors):
    axes[1,1].plot(GCNList, one_pi_list[count], color + 'o')
axes[1,1].legend([r'${H}_{dyz,2\pi}$ = %.2fGCN + %.2f eV' %(orbitalfit[0][0], orbitalfit[0][1])
,r'${H}_{dxz,2\pi}$ = %.2fGCN + %.2f eV' %(orbitalfit[1][0], orbitalfit[1][1])]
,loc='best',frameon=False, handlelength=1)
axes[1,1].text(-0.05,1.05,'(d)',transform=axes[1,1].transAxes)
fig.text(0.001, 0.5, 'Interaction energy [eV]', va='center', rotation='vertical')
fig.text(0.5, 0.01, 'Generalized coordination number (GCN)', ha='center')
figure_path = os.path.join(Downloads_folder,'CO_interaction.jpg')
fig.set_tight_layout({'pad':2,'w_pad':0.5,'h_pad':0.5})
plt.savefig(figure_path, format='jpg')
plt.close()

