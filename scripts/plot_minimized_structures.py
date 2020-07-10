# -*- coding: utf-8 -*-
"""
Created on Tue Apr 25 11:46:31 2017

@author: lansford
"""

from __future__ import absolute_import, division, print_function
import os
import json
import numpy as np
import matplotlib.pyplot as plt
from pdos_overlap.plotting_tools import set_figure_settings
from pdos_overlap.plotting_tools import adjust_text
set_figure_settings('paper')
Downloads_folder = os.path.join(os.path.expanduser("~"),'Downloads')
#from jl_spectra_2_structure.primary_data_creation.dft_2_data import Primary_DATA
#CO_path = (r'C:\Users\lansf\Documents\Data\IR_Materials_Gap\Lansford'
#           ' - Nature Communications 2020\VASP files\CO_nanoparticle\spin')
#CO_data = Primary_DATA(create_new_vasp_files=False)
#CO_data.generate_primary_data(CO_path,Downloads_folder+'minimized_spin.json')
CO_path = r'C:\Users\lansf\Documents\Data\PROBE_PDOS\minimized_spin.json'
with open(CO_path, 'r') as infile:
    CO_dict = json.load(infile)
dict_filter = np.all((np.array(CO_dict['CN_ADSORBATE'])==1
                      , np.min(CO_dict['FREQUENCIES'],axis=1) > 0
                      , np.array(CO_dict['NUM_METAL']) > 1),axis=0)
frequencies = np.array(CO_dict['FREQUENCIES'])[dict_filter]
GCN = np.array(CO_dict['GCN'])[dict_filter]
NumPt = np.array(CO_dict['NUM_METAL'])[dict_filter]
freqfit = np.polyfit(GCN,frequencies.max(axis=-1), 1)
plt.figure(figsize=(7.2,4), dpi=400)
plt.plot(np.sort(GCN), np.poly1d(freqfit)
             (np.sort(GCN)), 'k--')
plt.legend([r'${\nu}_{CO}$=%.2fGCN + %.2f cm$^{-1}$' %(freqfit[0],freqfit[1])]
,loc='best',frameon=False)
plt.plot(GCN,frequencies.max(axis=-1),'o')
texts = []
for x, y, s in zip(GCN, frequencies.max(axis=-1), NumPt.astype(int)):
    texts.append(plt.text(x, y, s, bbox={'pad':0, 'alpha':0}))
adjust_text(texts,autoalign=True,only_move={'points':'y','text':'y'}, arrowprops=dict(arrowstyle="-", color='k', lw=2))
plt.xlabel('GCN')
plt.ylabel('C-O frequency [cm$^{-1}$]')
plt.savefig(os.path.join(Downloads_folder,'GCN_scale.jpg'), format='jpg')
plt.close()
