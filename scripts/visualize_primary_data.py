# -*- coding: utf-8 -*-
"""
Created on Tue Apr 25 11:46:31 2017

@author: lansford
"""

from __future__ import division
import os
import json_tricks
import numpy as np
import matplotlib.pyplot as plt
from matplotlib import rcParams
from jl_spectra_2_structure import get_default_data_paths
from jl_spectra_2_structure.plotting_tools import set_figure_settings
from sklearn.decomposition import PCA
from sklearn.preprocessing import StandardScaler
from sklearn.pipeline import Pipeline
set_figure_settings('paper',**{'lines.markersize': 3})

#Folder where figures are to be saved
Downloads_folder = os.path.join(os.path.expanduser("~"),'Downloads')

#Generate Primary Set of Data from DFT for CO

#Only change this code below if you wish to use a different set of adsorbates
#and metal atoms

Adsorbates_list = ['C2H4']
pipeline = Pipeline([('scaling', StandardScaler()), ('pca', PCA(n_components=2))])
Xlim_list = [[-5,17]]
xlabel = ['Frequency PC-1 [cm$^{-1}$]']
ylabel1 = ['Frequency PC-2 [cm$^{-1}$]']
ylabel2 = ['Intensity PC-1']
marker = ['o','s','^','D']
CN_POC1 = [1,2,3,4]
color = ['green','red','orange','blue']
params = {'figure.autolayout': False,
              'axes.labelpad':2}
rcParams.update(params)
for count, adsorbate in enumerate(Adsorbates_list):

    nanoparticle_path, isotope_path, high_coverage_path\
    , coverage_scaling_path = get_default_data_paths(adsorbate)
    
    with open(nanoparticle_path, 'r') as infile:
        SingleCO= json_tricks.load(infile, preserve_order=False)
    
    BINDING_TYPES_unfiltered = np.array(SingleCO['CN_ADSORBATE'])
    cluster = np.array(SingleCO['CN_METAL']) >0
    max_forces = np.array(SingleCO['MAX_FORCE'])
    is_adsorbed = BINDING_TYPES_unfiltered >0
    correct_coordination = BINDING_TYPES_unfiltered <= 4
    small_force = max_forces < 0.05
    select_files = np.all((is_adsorbed,small_force,correct_coordination,cluster),axis=0)
    Xfreq_ALL = np.array(SingleCO['FREQUENCIES'])[select_files]
    Xint_ALL = np.array(SingleCO['INTENSITIES'])[select_files]
    is_local_minima = np.min(Xfreq_ALL,axis=1)>0
    CNCOList = np.array(SingleCO['CN_ADSORBATE'])[select_files]
    GCNList = np.array(SingleCO['GCN'])[select_files]
    if adsorbate == 'C2H4':
        Xfreq_PCA = pipeline.fit_transform(np.array(Xfreq_ALL))
        Xint_PCA = pipeline.fit_transform(np.array(Xint_ALL))
        COfrequencies = Xfreq_PCA[:,0]
        COintensities = Xint_PCA[:,0]
        MCfrequencies = Xfreq_PCA[:,1]
        MCintensities = Xint_PCA[:,1]
        site_filter = np.any((CNCOList==1,CNCOList==2),axis=0)
    else:
        COfrequencies = np.max(Xfreq_ALL,axis=1)
        COintensities = np.max(Xint_ALL,axis=1)
        MCfrequencies = np.max(Xfreq_ALL[:,0:5],axis=1)
        MCintensities = Xint_ALL[np.arange(len(Xint_ALL)),np.argmax(Xfreq_ALL[:,0:5],axis=1)]
        site_filter = CNCOList==1    
    fig = plt.figure(0,figsize=(7.2,4),dpi=400)
    axes = fig.subplots(nrows=2, ncols=2)
    for i in range(4):
        axes[0][0].plot(COfrequencies[CNCOList==CN_POC1[i]],GCNList[CNCOList==CN_POC1[i]],marker=marker[i],color=color[i],linestyle='None')
    axes[0][0].set_xticks([])
    axes[0][0].set_xlim(Xlim_list[count])
    axes[0][0].set_ylabel('Generalized\n Coordination Number')
    axes[0][0].text(0.01,0.93,'(a)',name ='Calibri',transform=axes[0][0].transAxes)
    
    for i in range(4):
        axes[1][0].plot(COfrequencies[np.all((CNCOList==CN_POC1[i],is_local_minima==True),axis=0)]\
            ,GCNList[np.all((CNCOList==CN_POC1[i],is_local_minima==True),axis=0)],marker=marker[i],color=color[i],linestyle='None')
    axes[1][0].set_xlabel(xlabel[count])
    axes[1][0].set_xlim(Xlim_list[count])
    axes[1][0].set_ylabel('Generalized\n Coordination Number')
    axes[1][0].text(0.01,0.93,'(b)',name ='Calibri',transform=axes[1][0].transAxes)
    
    im = axes[0][1].scatter(x=COfrequencies[np.all((site_filter,is_local_minima==True),axis=0)],\
                y=MCfrequencies[np.all((site_filter,is_local_minima==True),axis=0)],\
                c=GCNList[np.all((site_filter,is_local_minima==True),axis=0)])
    axes[0][1].set_xticks([])
    axes[0][1].set_ylabel(ylabel1[count])
    axes[0][1].text(0.01,0.93,'(c)',name ='Calibri',transform=axes[0][1].transAxes)
    
    im=axes[1][1].scatter(x=COfrequencies[np.all((site_filter,is_local_minima==True),axis=0)],\
                y=COintensities[np.all((site_filter,is_local_minima==True),axis=0)],\
                c=GCNList[np.all((site_filter,is_local_minima==True),axis=0)])
    if adsorbate != 'C2H4':
        axes[1][1].set_yscale('log')
    axes[1][1].set_xlabel(xlabel[count])
    axes[1][1].set_ylabel(ylabel2[count])
    axes[1][1].text(0.01,0.93,'(d)',name ='Calibri',transform=axes[1][1].transAxes)
    plt.gcf().subplots_adjust(bottom=0.09,top=0.98,left=0.07,right=0.98,wspace=0.2,hspace=0.05)
    cbar = fig.colorbar(im, ax=axes.ravel().tolist(),fraction=0.07,pad=0.01)
    cbar.set_label('Generalized Coordination Number')
    fig.savefig(os.path.join(Downloads_folder,'Removing Outliers_'+adsorbate+'.jpg'), format='jpg')
    plt.close()
rcParams.update({'figure.autolayout': True})