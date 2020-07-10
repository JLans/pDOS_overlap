# -*- coding: utf-8 -*-
"""
Created on Wed Oct 23 10:17:05 2019

@author: lansf
"""
import os
from jl_spectra_2_structure.plotting_tools import set_figure_settings
from jl_spectra_2_structure.cross_validation import LOAD_CROSS_VALIDATION
set_figure_settings('paper')

#Folder where figures are to be saved
Downloads_folder = os.path.join(os.path.expanduser("~"),'Downloads')
CV_path = (r'C:\Users\lansf\Documents\Data\PROBE_PDOS\NN_files\cv_BW'
'/learning_curves\include_low_frequencies\cross_validation_C2H4_GCN_low')
C2H4_CV = LOAD_CROSS_VALIDATION(cross_validation_path=CV_path)
C2H4_CV.load_CV_class(0)
C2H4_CV.MAINconv.set_GCNlabels(Minimum=C2H4_CV.MIN_GCN_PER_LABEL
                               , showfigures=True
                               , figure_directory=Downloads_folder
                               , BINDING_TYPE_FOR_GCN=[1,2])

