# -*- coding: utf-8 -*-
"""
Created on Tue Apr 25 11:46:31 2017

@author: lansford
"""

from __future__ import absolute_import, division, print_function
import os
from jl_spectra_2_structure.cross_validation import LOAD_CROSS_VALIDATION
from jl_spectra_2_structure.plotting_tools import set_figure_settings

Downloads_folder = os.path.join(os.path.expanduser("~"),'Downloads')
cv_paper = (r'C:\Users\lansf\Documents\Data\PROBE_PDOS\NN_files\cv_BW'
'/learning_curves/data_for_paper')
set_figure_settings('paper')
CV_class = LOAD_CROSS_VALIDATION(cross_validation_path=cv_paper)
CV_class.load_all_CV_data()
print(CV_class.NUM_TRAIN)
print(CV_class.TRAINING_ERROR)
print(CV_class.NN_PROPERTIES)

        
CV_class.plot_models(CV_class.CV_RESULTS,figure_directory=Downloads_folder
                     ,model_list=[0,1,2,3]
                     ,xlim=[0,1200,], ylim1=[0.008,0.082], ylim2=[0.025,0.2])
CV_class.plot_parity_plots(figure_directory=Downloads_folder,model_list=[0,1,2,3])