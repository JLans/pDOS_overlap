# -*- coding: utf-8 -*-
"""
Created on Tue May 19 22:28:09 2020

@author: lansf
"""

from __future__ import absolute_import, division, print_function
import os
import pkg_resources
import numpy as np

def get_data_path():
    """ Get default paths to experimental data.
    
    Returns
    -------
    data_path : str
        path to package data
    """
    data_path = pkg_resources.resource_filename(__name__, 'data/')
    return data_path

def get_example_data():
    """ Get default paths to experimental data.
    
    Returns
    -------
    example_data_path : str
        path to example VASP data
    """
    data_path = pkg_resources.resource_filename(__name__, 'data/')
    example_data_path = os.path.join(data_path, 'example_data')
    return example_data_path

class VASP_DOS:
    """Class for generating complex synthetic IR spectra"""
    def __init__(self):
        """ 
        Parameters
        ----------
        Blank : str
            some text
        
        Attributes
        ----------
        Blank : list
            some_text.
        """
        