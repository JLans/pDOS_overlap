# -*- coding: utf-8 -*-
"""
Created on Tue May 19 22:28:09 2020

@author: lansf
"""

from __future__ import absolute_import, division, print_function
import os
import pkg_resources
import numpy as np

def get_example_data():
    """ Get default paths to experimental data.
    
    Returns
    -------
    data_path : str
        path to example lobster data
    """
    data_path = pkg_resources.resource_filename(__name__, 'data/lobster')
    return data_path

def get_all_lobster_files(directory, file_type='COOPCAR.lobster'):
    """ Get all DOSCAR and CONTCAR file paths
    
    Parameters
    ----------
    directory : str
        root directory to look for DOSCAR files
        
    file_type : string
        file type for which to search and return file path
    
    Returns
    -------
    lobster_files : list[str]
        list of paths to lobster files of type file_type
    """
    lobster_directories = [os.path.join(r,subdirectory) for r,d,f in os.walk(directory) \
              for subdirectory in d \
              if file_type in os.listdir(os.path.join(r,subdirectory))]
    lobster_files = [os.path.join(d,file_type) for d in lobster_directories]
    return lobster_files

class OVERLAP_POPULATION:
    """Class for analyzing overlap population bonding data"""
    def __init__(self, file_name="COOPCAR.lobster"):
        """ 
        Parameters
        ----------
        file_name : str
            full lobster file location
        
        Attributes
        ----------
        emax : float
            maximum energy level
            
        emin : float
            minimum energy level
            
        ndos : int
            number of descritized energy levels
            
        e_fermi : float
            highest occupied energy level
            
        is_spin : bool
            indicates if projected density is spin resolved
            
        num_interactions : int
            number of interactions in COOP
        
        interactions : list[str]
            list of the interactions included in the file
            
        """
        
        #conditional read statements can be added
        num_interactions, interactions, is_spin, ndos, e_fermi, e_min, e_max\
            = self._read_coopcar(file_name=file_name)
        self.num_interactions = num_interactions
        self.interactions = interactions
        self.is_spin = is_spin
        self.ndos = ndos
        self.e_fermi = e_fermi
        self.e_min = e_min
        self.e_max = e_max
        
    def get_energies(self):
        """ method for obtaining energy levels
        
        Returns
        -------
        energies : numpy.ndarray
            1-D array of energies    
        """
        energies = self._pcoop[0,:].copy()
        return energies
    
    def get_average_pcoop(self, sum_spin=True):
        """ method for obtaining total density of states
        
        Returns
        -------
        average_pcoop : numpy.ndarray
            1-D or 2-D array of average pcoop for all interactions
        """
        if self.is_spin == True:
            if sum_spin == False:
                average_pcoop = self._pcoop[[1, 2 * self.num_interactions + 3], :]
            else:
                average_pcoop = self._pcoop[1, :] + self._pcoop[2 * self.num_interactions + 3, :]
        else:
            average_pcoop = self._pcoop[1, :]
        return average_pcoop

    def get_average_int_pcoop(self, sum_spin=True):
        """ method for obtaining total integrated density of states
        
        Returns
        -------
        average_int_pcoop : numpy.ndarray
            1-D or 2-D array of integrated pcoop for all interactions
        """
        if self.is_spin == True:
            if sum_spin == False:
                average_int_pcoop = self._pcoop[[2, 2 * self.num_interactions + 4], :]
            else:
                average_int_pcoop = self._pcoop[2, :] + self._pcoop[2 * self.num_interactions + 4, :]
        else:
            average_int_pcoop = self._pcoop[2, :]
        return average_int_pcoop
    
    def get_integrated_pcoop(self, interactions=[], sum_pcoop=False, sum_spin=True):
        if len(interactions) == 0:
            interactions = list(range(self.num_interactions))
        if self.is_spin == True:
            spin_up = self._pcoop[4:2 * self.num_interactions + 3:2, :][interactions]
            spin_down = self._pcoop[3 * self.num_interactions + 4::2, :][interactions]
            if sum_spin == True:
                integrated_pcoop = spin_up + spin_down
            else:
                integrated_pcoop = np.array([spin_up, spin_down])
        else:
            integrated_pcoop = self._pcoop[4::2, :][interactions]
        if sum_pcoop == True or len(interactions) == 1:
            axis = len(integrated_pcoop.shape) - 2
            integrated_pcoop = integrated_pcoop.sum(axis=axis)
        return integrated_pcoop
    
    def get_pcoop(self, interactions=[], sum_pcoop=False, sum_spin=True):
        if len(interactions) == 0:
            interactions = list(range(self.num_interactions))
        if self.is_spin == True:
            spin_up = self._pcoop[3:2 * self.num_interactions + 3:2, :][interactions]
            spin_down = self._pcoop[2 * self.num_interactions + 5::2, :][interactions]
            if sum_spin == True:
                pcoop = spin_up + spin_down
            else:
                pcoop = np.array([spin_up, spin_down])
        else:
            pcoop = self._pcoop[3::2, :][interactions]
        if sum_pcoop == True or len(interactions) == 1:
            axis = len(pcoop.shape) - 2
            pcoop = pcoop.sum(axis=axis)
        return pcoop
        
    def _read_coopcar(self, file_name="COOPCAR.lobster"):
        """Read lobster COOPCAR and extract projected overlap
        
        Parameters
        ----------
        file_name : str
            file location of the COOPCAR.lobster file file
            
        Attributes
        ----------
        _pcoop : numpy.ndarray
            numpy array that contains the energy of levels and the projected
            crystal orbital overlap population densities
            
        Returns
        -------
        emax : float
            maximum energy level
            
        emin : float
            minimum energy level
            
        ndos : int
            number of descritized energy levels
            
        e_fermi : float
            highest occupied energy level
            
        is_spin : bool
            indicates if projected density is spin resolved
            
        num_interactions : int
            number of interactions in COOP
        
        interactions : list[str]
            list of the interactions included in the file
        """
        f = open(file_name)
        f.readline() #skip the first line
        descriptive_line = f.readline().split()
        num_interactions = int(descriptive_line[0]) - 1
        is_spin = int(descriptive_line[1])
        if is_spin == 2:
            is_spin = True
        else:
            is_spin = False
        ndos = int(descriptive_line[2])
        e_fermi = float(descriptive_line[5])
        e_min = float(descriptive_line[3])
        e_max = float(descriptive_line[4])
        f.readline() #skip the line saying average
        interactions = []
        for i in range(num_interactions):
            interactions += f.readline().split()
        line = f.readline().split()
        pcoop = np.zeros((ndos,len(line)))
        pcoop[0] = np.array(line)
        for nd in range(1,ndos):
            line = f.readline().split()
            pcoop[nd] = np.array(line)
        pcoop = pcoop.T
        pcoop[0] += e_fermi
        self._pcoop = pcoop
        return num_interactions, interactions, is_spin, ndos, e_fermi, e_min, e_max
        
        