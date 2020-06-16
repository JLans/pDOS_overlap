# -*- coding: utf-8 -*-
"""
Created on Tue May 19 22:28:09 2020

@author: lansf
"""

from __future__ import absolute_import, division, print_function
import numpy as np

class LOBSTER:
    """Class for analyzing LOBSTER bonding data"""
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
        
        num_interactions, interactions, is_spin, ndos, e_fermi, e_min, e_max\
            = self._read_coopcar(file_name=file_name)
        self.num_interactions = num_interactions
        self.interactions = interactions
        self.is_spin 
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
            spin_up = self._pcoop[3:2 * self.num_interactions + 3:2, :][interactions, :]
            spin_down = self._pcoop[3 * self._num_interactions + 3::2, :][interactions, :]
            if sum_spin == True:
                integrated_pcoop = spin_up + spin_down
            else:
                integrated_pcoop = np.array([spin_up, spin_down])
        else:
            integrated_pcoop = self._pcoop[3::2, :][interactions, :]
        if sum_pcoop == True:
            integrated_pcoop = integrated_pcoop.sum(axis=-1)
        return integrated_pcoop
    
    def get_pcoop(self, interactions=[], sum_pcoop=False, sum_spin=True):
        if len(interactions) == 0:
            interactions = list(range(self.num_interactions))
        if self.is_spin == True:
            spin_up = self._pcoop[3:2 * self.num_interactions + 3:2, :][interactions, :]
            spin_down = self._pcoop[3 * self._num_interactions + 3::2, :][interactions, :]
            if sum_spin == True:
                pcoop = spin_up + spin_down
            else:
                pcoop = np.array([spin_up, spin_down])
        else:
            pcoop = self._pcoop[3::2, :][interactions, :]
        if sum_pcoop == True:
            pcoop = pcoop.sum(axis=-1)
        return pcoop
        