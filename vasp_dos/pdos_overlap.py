# -*- coding: utf-8 -*-
"""
Created on Tue May 19 22:28:09 2020

@author: lansf
"""

from __future__ import absolute_import, division, print_function
import numpy as np
from vasp_dos import get_band_center



class PDOS_OVERLAP:
    """Class for calculating adsorbate-surface relative orbital overlap"""
    def __init__(self, GAS_PDOS, REFERENCE_PDOS, adsorbate_index=[], site_index=[]):
        """ 
        Parameters
        ----------
        GAS_PDOS : vasp_dos.VASP_DOS
            VASP_DOS object of gas phase calculation of adsorbate
            
        REFERENCE_PDOS : vasp_dos.VASP_DOS
            VASP_DOS object of reference adsorbate+surface system
            
        adsorbate_index : list[int]
            index (indices) of adsorbate atom(s) in REFERENCE_PDOS
            
        site_index : list[int]
            index (indices) of adsorbate atom(s) in REFERENCE_PDOS
        
        Attributes
        ----------
        emax : float
            maximum energy level
            
        Notes
        -----
        GAS_PDOS is used only to determine the number of orbitals that can
        interact with the surface and to calculate relative orbital overlap with 
        projected density of adsorption sites without adsorbates.
        """
        gas_index = list(range(GAS_PDOS.natoms))
        
        
        gas_peak_energies, gas_peak_densities, gas_orbitals\
            , gas_orbital_indices = self._get_mol_orbitals(GAS_PDOS, gas_index)
        
        adsorbate_peak_energies, adsorbate_peak_densities, adsorbate_orbitals \
            , adsorbate_orbital_indices \
            = self._get_mol_orbitals(REFERENCE_PDOS, adsorbate_index)
            
        gas_orbital_features\
            = self._get_orbital_features(GAS_PDOS, gas_orbital_indices, gas_index)
        
        adsorbate_orbital_features\
            = self._get_orbital_features(REFERENCE_PDOS, adsorbate_orbital_indices\
                              , adsorbate_index)
        
        orbital_errors = self._get_orbital_scores(gas_orbital_features\
                                                 , adsorbate_orbital_features)
            
        orbital_scores = self._get_orbital_scores(orbital_errors)
        
        self.gas_orbitals = gas_orbitals
        self.gas_peak_energies = gas_peak_energies
        self.gas_peak_densities = gas_peak_densities
        self.orbital_scores = orbital_scores
        
        
        #Use distance between similarity to get probability
        #take logit of probability https://en.wikipedia.org/wiki/Logit
        #take softmax of logit functions https://en.wikipedia.org/wiki/Softmax_function
        #find which adsorbate orbital matches a give gas orbital (give adsorbate orbital can match more than one gas orbital)
    
    def _get_mol_orbitals(self, PDOS, atom_indices):
        """ molecular orbitals as an M x ndos array
        
        Parameters
        ----------
        PDOS : vasp_dos.VASP_DOS
            VASP_DOS object of gas or surface phase calculation of adsorbate
            
        atom_indices : list
            indices of atoms to include in the TOTAL DOS. If empty, then all
            atoms are included (as for the gas)
        
        Returns
        -------
        peak_energies : float or numpy.ndarray
            peak energy(s) of molecular orbital(s)
            
        peak_densities : float or numpy.ndarray
            peak density(ies) of molecular orbital(s)
            
        mol_orbitals : numpy.ndarray
            M x ndos 2D array where M is the number of molecular orbitals and
            ndos is the number of gridpoints for the density of states
        
        orbital_indices : list[numpy.ndarray]
            Length M list of non-zero orbital indices

        """
        
        TOTAL_DOS = np.zeros(PDOS.ndos)
        #sum over projected density of states for each atom
        orbital_list = [key for key in PDOS.orbital_dictionary.keys()]
        for atom in atom_indices:
            orbital_list, projected_dos = PDOS.get_site_dos(atom, orbital_list)
            TOTAL_DOS += projected_dos.sum(axis=0)
        #obtain all local maxima (including endpoints)
        peaks = np.r_[True, TOTAL_DOS[1:] > TOTAL_DOS[:-1]]\
                & np.r_[TOTAL_DOS[:-1] > TOTAL_DOS[1:], True]
        #endpoints should be zero density
        all_indices = np.arange(PDOS.ndos)
        peak_indices = all_indices[peaks]
        peak_energies = PDOS.get_energies()[peak_indices]
        peak_densities = TOTAL_DOS[peak_indices]
        if peak_indices.size == 1:
            mol_orbitals = TOTAL_DOS
            orbital_indices = all_indices
        else:
            mol_orbitals = np.zeros((peak_indices.size,TOTAL_DOS.size))
            orbital_indices = []
            midpoints = (0.5 * peak_indices[0:-1] + 0.5 * peak_indices[1:]).astype(int)
            mol_orbitals[0][0:midpoints[0]] = TOTAL_DOS[0:midpoints[0]]
            orbital_indices.append(all_indices[0:midpoints[0]])
            for i in range(midpoints.size-1):
                mol_orbitals[i+1][midpoints[i]:midpoints[i+1]] = TOTAL_DOS[midpoints[i]:midpoints[i+1]]
                orbital_indices.append(all_indices[midpoints[i]:midpoints[i+1]])
            mol_orbitals[-1][midpoints[-1]:] = TOTAL_DOS[midpoints[-1]:]
            orbital_indices.append(all_indices[midpoints[-1]:])
               
        return peak_energies, peak_densities, mol_orbitals, orbital_indices
    
    def _get_orbital_features(self, PDOS, orbital_indices, atom_indices=[]):
        """ molecular orbitals as an M x ndos array
        
        Parameters
        ----------
        PDOS : vasp_dos.VASP_DOS
            VASP_DOS object of gas or surface phase calculation of adsorbate
            
        atom_indices : list
            indices of atoms to include in the TOTAL DOS. If empty, then all
            atoms are included (as for the gas)
            
        orbital_indices : list[numpy.ndarray]
            Length M list of non-zero orbital indices
        
        Returns
        -------
        orbital_features : numpy.ndarray
            M x n_features 2D array where M is the number of molecular orbitals
            and n_features is the number of orbital features for calculating
            orbital similarity. Includes molecular orbital energy center and
            integrated s, p, and d molecular orbital density of states

        """
        
        num_orbitals = len(orbital_indices)
        energies = PDOS.get_energies()
        
        orbital_list = [key for key in ['s', 'p', 'd', 'f',] if key in PDOS.orbital_dictionary.keys()]
        feature_levels = ['s','p','d']
        TOTAL_DOS = np.zeros(PDOS.ndos)
        TOTAL_PDOS = np.zeros((3,PDOS.ndos))
        #sum over projected density of states for each atom
        for atom in atom_indices:
            orbital_list, projected_dos = PDOS.get_site_dos(atom, orbital_list, sum_density=True)
            TOTAL_DOS += projected_dos.sum(axis=0)
            for count, value in enumerate(feature_levels):
                if value in orbital_list:
                    TOTAL_PDOS[count] += projected_dos[orbital_list.index(value)]
        
        orbital_features = np.zeros((num_orbitals,4))
        for count, index_values in enumerate(orbital_indices):
            orbital_features[count][0] = get_band_center(energies[index_values],TOTAL_DOS[index_values])
            for count2, feature_level in TOTAL_PDOS:
                orbital_features[count][count2+1] = np.trapz(feature_level[index_values], energies[index_values])
                       
        return orbital_features
    
    def _get_orbital_errors(self, gas_orbital_features, adsorbate_orbital_features):
        """ Relative molecular orbital feature errors as an M x N array
        
        Parameters
        ----------
        gas_orbital_features : numpy.ndarray
            gas orbital features
            
        adsorbate_orbital_features : numpy.ndarray
            adsorbate orbital features
                   
        Returns
        -------
        orbital_errors : numpy.ndarray
            M x N 2D array where M is the number of gas molecular orbitals
            and N is the number of adsorbate molecular orbitals

        """
               
        num_gas_orbitals = gas_orbital_features.shape[0]
        num_ad_orbitals = adsorbate_orbital_features.shape[0]
        num_features = num_gas_orbitals.shape[1]
        orbital_errors = np.zeros((num_gas_orbitals,num_ad_orbitals))
        
        min_features = np.min(np.concatenate((gas_orbital_features\
                                              , adsorbate_orbital_features)\
                                             , axis=0), axis=0)
            
        max_features = np.max(np.concatenate((gas_orbital_features\
                                              , adsorbate_orbital_features)\
                                             , axis=0), axis=0)
            
        feature_range = max_features - min_features
        feature_range[feature_range[...] == 0] =1
        
        for count in num_gas_orbitals:
            diff = adsorbate_orbital_features - gas_orbital_features[count]
            norm_diff = diff / feature_range
            RMSE = np.sum(norm_diff**2/num_features, axis=1)**0.5
            orbital_errors[count] = RMSE
            
        return orbital_errors
            
    def _get_orbital_scores(self, orbital_errors):
        """ orbital scores which represent probabilities of the gas molecular
        orbital matching each adsorbate molecular orbital
        
        Parameters
        ----------
        orbital_errors : numpy.ndarray
            M x N 2D array where M is the number of gas molecular orbitals
            and N is the number of adsorbate molecular orbitals
            
        Returns
        -------
        orbital_scores : numpy.ndarray
            M x N 2D array where M is the number of gas molecular orbitals
            and N is the number of adsorbate molecular orbitals
            
            
        Notes
        
        p = 1 - orbital_errors
        
        To get probability that an adsorbate molecular orbital matches a gas
        molecular orbital:
        
        take logit of probability https://en.wikipedia.org/wiki/Logit
        take softmax of logit functions https://en.wikipedia.org/wiki/Softmax_function
        logit(p) = ln(p/(1-p)) --> maps p to -infinity to +infinity
        softmax = exp(zi)/sum(exp(zj))
        softmax(logits) = p_i/(1-p_i) / sum(p_j/(1-p_j))

        """
           
        pairwise_prob = 1 - orbital_errors
        orbital_errors[orbital_errors[...] == 0] = 10**-8 # some small number
        odds = pairwise_prob/orbital_errors
        orbital_scores = odds/odds.sum(axis=1)

        return orbital_scores
    