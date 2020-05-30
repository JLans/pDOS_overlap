# -*- coding: utf-8 -*-
"""
Created on Tue May 19 22:28:09 2020

@author: lansf
"""

from __future__ import absolute_import, division, print_function
import numpy as np
from vasp_dos import get_band_center
from ase.io import read

def get_adsorbate_indices(GAS_CONTCAR, ADSORBATE_CONTCAR):
    gas_symbols = list(set(read(GAS_CONTCAR).get_chemical_symbols()))
    adsorbate_symbols = read(ADSORBATE_CONTCAR).get_chemical_symbols()
    adsorbate_indices = np.arange(len(adsorbate_symbols))[np.isin(adsorbate_symbols,gas_symbols)]
    return adsorbate_indices

class PDOS_OVERLAP:
    """Class for calculating adsorbate-surface relative orbital overlap"""
    def __init__(self, GAS_PDOS, REFERENCE_PDOS, adsorbate_indices=[]):
        """ 
        Parameters
        ----------
        GAS_PDOS : vasp_dos.VASP_DOS
            VASP_DOS object of gas phase calculation of adsorbate
            
        REFERENCE_PDOS : vasp_dos.VASP_DOS
            VASP_DOS object of reference adsorbate+surface system
            
        adsorbate_indices : list[int]
            index (indices) of adsorbate atom(s) in REFERENCE_PDOS
            
        site_indices : list[int]
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
            = self._get_mol_orbitals(REFERENCE_PDOS, adsorbate_indices)
            
        gas_orbital_features\
            = self._get_orbital_features(GAS_PDOS, gas_orbital_indices, gas_index)
        
        adsorbate_orbital_features\
            = self._get_orbital_features(REFERENCE_PDOS, adsorbate_orbital_indices\
                              , adsorbate_indices)
        
        orbital_scores = self._get_orbital_scores(gas_orbital_features\
                                                 , adsorbate_orbital_features)
        
        self.gas_orbitals = gas_orbitals
        self.gas_peak_energies = gas_peak_energies
        self.gas_peak_densities = gas_peak_densities
        self.adsorbate_orbitals = adsorbate_orbitals
        self.adsorbate_peak_energies = adsorbate_peak_energies
        self.adsorbate_peak_densities = adsorbate_peak_densities
        self.orbital_scores = orbital_scores
    
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
        #obtain all local maxima (including endpoints) - True or False
        peaks = np.r_[True, TOTAL_DOS[1:] > TOTAL_DOS[:-1]]\
                & np.r_[TOTAL_DOS[:-1] > TOTAL_DOS[1:], True]
                
        #endpoints should be zero density
        all_indices = np.arange(PDOS.ndos)
        peak_indices = all_indices[peaks]
        num_peaks = peak_indices.size
        peak_energies = PDOS.get_energies()[peak_indices]
        peak_densities = TOTAL_DOS[peak_indices]
        if num_peaks == 1:
            mol_orbitals = TOTAL_DOS
            orbital_indices = all_indices
        else:
            mol_orbitals = np.zeros((num_peaks,TOTAL_DOS.size))
            orbital_indices = []
            midpoints = (0.5 * peak_indices[0:-1] + 0.5 * peak_indices[1:]).astype(int)
            mol_orbitals[0][0:midpoints[0]] = TOTAL_DOS[0:midpoints[0]]
            orbital_indices.append(all_indices[0:midpoints[0]])
            for i in range(midpoints.size-1):
                mol_orbitals[i+1][midpoints[i]:midpoints[i+1]] = TOTAL_DOS[midpoints[i]:midpoints[i+1]]
                orbital_indices.append(all_indices[midpoints[i]:midpoints[i+1]])
            mol_orbitals[-1][midpoints[-1]:] = TOTAL_DOS[midpoints[-1]:]
            orbital_indices.append(all_indices[midpoints[-1]:])
        
        orbital_occupations = np.trapz(mol_orbitals, PDOS.get_energies(), axis=1)
        new_mol_orbitals = []
        new_orbital_indices = []
        new_occupations = []
        orbital_num = 0
        #theoreticly new_occupations < 1 but account for some loss of electrons
        while orbital_num < num_peaks:
            new_mol_orbitals.append(mol_orbitals[orbital_num])
            new_orbital_indices.append(orbital_indices[orbital_num])
            new_occupations.append(orbital_occupations[orbital_num])
            orbital_num += 1
            while orbital_num < num_peaks\
            and orbital_occupations[orbital_num - 1] < 0.9\
            and orbital_occupations[orbital_num] < 0.9:
                new_mol_orbitals[-1] += mol_orbitals[orbital_num]
                new_orbital_indices[-1]\
                                    = np.concatenate( ( new_orbital_indices[-1]\
                                    , orbital_indices[orbital_num] ) )
                new_occupations[-1] += orbital_occupations[orbital_num]
                orbital_num += 1
                
        for count, i in enumerate(new_mol_orbitals):
            plt.figure(count)
            plt.plot(PDOS.get_energies(), i)
            plt.show()
               
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
        orbital_list = list(set([key[0] for key in PDOS.orbital_dictionary.keys()]))
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
            for count2, feature_level in enumerate(TOTAL_PDOS):
                orbital_features[count][count2+1] = np.trapz(feature_level[index_values], energies[index_values])
                       
        return orbital_features
    
    def _get_orbital_scores(self, gas_orbital_features, adsorbate_orbital_features):
        """ orbital scores which represent probabilities of the gas molecular
            orbital matching each adsorbate molecular orbital
        
        Parameters
        ----------
        gas_orbital_features : numpy.ndarray
            gas orbital features
            
        adsorbate_orbital_features : numpy.ndarray
            adsorbate orbital features
                   
        Returns
        -------
        orbital_scores : numpy.ndarray
            M x N 2D array where M is the number of gas molecular orbitals
            and N is the number of adsorbate molecular orbitals

        """
               
        num_gas_orbitals = gas_orbital_features.shape[0]
        num_adsorbate_orbitals = adsorbate_orbital_features.shape[0]
        orbital_scores_gas = np.zeros((num_gas_orbitals,num_adsorbate_orbitals))
        orbital_scores_adsorbate = np.zeros((num_adsorbate_orbitals,num_gas_orbitals))
                          
        for count in range(num_gas_orbitals):
            diff_squared = (adsorbate_orbital_features - gas_orbital_features[count])**2
            var_feature = np.sum(diff_squared,axis=0)/num_adsorbate_orbitals
            var_feature[var_feature[...] == 0] = 10**-10 #prevent divide by zero error
            orbital_scores_gas[count] = np.prod( np.exp(-1 * diff_squared / var_feature), axis=1)
            orbital_scores_gas[count] /= orbital_scores_gas.sum()
            
        for count in range(num_adsorbate_orbitals):
            diff_squared = (gas_orbital_features - adsorbate_orbital_features[count])**2
            var_feature = np.sum(diff_squared,axis=0)/num_gas_orbitals
            var_feature[var_feature[...] == 0] = 10**-10 #prevent divide by zero error
            orbital_scores_adsorbate[count] = np.prod( np.exp(-1 * diff_squared / var_feature), axis=1)
            orbital_scores_adsorbate[count] /= orbital_scores_adsorbate.sum()
        
        orbital_scores = orbital_scores_gas * orbital_scores_adsorbate.T
        
        orbital_scores /= orbital_scores.sum(axis=1, keepdims=True)
        
        return orbital_scores
    