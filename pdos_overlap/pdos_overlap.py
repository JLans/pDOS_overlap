# -*- coding: utf-8 -*-
"""
Created on Tue May 19 22:28:09 2020

@author: lansf
"""

from __future__ import absolute_import, division, print_function
import numpy as np
import matplotlib.pyplot as plt
from .vasp_dos import get_band_center
from ase.io import read
from .coordination import Coordination
from .plotting_tools import set_figure_settings

def get_adsorbate_indices(GAS_CONTCAR, ADSORBATE_CONTCAR):
    gas_symbols = list(set(read(GAS_CONTCAR).get_chemical_symbols()))
    adsorbate_symbols = read(ADSORBATE_CONTCAR).get_chemical_symbols()
    adsorbate_indices = np.arange(len(adsorbate_symbols))[np.isin(adsorbate_symbols,gas_symbols)]
    CN = Coordination(read(ADSORBATE_CONTCAR),exclude=adsorbate_indices,cutoff=1.25)
    bonded = CN.get_bonded()
    site_indices = []
    for i in adsorbate_indices:
        site_indices += bonded[i]
    return adsorbate_indices, site_indices

class PDOS_OVERLAP:
    """Class for calculating adsorbate-surface relative orbital overlap"""
    def __init__(self, GAS_PDOS, REFERENCE_PDOS, adsorbate_indices, site_indices\
                 , min_occupation=0.9, sum_density=False, sum_spin=True):
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
        gas_indices = list(range(GAS_PDOS.natoms))
        #GCN = [CN.get_gcn(CN.bonded[i]) for i in Catoms]
        gas_peak_energies, gas_peak_densities, gas_orbitals\
            , gas_orbital_indices, gas_occupations\
        = self._get_mol_orbitals(GAS_PDOS, gas_indices, min_occupation=min_occupation)
        
        adsorbate_peak_energies, adsorbate_peak_densities, adsorbate_orbitals \
            , adsorbate_orbital_indices, adsorbate_occupations \
            = self._get_mol_orbitals(REFERENCE_PDOS, adsorbate_indices, min_occupation=min_occupation)
            
        overlap_orbitals, energy_overlap \
        = self._calculate_overlap(adsorbate_orbitals, REFERENCE_PDOS, site_indices
                                 , sum_density=sum_density, sum_spin=sum_spin) 
        
        if gas_orbitals.shape[0] == adsorbate_orbitals.shape[0]:
            feature_type = 'states'
        else:
            feature_type = 'state fractions'
        
        gas_features\
            = self._get_orbital_features(GAS_PDOS, gas_orbital_indices\
                                , gas_indices, feature_type, sum_spin=sum_spin)
        
        adsorbate_features\
            = self._get_orbital_features(REFERENCE_PDOS, adsorbate_orbital_indices\
                              , adsorbate_indices, feature_type, sum_spin=sum_spin)
        
        orbital_scores_gas, orbital_scores_adsorbate, orbital_scores_averaged\
        = self._get_orbital_scores(gas_features, adsorbate_features)
        
        if gas_orbitals.shape[0] == adsorbate_orbitals.shape[0]:
            orbital_scores = self.get_pair_scores(orbital_scores_averaged)
        elif gas_orbitals.shape[0] > adsorbate_orbitals.shape[0]:
            orbital_scores = orbital_scores_gas
        else:
            orbital_scores = orbital_scores_adsorbate
        
        self.gas_band_centers = get_band_center(GAS_PDOS.get_energies(), gas_orbitals)
        self.adsorbate_band_centers = get_band_center(\
                                    REFERENCE_PDOS.get_energies(), adsorbate_orbitals)
        self.gas_2_adsorbate = self._assign_orbitals(orbital_scores\
                                    , self.gas_band_centers, self.adsorbate_band_centers)
        
        self.gas_features = gas_features
        self.adsorbate_features = adsorbate_features
        self.gas_orbitals = gas_orbitals
        self.gas_occupations = gas_occupations
        self.adsorbate_orbitals = adsorbate_orbitals
        self.adsorbate_occupations = adsorbate_occupations
        self.orbital_scores = orbital_scores
        self.overlap_orbitals = overlap_orbitals
        self.energy_overlap = energy_overlap
        self.GAS_PDOS = GAS_PDOS
        self.REFERENCE_PDOS = REFERENCE_PDOS
        self.gas_indices = gas_indices
        self.adsorbate_indices = adsorbate_indices
        self.site_indices = site_indices
        self.min_occupation = min_occupation
        self.sum_density = sum_density
        self.sum_spin = sum_spin
    
    @staticmethod
    def _get_mol_orbitals(PDOS, atom_indices, min_occupation=0.9):
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
        orbital_list, projected_dos = PDOS.get_site_dos(atom_indices, orbital_list)
        TOTAL_DOS = projected_dos.sum(axis=0)
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
            and orbital_occupations[orbital_num - 1] < min_occupation\
            and orbital_occupations[orbital_num] < min_occupation:
                new_mol_orbitals[-1] += mol_orbitals[orbital_num]
                new_orbital_indices[-1]\
                                    = np.concatenate( ( new_orbital_indices[-1]\
                                    , orbital_indices[orbital_num] ) )
                new_occupations[-1] += orbital_occupations[orbital_num]
                orbital_num += 1
        num_new_orbitals = len(new_mol_orbitals)
        band_centers = get_band_center(PDOS.get_energies(),np.array(new_mol_orbitals))
        
        num_orbitals_final = len([value for value in new_occupations if value > min_occupation])
        mol_orbitals_final = np.zeros((num_orbitals_final,TOTAL_DOS.size))
        orbital_indices_final = [np.array([]).astype(int) for _ in range(num_orbitals_final)]
        occupations_final = np.zeros(num_orbitals_final)
        
        j = 0
        for i in range(num_new_orbitals):
            if new_occupations[i] > min_occupation or i == 0:
                mol_orbitals_final[j] += new_mol_orbitals[i]
                orbital_indices_final[j] = np.concatenate(
                    (orbital_indices_final[j],new_orbital_indices[i]) )
                occupations_final[j] += new_occupations[i]
                #only update if a complete orbital is found    
                if new_occupations[i] > min_occupation:
                    j += 1
            elif i == num_new_orbitals - 1:
                mol_orbitals_final[-1] += new_mol_orbitals[i]
                orbital_indices_final[-1] = np.concatenate(
                    (orbital_indices_final[-1],new_orbital_indices[i]) )
                occupations_final[-1] += new_occupations[i]
            elif abs(band_centers[i] - band_centers[i-1])\
                <= abs(band_centers[i] - band_centers[i+1]):
                mol_orbitals_final[j-1] += new_mol_orbitals[i]
                orbital_indices_final[j-1] = np.concatenate(
                    (orbital_indices_final[j-1],new_orbital_indices[i]) )
                occupations_final[j-1] += new_occupations[i]
            elif abs(band_centers[i] - band_centers[i-1])\
                > abs(band_centers[i] - band_centers[i+1]):
                mol_orbitals_final[j] += new_mol_orbitals[i]
                orbital_indices_final[j] = np.concatenate(
                    (orbital_indices_final[j],new_orbital_indices[i]) )
                occupations_final[j] += new_occupations[i]
                
        molecular_orbitals = mol_orbitals_final
        orbital_indices = orbital_indices_final
        occupations = occupations_final
               
        return peak_energies, peak_densities, molecular_orbitals, orbital_indices, occupations
    
    
    @staticmethod
    def _get_orbital_features(PDOS, orbital_indices, atom_indices=[]\
                              , feature_type='states', sum_spin=True):
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
        energies = PDOS.get_energies() - PDOS.e_fermi
        orbital_list = list(PDOS.orbital_dictionary.keys())
        #feature_list = ['s','py', 'pz', 'px', 'dxy', 'dyz', 'dz2', 'dxz', 'dx2-y2']
        feature_list = ['s','py', 'pz', 'px']
        if sum_spin == False:
            orbital_list=['s+', 's-', 'py+', 'py-', 'pz+', 'pz-', 'px+', 'px-']
        num_features = len(feature_list) + 1
        TOTAL_DOS = np.zeros(PDOS.ndos)
        TOTAL_PDOS = np.zeros((num_features - 1,PDOS.ndos))
        #sum over projected density of states for each atom
        orbital_list, projected_dos = PDOS.get_site_dos(atom_indices\
                    , orbital_list, sum_density=False, sum_spin=sum_spin)
        TOTAL_DOS = projected_dos.sum(axis=0)
        for count, value in enumerate(feature_list):
            if value in orbital_list:
                TOTAL_PDOS[count] += projected_dos[orbital_list.index(value)]
                
        orbital_features = np.zeros((num_orbitals,num_features))
        for count, index_values in enumerate(orbital_indices):
            orbital_features[count][0] = get_band_center(energies[index_values],TOTAL_DOS[index_values])
            for count2, feature_level in enumerate(TOTAL_PDOS):
                if feature_type == 'states':
                    orbital_features[count][count2+1]\
                    = np.trapz(feature_level[index_values], energies[index_values])
                    
                elif feature_type == 'state fractions':
                    orbital_features[count][count2+1]\
                    = np.trapz(feature_level[index_values], energies[index_values])\
                      / np.trapz(TOTAL_DOS[index_values], energies[index_values])
        return orbital_features
    
    @staticmethod
    def _get_orbital_scores(gas_orbital_features, adsorbate_orbital_features):
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

        Notes
        -----
        Matrix scaling perfomed using RAS method described in 
        DOI: 10.1109/FOCS.2017.87
        """
        def get_scores(reference_features, comparison_features):
            num_reference_orbitals = reference_features.shape[0]
            num_comparison_orbitals = comparison_features.shape[0]
            orbital_scores = np.zeros((num_reference_orbitals, num_comparison_orbitals))
            for count in range(num_reference_orbitals):
                diff_squared = (comparison_features - reference_features[count])**2
                var_feature = np.sum(diff_squared,axis=0) / num_comparison_orbitals
                var_feature[var_feature[...] <= 10**-6] = 10 #prevent divide by zero error
                #gaussian based likelihood
                diff_squared[:,0] *=3
                univariate = np.exp(-1 * diff_squared / var_feature)
                orbital_scores[count] = np.prod(univariate, axis=1)
            return orbital_scores
        
        orbital_scores_gas = get_scores(gas_orbital_features, adsorbate_orbital_features)
        orbital_scores_adsorbate = get_scores(adsorbate_orbital_features, gas_orbital_features)
        
        #normalize
        orbital_scores_gas /= orbital_scores_gas.sum(axis=1, keepdims=True)
        orbital_scores_adsorbate /= orbital_scores_adsorbate.sum(axis=1, keepdims=True)
        orbital_scores_averaged = (orbital_scores_gas * orbital_scores_adsorbate.T)**0.5
        orbital_scores_averaged /= orbital_scores_averaged.sum(axis=1, keepdims=True)
        
        return orbital_scores_gas, orbital_scores_adsorbate, orbital_scores_averaged
    
    @staticmethod
    def get_pair_scores(orbital_scores, max_iterations = 500, verbose=False):
        error = 1
        iteration = 0
        orbital_scores = orbital_scores.copy()
        while error > 10**-7 and iteration < max_iterations:
            orbital_scores = (orbital_scores * orbital_scores.T)**0.5
            orbital_scores /= orbital_scores.sum(axis=0, keepdims=True)
            orbital_scores /= orbital_scores.sum(axis=1, keepdims=True)
            error = (abs(1 - orbital_scores.sum(axis=0))).max()
            iteration += 1
        if verbose == True:
            print('The max error is ' + str(error))
            print('Number of iterations is ' + str(iteration))
        return orbital_scores
    
    def set_orbital_scores(self,orbital_scores):
        self.orbital_scores = orbital_scores
        
    @staticmethod
    def _calculate_overlap(molecular_orbitals, PDOS, site_indices\
                           , sum_density=False, sum_spin=True):            
        
        num_molecular_orbitals = molecular_orbitals.shape[0]
        atomic_orbitals = list(PDOS.orbital_dictionary.keys())
        orbital_list, TOTAL_PDOS = PDOS.get_site_dos(site_indices, atomic_orbitals\
                                        , sum_density=sum_density, sum_spin=sum_spin)
        #calculate sumj(Mij*projected_dosj) for all i  
        overlap_orbitals = orbital_list
        energy_overlap = np.zeros((num_molecular_orbitals, TOTAL_PDOS.shape[0]))
        for i in range(num_molecular_orbitals):
            energy_overlap[i] = np.trapz(molecular_orbitals[i]*TOTAL_PDOS,PDOS.get_energies(),axis=1)
        return overlap_orbitals, energy_overlap
    
    def calculate_orbital_interaction(self,PDOS, atoms, gas_orbital_index, atomic_orbitals\
                                      ,sum_density=False, sum_spin=True):
        GAS_PDOS = self.GAS_PDOS
        gas_orbital = self.gas_orbitals[gas_orbital_index]
        orbitals, projected_density = PDOS.get_site_dos(atoms\
                                    , atomic_orbitals, sum_density=sum_density\
                                    , sum_spin=sum_spin)
        if PDOS.ndos != GAS_PDOS.ndos:
            projected_density = np.zeros((len(orbitals),GAS_PDOS.ndos))
            old_x = PDOS.get_energies()
            new_x = np.linspace(old_x.min(), old_x.max(), GAS_PDOS.ndos)
            for i in range(len(orbitals)):
                _, temp_density = PDOS.get_site_dos(list(range(PDOS.natoms))\
                                    , [atomic_orbitals[i]], sum_density=sum_density\
                                    , sum_spin=sum_spin)
                projected_density[i] = np.interp(new_x, old_x, temp_density)
        gas_band_center = get_band_center(GAS_PDOS.get_energies(), gas_orbital)
        
        orbital_interaction = np.trapz(projected_density\
                          / np.abs(gas_band_center - new_x)\
                              ,PDOS.get_energies(), axis=1)
        return orbital_interaction
    
    @staticmethod
    def _assign_orbitals(orbital_scores, gas_band_centers, adsorbate_band_centers):
        orbital_assignments = np.argmax(orbital_scores,axis=1).astype('int')
        indices = np.arange(orbital_assignments.size).astype('int')
        if gas_band_centers.shape >= adsorbate_band_centers.shape:
            gas_2_adsorbate = np.vstack((indices,orbital_assignments\
            , gas_band_centers, adsorbate_band_centers[orbital_assignments])).T
        else:
            gas_2_adsorbate = np.vstack((orbital_assignments,indices\
            , gas_band_centers[orbital_assignments], adsorbate_band_centers)).T
                
        return gas_2_adsorbate
                
    def plot_energy_overlap(self, indices = [...],settings='print'):
        if settings == 'presentation':
            set_figure_settings('presentation')
        else:
            set_figure_settings('paper')
        energy_overlap = self.energy_overlap[indices]
        if len(energy_overlap.shape) == 1:
            energy_overlap = [energy_overlap]
        for overlap in energy_overlap:
            plt.figure()
            xticks = np.arange(1,len(self.overlap_orbitals)+1)
            plt.bar(xticks,overlap)
            plt.xticks(xticks,labels=self.overlap_orbitals)
            plt.xlabel('Metal states projected onto atomic orbitals')
            plt.ylabel('Overlap with adsorbate molecular orbitals')
            if settings == 'print':
                plt.show()
            
    def plot_projected_density(self,sum_density=True, sum_spin=True, settings='print'):
        GAS_PDOS = self.GAS_PDOS
        gas_indices = self.gas_indices
        REFERENCE_PDOS = self.REFERENCE_PDOS
        adsorbate_indices = self.adsorbate_indices
        site_indices = self.site_indices
        if settings == 'presentation':
            set_figure_settings('presentation')
        else:
            set_figure_settings('paper')
        orbital_list=['s', 'p', 'd']
        colors = ['b-','g-','r-']
        zorder = [2,3,4]
        if sum_spin == False:
            orbital_list=['s+', 's-', 'p+', 'p-', 'd+', 'd-']
            colors = ['b-', 'b-', 'g-', 'g-', 'r-', 'r-']
            zorder = [2,3,4,5,6,7]
        
        #plotting function
        def plot_density(PDOS, indices, title):
            orbitals, projected_density = PDOS.get_site_dos(\
                    atom_list=indices, orbital_list=orbital_list\
                                 ,sum_density=sum_density, sum_spin=sum_spin)
            if sum_spin == False:
                projected_density[1::2,:] *= -1
            if settings == 'print':
                plt.figure()
            for count, density in enumerate(projected_density):
                plt.plot(density, PDOS.get_energies(), colors[count], zorder=zorder[count])
            plt.plot([np.min(projected_density), np.max(projected_density)]\
                     ,[PDOS.e_fermi, PDOS.e_fermi]\
                         ,'k--', zorder=1,linewidth=5) 
            plt.legend([i for i in orbitals]+ ['fermi level'])
            plt.xlabel('State density')
            plt.ylabel('Energy [eV]')
            plt.title(title)
            if settings == 'print':
                plt.show()
            
        #plot gas density
        plot_density(GAS_PDOS, gas_indices, 'Gas')
        #plot adsorbate density
        plot_density(REFERENCE_PDOS, adsorbate_indices, 'Adsorbate')
        #plot adsorption-site density
        plot_density(REFERENCE_PDOS, site_indices, 'Adsorbation-site')
