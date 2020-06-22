# -*- coding: utf-8 -*-
"""
Created on Tue May 19 22:28:09 2020

@author: lansf
"""

from __future__ import absolute_import, division, print_function
import os
import pkg_resources
import numpy as np
from ase.io import read
from .vasp_dos import VASP_DOS
from .coordination import Coordination
import itertools

def write_lobsterin(directory='.', adsorbate_atoms=['C','H','O','N']\
                    ,basisSet='pbeVaspFit2015'\
                    , basisfunctions={'C':'2s 2p', 'O':'2s 2p'\
                                      , 'N':'2s 2p', 'H':'1s'\
                                      ,'Pt':'5d 6s'}\
                    , gaussianSmearingWidth=0.003):
    """ Write a lobster input file
    
    Parameters
    ----------
    directory : str
        directory to write the lobsterin file. Must also contain a CONTCAR and
        DOSCAR file
        
    adsorbate_atoms : list[str or int]
        adsorbate atom symbols or indices
        
    basisSet : str
        basis used to projected density onto orbitals
        
    basisfunctions : dict
        dictionary of atomic symbols and corresponding basis functions
        
    gaussianSmearingWidth : float
        float if Gaussian smearing is used. If tetrahedron method is used then
        set to None
    """
    if type(adsorbate_atoms[0]) == str:
        CONTCAR = read(os.path.join(directory,'CONTCAR'))
        all_symbols = CONTCAR.get_chemical_symbols()
        adsorbate_indices = np.arange(len(all_symbols))[np.isin(all_symbols,adsorbate_atoms)]
    else:
        adsorbate_indices = adsorbate_atoms
    CN = Coordination(CONTCAR,exclude=adsorbate_indices,cutoff=1.25)
    bonded = CN.get_bonded()
    site_indices = []
    for i in adsorbate_indices:
        site_indices += bonded[i]
    site_indices = list(set(site_indices))
    atom_pairs = []
    for site_index in site_indices:
        for adsorbate_index in adsorbate_indices:
            atom_pairs.append((site_index, adsorbate_index))
    atom_pairs += list(itertools.combinations(adsorbate_indices,2))
    atom_pairs = np.array(atom_pairs) + 1
    DOSCAR = os.path.join(directory,'DOSCAR')
    DOS = VASP_DOS(DOSCAR)
    COHPstartEnergy = DOS.emin - DOS.e_fermi
    COHPendEnergy = DOS.emax - DOS.e_fermi
    COHPSteps = int(DOS.ndos - 1)
    basisSet = basisSet
    basisfunctions = basisfunctions
    gaussianSmearingWidth = gaussianSmearingWidth
    file_name = os.path.join(directory,'lobsterin')
    file = open(file_name,'w')
    file.write('COHPstartEnergy ' + str(COHPstartEnergy))
    file.write('\n')
    file.write('COHPendEnergy ' + str(COHPendEnergy))
    file.write('\n')
    file.write('COHPSteps ' + str(COHPSteps))
    file.write('\n')
    file.write('basisSet ' + str(basisSet))
    file.write('\n')
    if gaussianSmearingWidth is not None:
        file.write('gaussianSmearingWidth ' + str(gaussianSmearingWidth))
        file.write('\n')
    for key in basisfunctions.keys():
        file.write('basisfunctions ' + key + ' ' + basisfunctions[key])
        file.write('\n')
    for pair in atom_pairs:
        file.write('cohpbetween ' + 'atom ' + str(pair[0]) + ' atom ' + str(pair[1]))
        file.write('\n')
    file.close()

def get_example_data():
    """ Get default path to experimental crystal ovelap populaiton data
    
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

def get_bonding_fraction(orbital_indices, dos_energies, pcoop, pcoop_energies\
                             , set_antibonding_zero=False, emax=float('inf')):
        """ method for obtaining bonding fraction of dos or dos-like array
        
        Parameters
        ----------
        orbital_indices : list[list]
            list of orbital indices for each molecular orbital
            
        dos_energies : numpy.ndarray
            energies at which dos is calculated
            
        pcoop : numpy.ndarray
             projected orbital overlap populations
        
        pcoop_energies : numpy.ndarray
            energies at which the pcoop is calculated
            
        set_antibonding_zero : bool
            if true, set antibonding populations to zero to look at total
            instead of net bonding characteristics
        
        Returns
        -------
        dos_bonding : numpy.ndarray
            1-D array of the relative bonding for either dos-like array
        """
        if set_antibonding_zero == True:
            pcoop[pcoop[...] < 0] = 0
        pcoop = np.interp(dos_energies, pcoop_energies,pcoop)
        pcoop[dos_energies[...] > emax] = 0
        dos_bonding = []
        for i in range(len(orbital_indices)):
            dos_bonding.append(np.trapz(pcoop[orbital_indices[i]],dos_energies[orbital_indices[i]]))
        return dos_bonding

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
        """ obtain average overlap population for all atom pairs
        
        Parameters
        ----------
        sum_spin : bool
            indicates whether data of different spins should be summed
        
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
        """obtain average integrated overlap population for all atom pairs
        
        Parameters
        ----------
        sum_spin : bool
            indicates whether data of different spins should be summed
        
        Returns
        -------
        average_int_pcoop : numpy.ndarray
            1-D or 2-D array of average integrated pcoop for all interactions
        """
        if self.is_spin == True:
            if sum_spin == False:
                average_int_pcoop = self._pcoop[[2, 2 * self.num_interactions + 4], :]
            else:
                average_int_pcoop = self._pcoop[2, :] + self._pcoop[2 * self.num_interactions + 4, :]
        else:
            average_int_pcoop = self._pcoop[2, :]
        return average_int_pcoop
    
    def get_integrated_pcoop(self, interactions=[], sum_pcoop=False, sum_spin=True
                             , set_antibonding_zero=False):
        """ obtain integrated projected crystal orbital overlap populations
        
        Parameters
        ----------
        interactions : list
            indices of interactions for which to find the integrated pcoop
            
        sum_pcoop : bool
            indicates whether all pcoop should be summed
        
        sum_spin : bool
            indicates whether data of different spins should be summed
            
        set_antibonding_zero : bool
            if true, set antibonding populations to zero to look at total
            instead of net bonding characteristics
        
        Returns
        -------
        integrated_pcoop : numpy.ndarray
            1-D or 2-D array of integrated pcoop for all interactions
        """
        if len(interactions) == 0:
            interactions = list(range(self.num_interactions))
        if self.is_spin == True:
            spin_up = self._pcoop[4:2 * self.num_interactions + 3:2, :][interactions]
            spin_down = self._pcoop[2 * self.num_interactions + 6::2, :][interactions]
            if set_antibonding_zero == True:
                spin_up[spin_up[...] < 0] = 0
                spin_down[spin_down[...] < 0] = 0
            if sum_spin == True:
                integrated_pcoop = spin_up + spin_down
            else:
                integrated_pcoop = np.array([spin_up, spin_down])
        else:
            integrated_pcoop = self._pcoop[4::2, :][interactions]
            if set_antibonding_zero == True:
                integrated_pcoop[integrated_pcoop[...] < 0] = 0
        if sum_pcoop == True or len(interactions) == 1:
            axis = len(integrated_pcoop.shape) - 2
            integrated_pcoop = integrated_pcoop.sum(axis=axis)
        return integrated_pcoop
    
    def get_pcoop(self, interactions=[], sum_pcoop=False, sum_spin=True\
                  , set_antibonding_zero=False):
        """ method for obtaining projected crystal orbital overlap populations
        
        Parameters
        ----------
        interactions : list
            indices of interactions for which to find the integrated pcoop
            
        sum_pcoop : bool
            indicates whether all pcoop should be summed
        
        sum_spin : bool
            indicates whether data of different spins should be summed
            
        set_antibonding_zero : bool
            if true, set antibonding populations to zero to look at total
            instead of net bonding characteristics
        
        Returns
        -------
        pcoop : numpy.ndarray
            1-D or 2-D array of pcoop for all interactions
        """
        if len(interactions) == 0:
            interactions = list(range(self.num_interactions))
        if self.is_spin == True:
            spin_up = self._pcoop[3:2 * self.num_interactions + 3:2, :][interactions]
            spin_down = self._pcoop[2 * self.num_interactions + 5::2, :][interactions]
            if set_antibonding_zero == True:
                spin_up[spin_up[...] < 0] = 0
                spin_down[spin_down[...] < 0] = 0
            if sum_spin == True:
                pcoop = spin_up + spin_down
            else:
                pcoop = np.array([spin_up, spin_down])
        else:
            if set_antibonding_zero == True:
                pcoop[pcoop[...] < 0] = 0
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
    
    def get_bonding_fraction(self, orbital_indices, dos_energies\
                             , set_antibonding_zero=False\
                             , sum_pcoop=True, sum_spin=True, emax=float('inf')):
        """ method for obtaining bonding fraction of dos or dos-like array
        
        Parameters
        ----------
        orbital_indices : list[list]
            list of orbital indices for each molecular orbital
            
        dos_energies : numpy.ndarray
            energies at which dos is calculated
            
        set_antibonding_zero : bool
            if true, set antibonding populations to zero to look at total
            instead of net bonding characteristics

        sum_pcoop : bool
            indicates whether all pcoop should be summed
        
        sum_spin : bool
            indicates whether data of different spins should be summed
        
        Returns
        -------
        dos_bonding : numpy.ndarray
            1-D array of the relative bonding for either dos-like array
        """
        pcoop = self.get_pcoop(sum_pcoop=sum_pcoop, sum_spin=sum_spin\
                               , set_antibonding_zero=set_antibonding_zero)
        bonding_fraction = get_bonding_fraction(orbital_indices, dos_energies, pcoop\
                                        , self.get_energies()\
                                        , set_antibonding_zero=set_antibonding_zero\
                                        , emax = emax)
        return bonding_fraction
        
        