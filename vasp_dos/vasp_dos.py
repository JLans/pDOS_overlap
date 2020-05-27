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
    """ Get default paths to package data.
    
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

def get_all_VASP_files(directory):
    """ Get all DOSCAR and CONTCAR file paths
    
    Parameters
    ----------
    directory : str
        root directory to look for DOSCAR files
    
    Returns
    -------
    example_data_path : str
        path to example VASP data
    """
    
    DOSCAR_directories = [os.path.join(r,subdirectory) for r,d,f in os.walk(directory) \
              for subdirectory in d \
              if 'DOSCAR' in os.listdir(os.path.join(r,subdirectory))]
    DOSCAR_files = [os.path.join(d,'DOSCAR') for d in DOSCAR_directories]
    CONTCAR_files = [os.path.join(d,'CONTCAR') for d in DOSCAR_directories]
    return DOSCAR_files, CONTCAR_files

def get_band_center(energies, densities, max_energy=None, axis=-1):
        """ Get band center given energies and densities
        
        Parameters
        ----------
        energies : numpy.ndarray
            discretized orbital energies
        
        densities : numpy.ndarray
            projected state densities
            
        max_energy : float
            cutoff energy (often the fermi energy)
            
        axis : int
            axis of densities on which to integrate
        
        Returns
        -------
        band_center : float or numpy.ndarray
            center of the band(s) up to max_energy
            
        Notes
        -----
        trapezoidal rule is better for narrow gaussian peaks and for "rough" functions
        https://doi.org/10.1016/j.chemolab.2018.06.001
        http://emis.icm.edu.pl/journals/JIPAM/v3n4/031_02.html
        
        """
        
        if len(densities.shape) == 1:
            densities = np.array([densities.copy()])
        if max_energy is None:
            idx_stop = len(energies)
        else:
            idx_stop = (np.abs(energies-max_energy)).argmin()+1
        Integrated_Energy = np.trapz(densities[:,0:idx_stop]*energies[0:idx_stop], energies[0:idx_stop], axis=axis)
        Integrated_Filling = np.trapz(densities[:,0:idx_stop], energies[0:idx_stop], axis=axis)
        band_center = Integrated_Energy/Integrated_Filling
        return band_center

class VASP_DOS:
    """Class for extracting projected density of states from VASP"""
    def __init__(self, file_name="DOSCAR"):
        """ 
        Parameters
        ----------
        file_name : str
            full DOSCAR file location
        
        Attributes
        ----------
        natoms : int
            number of atoms in the system
            
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
            
        m_projected : bool
            indicates if projected density is orbital resolved
            
        orbital_dictionary: dict
            dictionary that maps resolved sublevels/orbitals to indices
        """
        
        natoms, emax, emin, ndos, e_fermi, is_spin, m_projected\
            , orbital_dictionary = self._read_doscar(file_name=file_name)
        self.natoms = natoms
        self.emax = emax
        self.emin = emin
        self.ndos = ndos
        self.e_fermi = e_fermi
        self.is_spin = is_spin
        self.m_projected = m_projected
        self.orbital_dictionary = orbital_dictionary
        
    def get_band_center(self, atom, orbital_list, sum_density=False, max_energy=None, axis=-1):
        """ Get band center for a given atom and list of orbitals
        
        Parameters
        ----------
        atom : int
            Atom index
            
        orbital_list : list[str]
            Which orbitals to return
            
        sum_density : bool
            if a sub-level is provided instead of an orbital, sum_density
            indicates if the individual sub-level densities should be summed
            
        max_energy : float
            cutoff energy (often the fermi energy)
            
        axis : int
            axis of densities on which to integrate
        
        Returns
        -------
        band_center : float or numpy.ndarray
            center of the band(s) up to max_energy

        """
        
        energies = self.get_energies()
        get_site_dos = self.get_site_dos
        orbitals, density = get_site_dos(atom,orbital_list, sum_density=sum_density)
        band_center = get_band_center(energies, density, max_energy=max_energy, axis=-1)
        return band_center
        
    def get_energies(self):
        """ method for obtaining energy levels
        
        Returns
        -------
        energies : numpy.ndarray
            1-D array of energies
            
        """
        
        energies = self._total_dos[0,:].copy()
        return energies
    
    def get_total_dos(self, sum_spin=False):
        """ method for obtaining total density of states
        
        Returns
        -------
        total_dos : numpy.ndarray
            1-D or 2-D array of state densities   
        """
        
        if self._total_dos.shape[0] == 3:
            total_dos = self._total_dos[1, :]
        elif self._total_dos.shape[0] == 5:
            if sum_spin == True:
                total_dos = self._total_dos[1:3, :].sum(axis=0)
            else:
                total_dos = self._total_dos[1:3, :]
                
        return total_dos

    def get_integrated_dos(self, sum_spin=False):
        """ method for obtaining total integrated density of states
        
        Returns
        -------
        integrated_dos : numpy.ndarray
            1-D or 2-D array of state integrated densities  
        """
        
        if self._total_dos.shape[0] == 3:
            integrated_dos = self._total_dos[2, :]
        elif self._total_dos.shape[0] == 5:
            if sum_spin == True:
                integrated_dos = self._total_dos[3:5, :].sum(axis=0)
            else:
                integrated_dos = self._total_dos[3:5, :].sum(axis=0)
        return integrated_dos
        
    def get_site_dos(self, atom, orbital_list, sum_density = False):
        """Return an NDOSxM array with dos for the chosen atom and orbital(s).
        
        Parameters
        ----------
        atom : int
            Atom index
            
        orbital_list : list[str]
            Which orbitals to return
            
        sum_density : bool
            if a sub-level is provided instead of an orbital, sum_density
            indicates if the individual sub-level densities should be summed
        
        Returns
        -------
        new_orbital_list : list[str]
            If sum_density is True, new_orbital_list is the same as
            orbital_list. If sum_density is False, new_orbital_list is resolved
            by both orbital and spin (if available)
            
        projected_density : np.array
            Array of shape (len(new_orbital_list), ndos)

        """

        # Integer indexing for orbitals starts from 1 in the _site_dos array
        # since the 0th column contains the energies
        orbital_dictionary = self.orbital_dictionary
        is_spin = self.is_spin
        m_projected = self.m_projected
        _site_dos = self._site_dos
        ndos = self.ndos
        def get_orbitals(orbital):
            #case where spin polarization is false and m level is resloved
            if is_spin == False and m_projected == True:
                if orbital == 'p':
                    orbitals = ['py', 'pz', 'px']
                elif orbital == 'd':
                    orbitals = ['dxy', 'dyz', 'dz2', 'dxz', 'dx2-y2']
                elif orbital == 'f':
                    orbitals = ['fy(3x2-y2)', 'fxyz', 'fyz2', 'fz3', 'fxz2',
                              'fz(x2-y2)', 'fx(x2-3y2)']
                else:
                    orbitals = [orbital]
            #case where spin polarization is true and m level is resloved
            elif is_spin == True and m_projected == True:
                #case where neither spin or m-level are provided
                if orbital == 's':
                    orbitals = ['s+','s-']
                elif orbital == 'p':
                    orbitals = ['py+', 'py-', 'pz+', 'pz-', 'px+', 'px-']
                elif orbital == 'd':
                    orbitals = ['dxy+', 'dxy-', 'dyz+', 'dyz-', 'dz2+', 'dz2-',
                              'dxz+', 'dxz-', 'dx2-y2+', 'dx2-y2-']
                elif orbital == 'f':
                    orbitals = ['fy(3x2-y2)+', 'fy(3x2-y2)-','fxyz+', 'fxyz-',
                              'fyz2+', 'fyz2-', 'fz3+', 'fz3-', 'fxz2+',
                              'fxz2-', 'fz(x2-y2)+', 'fz(x2-y2)-',
                              'fx(x2-3y2)+', 'fx(x2-3y2)-']
                #case where spin is provided but m-level is not
                #up spin
                elif orbital == 'p+':
                    orbitals = ['py+', 'pz+', 'px+']
                elif orbital == 'd+':
                    orbitals = ['dxy+', 'dyz+', 'dz2+', 'dxz+', 'dx2-y2+']
                elif orbital == 'f+':
                    orbitals = ['fy(3x2-y2)+', 'fxyz+', 'fyz2+', 'fz3+', 'fxz2+',
                              'fz(x2-y2)+', 'fx(x2-3y2)+']
                # down spin
                elif orbital == 'p-':
                    orbitals = ['py-', 'pz-', 'px-']
                elif orbital == 'd-':
                    orbitals = ['dxy-', 'dyz-', 'dz2-', 'dxz-', 'dx2-y2-']
                elif orbital == 'f-':
                    orbitals = ['fy(3x2-y2)-', 'fxyz-', 'fyz2-', 'fz3-', 'fxz2-',
                              'fz(x2-y2)-', 'fx(x2-3y2)-']
                else:
                    orbitals = [orbital]
            else:
                orbitals = [orbital]
            return orbitals
        
        if sum_density == True:
            projected_density = np.zeros((len(orbital_list),ndos))
            for count, orbital in enumerate(orbital_list):
                new_orbital_list = get_orbitals(orbital)
                indices = [orbital_dictionary[key] for key in new_orbital_list]
                projected_density[count] = _site_dos[atom, indices,:].sum(axis=0)
            #return the list of orbitals provided
            new_orbital_list = orbital_list
        else:
            #get complete set of orbitals
            new_orbital_list = []
            for orbital in orbital_list:
                new_orbital_list += get_orbitals(orbital)
            indices = [orbital_dictionary[key] for key in new_orbital_list]
            projected_density = _site_dos[atom, indices, :]
        return new_orbital_list, projected_density
        
    def _read_doscar(self, file_name="DOSCAR"):
        """Read VASP DOSCAR and extract projected densities
        
        Parameters
        ----------
        file_name : str
            file location of the DOSCAR file
            
        Attributes
        ----------
        _total_dos : numpy.ndarray
            numpy array that contains the energy of the orbitals and the
            total projected and integrated density
            
        _site_dos : numpy.ndarray
            numpy array that contains the energy of the orbitals and the
            site and orbital projected density of states. Only available if a
            site projected calculation was performed.
            
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
            
        m_projected : bool
            indicates if projected density is orbital resolved
            
        orbital_dictionary : dict
            dictionary that maps resolved sublevels/orbitals to indices

        """
        
        #Accepts a file and reads through to get the density of states
        def get_dos(f, ndos):
            #get first line of density
            line = f.readline().split()
            dos = np.zeros((ndos, len(line)))
            dos[0] = np.array(line)
            for nd in range(1, ndos):
                line = f.readline().split()
                dos[nd] = np.array([float(x) for x in line])
            return dos.T
            
        f = open(file_name)
        natoms = int(f.readline().split()[0])
        [f.readline() for lines in range(4)]  # Skip next 4 lines.
        # First we have a block with total and total integrated DOS
        descriptive_line = f.readline().split()
        emax = float(descriptive_line[0])
        emin = float(descriptive_line[1])
        ndos = int(descriptive_line[2])
        e_fermi = float(descriptive_line[3])
        _total_dos = get_dos(f,ndos)
        # Next we have one block per atom, if DOSCAR contains the stuff
        # necessary for generating site-projected DOS
        dos = []
        for na in range(natoms):
            line = f.readline()
            if line == '':
                # No site-projected DOS
                break
            #ndos = int(line.split()[2]) ndos does not change
            pdos = get_dos(f,ndos)
            dos.append(pdos)
        _site_dos = np.array(dos)
        
        
        # Integer indexing for orbitals starts from 1 in the _site_dos array
        # since the 0th column contains the energies
        norbs = _site_dos.shape[1] - 1
        if norbs == 3:
            is_spin = False
            m_projected = False
            orbitals = {'s': 1, 'p': 2, 'd': 3}
        elif norbs == 4:
            is_spin = False
            m_projected = False
            orbitals = {'s': 1, 'p': 2, 'd': 3, 'f': 4}
        elif norbs == 6:
            is_spin = True
            m_projected = False
            orbitals = {'s+': 1, 's-': 2, 'p+': 3, 'p-': 4, 'd+': 5, 'd-': 6}
        elif norbs == 8:
            is_spin = True
            m_projected = False
            orbitals = {
                's+': 1,
                's-': 2,
                'p+': 3,
                'p-': 4,
                'd+': 5,
                'd-': 6,
                'f+': 7,
                'f-': 8}
        elif norbs == 9:
            is_spin = False
            m_projected = True
            orbitals = {'s': 1, 'py': 2, 'pz': 3, 'px': 4,
                    'dxy': 5, 'dyz': 6, 'dz2': 7, 'dxz': 8, 'dx2-y2': 9}
        elif norbs == 16:
            is_spin = False
            m_projected = True
            orbitals = {
                's': 1,
                'py': 2,
                'pz': 3,
                'px': 4,
                'dxy': 5,
                'dyz': 6,
                'dz2': 7,
                'dxz': 8,
                'dx2': 9,
                'fy(3x2-y2)': 10,
                'fxyz': 11,
                'fyz2': 12,
                'fz3': 13,
                'fxz2': 14,
                'fz(x2-y2)': 15,
                'fx(x2-3y2)': 16}
        elif norbs == 18:
            is_spin = True
            m_projected = True
            orbitals = {
                's+': 1,
                's-': 2,
                'py+': 3,
                'py-': 4,
                'pz+': 5,
                'pz-': 6,
                'px+': 7,
                'px-': 8,
                'dxy+': 9,
                'dxy-': 10,
                'dyz+': 11,
                'dyz-': 12,
                'dz2+': 13,
                'dz2-': 14,
                'dxz+': 15,
                'dxz-': 16,
                'dx2-y2+': 17,
                'dx2-y2-': 18}
        elif norbs == 32:
            is_spin = True
            m_projected = True
            orbitals = {
                's+': 1,
                's-': 2,
                'py+': 3,
                'py-': 4,
                'pz+': 5,
                'pz-': 6,
                'px+': 7,
                'px-': 8,
                'dxy+': 9,
                'dxy-': 10,
                'dyz+': 11,
                'dyz-': 12,
                'dz2+': 13,
                'dz2-': 14,
                'dxz+': 15,
                'dxz-': 16,
                'dx2-y2+': 17,
                'dx2-y2-': 18,
                'fy(3x2-y2)+': 19,
                'fy(3x2-y2)-': 20,
                'fxyz+': 21,
                'fxyz-': 22,
                'fyz2+': 23,
                'fyz2-': 24,
                'fz3+': 25,
                'fz3-': 26,
                'fxz2+': 27,
                'fxz2-': 28,
                'fz(x2-y2)+': 29,
                'fz(x2-y2)-': 30,
                'fx(x2-3y2)+': 31,
                'fx(x2-3y2)-': 32}
        
        self._total_dos = _total_dos
        self._site_dos = _site_dos
        return natoms, emax, emin, ndos, e_fermi, is_spin, m_projected, orbitals
        
        