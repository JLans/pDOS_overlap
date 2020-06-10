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
    data_path = pkg_resources.resource_filename(__name__, 'data/example_data')
    return data_path

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

def get_band_center(energies, densities, max_energy=None, min_energy=None, axis=-1):
        """ Get band center given energies and densities
        
        Parameters
        ----------
        energies : numpy.ndarray
            discretized orbital energies
        
        densities : numpy.ndarray
            projected state densities
            
        max_energy : float
            cutoff energy (often the fermi energy)
            
        min_energy : float
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
        if min_energy is None:
            idx_start = 0
        else:
            idx_start = (np.abs(energies-min_energy)).argmin()
        Integrated_Energy = np.trapz(densities[:,idx_start:idx_stop]\
                                     * energies[idx_start:idx_stop]\
                                     , energies[idx_start:idx_stop], axis=axis)
        Integrated_Filling = np.trapz(densities[:,idx_start:idx_stop]\
                                      , energies[idx_start:idx_stop], axis=axis)
        band_center = Integrated_Energy/Integrated_Filling
        return band_center
    
def get_band_width(energies, densities, fraction=0):
        """ Get band width given energies and densities
        
        Parameters
        ----------
        energies : numpy.ndarray
            discretized orbital energies
        
        densities : numpy.ndarray
            projected state densities
            
        fraction : float
            maximum projected density considered part of a band
        
        Returns
        -------
        band_width : float or numpy.ndarray
            width of the band(s)
        """
        if len(densities.shape) == 1:
            densities = np.array([densities.copy()])
        index_values = np.arange(densities.shape[-1])
        min_index = np.zeros(densities.shape[0])
        max_index = np.zeros(densities.shape[0])
        for count, density in enumerate(densities):
            min_index[count] = np.min(index_values[density > fraction * density.max()])
            max_index[count] = np.max(index_values[density > fraction * density.max()])
        band_width = energies[max_index.astype(int)] - energies[min_index.astype(int)]
        return band_width
    
def get_center_width(energies, densities, energy):
        """ Get width between to band centers given a division
        
        Parameters
        ----------
        energies : numpy.ndarray
            discretized orbital energies
        
        densities : numpy.ndarray
            projected state densities
            
        energy : int
            energy that divides band centers
        
        Returns
        -------
        center_width : float or numpy.ndarray
            width of two band centers separated by some ennergy
            
        """
        center_lower = get_band_center(energies, densities, max_energy=energy)
        center_upper = get_band_center(energies, densities, min_energy=energy)
        center_width = center_upper - center_lower
        return center_width
    
def get_filling(energies, densities, max_energy=None, min_energy=None, axis=-1):
        """ Get band center given energies and densities
        
        Parameters
        ----------
        energies : numpy.ndarray
            discretized orbital energies
        
        densities : numpy.ndarray
            projected state densities
            
        max_energy : float
            cutoff energy (often the fermi energy)
            
        min_energy : float
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
        if min_energy is None:
            idx_start = 0
        else:
            idx_start = (np.abs(energies-min_energy)).argmin()
        Integrated_Filling = np.trapz(densities[:,idx_start:idx_stop]\
                                      , energies[idx_start:idx_stop], axis=axis)
        return Integrated_Filling

class VASP_DOS:
    """Class for extracting projected density of states from VASP"""
    def __init__(self, file_name="DOSCAR",no_negatives=True):
        """ 
        Parameters
        ----------
        file_name : str
            full DOSCAR file location
        
        no_negatives : bool
            Indicates wheather negative occupances will be converted to zero.
            Negative occupances can occur if Methfessel-Paxton is used.
        
        Attributes
        ----------
        no_negatives : bool
            Indicates wheather negative occupances will be converted to zero.
            Negative occupances can occur if Methfessel-Paxton is used.
        
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
            , orbital_dictionary = self._read_doscar(file_name=file_name, no_negatives=no_negatives)
        self.no_negatives = no_negatives
        self.natoms = natoms
        self.emax = emax
        self.emin = emin
        self.ndos = ndos
        self.e_fermi = e_fermi
        self.is_spin = is_spin
        self.m_projected = m_projected
        self.orbital_dictionary = orbital_dictionary
        
    def get_band_center(self, atom_indices, orbital_list, sum_density=False, sum_spin=True\
                        , max_energy=None, min_energy=None, axis=-1):
        """ Get band center for a given atom and list of orbitals
        
        Parameters
        ----------
        atom_indices : list[int]
            list of atom indices
            
        orbital_list : list[str]
            Which orbitals to return
            
        sum_density : bool
            if a sub-level is provided instead of an orbital, sum_density
            indicates if the individual sub-level densities should be summed
            
        sum_spin : bool
            different spin densities are summed.
            
        max_energy : float
            cutoff energy (often the fermi energy)
        
        min_energy : float
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
        try:
            len(atom_indices)
        except:
            atom_indices = [atom_indices]
        orbitals, density = get_site_dos(atom_indices,orbital_list\
                                         , sum_density=sum_density\
                                         , sum_spin=sum_spin)
        band_center = get_band_center(energies, density, max_energy=max_energy\
                                      , min_energy = min_energy, axis=-1)
        return band_center
    
    def get_filling(self, atom_indices, orbital_list, sum_density=False, sum_spin=True\
                        , max_energy=None, min_energy=None, axis=-1):
        """ Get band center for a given atom and list of orbitals
        
        Parameters
        ----------
        atom_indices : list[int]
            list of atom indices
            
        orbital_list : list[str]
            Which orbitals to return
            
        sum_density : bool
            if a sub-level is provided instead of an orbital, sum_density
            indicates if the individual sub-level densities should be summed
            
        sum_spin : bool
            different spin densities are summed.
            
        max_energy : float
            cutoff energy (often the fermi energy)
        
        min_energy : float
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
        try:
            len(atom_indices)
        except:
            atom_indices = [atom_indices]
        orbitals, density = get_site_dos(atom_indices,orbital_list\
                                         , sum_density=sum_density\
                                         , sum_spin=sum_spin)
        filling = get_filling(energies, density, max_energy=max_energy\
                                      , min_energy = min_energy, axis=-1)
        return filling

    def get_band_width(self, atom_indices, orbital_list, sum_density=False\
                       , sum_spin=True, fraction=0):
        """ Get band width given energies and densities
        
        Parameters
        ----------
        atom_indices : list[int]
            list of atom indices
            
        orbital_list : list[str]
            Which orbitals to return
            
        sum_density : bool
            if a sub-level is provided instead of an orbital, sum_density
            indicates if the individual sub-level densities should be summed
            
        sum_spin : bool
            different spin densities are summed.
            
        fraction : float
            maximum projected density considered part of a band
        
        Returns
        -------
        band_width : float or numpy.ndarray
            width of the band(s)
        """
        energies = self.get_energies()
        get_site_dos = self.get_site_dos
        try:
            len(atom_indices)
        except:
            atom_indices = [atom_indices]
        orbitals, densities = get_site_dos(atom_indices,orbital_list\
                                         , sum_density=sum_density\
                                         , sum_spin=sum_spin)
        band_width = get_band_width(energies, densities, fraction)
        return band_width
    
    def get_center_width(self, energy, atom_indices, orbital_list, sum_density=False\
                       , sum_spin=True):
        """ Get width between to band centers given a division
        
        Parameters
        ----------            
        energy : int
            energy that divides band centers
            
        atom_indices : list[int]
            list of atom indices
            
        orbital_list : list[str]
            Which orbitals to return
            
        sum_density : bool
            if a sub-level is provided instead of an orbital, sum_density
            indicates if the individual sub-level densities should be summed
            
        sum_spin : bool
            different spin densities are summed.
        
        Returns
        -------
        center_width : float or numpy.ndarray
            width of two band centers separated by some ennergy
            
        """
        energies = self.get_energies()
        get_site_dos = self.get_site_dos
        try:
            len(atom_indices)
        except:
            atom_indices = [atom_indices]
        orbitals, densities = get_site_dos(atom_indices,orbital_list\
                                         , sum_density=sum_density\
                                         , sum_spin=sum_spin)
        
        center_width = get_center_width(energies, densities, energy)
        return center_width
        
    def get_energies(self):
        """ method for obtaining energy levels
        
        Returns
        -------
        energies : numpy.ndarray
            1-D array of energies    
        """
        energies = self._total_dos[0,:].copy()
        return energies
    
    def get_total_dos(self, sum_spin=True):
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

    def get_integrated_dos(self, sum_spin=True):
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
                integrated_dos = self._total_dos[3:5, :]
        return integrated_dos
        
    def get_site_dos(self, atom_indices, orbital_list=[], sum_density=False, sum_spin=True):
        """Return an NDOSxM array with dos for the chosen atom and orbital(s).
        
        Parameters
        ----------
        atom_list : list[int]
            List of atom index
            
        orbital_list : list[str]
            Which orbitals to return
            
        sum_density : bool
            if a sub-level is provided instead of an orbital, sum_density
            indicates if the individual orbital densities should be summed
            
        sum_spin : bool
            different spin densities are summed.
        
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
        try:
            len(atom_indices)
        except:
            atom_indices = [atom_indices]
        if len(orbital_list) == 0:
            orbital_list = list(orbital_dictionary.keys())
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
                elif orbital == 's+':
                    orbitals = ['s+']
                elif orbital == 's-':
                    orbitals = ['s-']
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
                elif '+' not in orbital and '-' not in orbital:
                    orbitals = [orbital + '+', orbital + '-']
                else:
                    orbitals = [orbital]
            else:
                orbitals = [orbital]
            return orbitals
        #force list of atom_list is an int
        if sum_density == True:
            projected_density = np.zeros((len(orbital_list), ndos))
            for atom in atom_indices:
                for count, orbital in enumerate(orbital_list):
                    new_orbital_list = get_orbitals(orbital)
                    indices = [orbital_dictionary[key] for key in new_orbital_list]
                    projected_density[count]+= _site_dos[atom, indices,:].sum(axis=0)
            #return the list of orbitals provided
            new_orbital_list = orbital_list
        else:
            #get complete set of orbitals
            new_orbital_list = []
            for orbital in orbital_list:
                new_orbital_list += get_orbitals(orbital)
            projected_density = np.zeros((len(new_orbital_list), ndos))
            indices = [orbital_dictionary[key] for key in new_orbital_list]
            for atom in atom_indices:
                projected_density += _site_dos[atom, indices, :]
            
        if is_spin == True and sum_spin == True and sum_density == False:
            projected_density = projected_density[0::2,:] + projected_density[1::2,:]
            new_orbital_list = [key.rstrip('+') for key in new_orbital_list[::2]]
        return new_orbital_list, projected_density
        
    def _read_doscar(self, file_name="DOSCAR", no_negatives=True):
        """Read VASP DOSCAR and extract projected densities
        
        Parameters
        ----------
        file_name : str
            file location of the DOSCAR file
            
        no_negatives : bool
            Indicates wheather negative occupances will be converted to zero.
            Negative occupances can occur if Methfessel-Paxton is used.
            
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
        if no_negatives == True:
            _total_dos[1:][_total_dos[1:][...] < 0] = 0
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
        if no_negatives == True:
            _site_dos[:,1:,:][_site_dos[:,1:,:][...] < 0] = 0 
        
        
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
        
        