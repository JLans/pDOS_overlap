# -*- coding: utf-8 -*-
"""
Created on Thu Mar 02 17:26:23 2017

@author: lansford
"""

"""Vibrational modes."""
from ase.data import covalent_radii as CR
from ase.io import read
import numpy as np

def get_geometric_data(CONTCAR, cutoff=1.25, crystal_type='fcc'):
    """ Obtain important geometric data for all atoms
        
    Parameters
    ----------
    CONTCAR : str
        Location to VASP CONTCAR file for a gas molecule.
    
    cutoff : float
            Fudge factor multiplied by covalent radii used to determine
            connectivity. If distance between atoms is less than the summ of
            their radii times cutoff they are considered to be coordinated
    
    Returns
    -------
    indices : list[int]
        Atom indices
        
    GCNs : list[float]
        GCN values of all atoms
        
    atom_types : list[str]
        Indicates if atoms are of type 'surface' or 'bulk'
    """
    if crystal_type == 'fcc':
        bulk_CN = 12
    elif crystal_type == 'bcc':
        bulk_CN = 14
    else:
        bulk_CN = crystal_type
    
    nanoparticle = read(CONTCAR)
    CN = Coordination(nanoparticle,cutoff=cutoff)
    # read and return densityofstates object
    indices = np.arange(len(nanoparticle))
    atom_types = []
    GCNs = np.zeros_like(indices)
    for atom_index in indices:
        GCNs[atom_index] = CN.get_gcn([atom_index])
        if CN.cn[atom_index] < bulk_CN:
            atom_types.append('surface')
        else:
            atom_types.append('bulk')
    return indices, GCNs, np.array(atom_types)

class Coordination:
    """Class for calculating coordination and generalized coordination number."""

    def __init__(self, atoms, exclude=[], cutoff=1.25, cutoff_type = 'percent'):
        """
        Parameters
        ----------
        atoms : atoms
        	ASE atoms type.
            
        exclude : list of int
        	Indices of atoms to exclude in tabulating the coordination of
            each atom
            
        cutoff : float
            Fudge factor multiplied by covalent radii used to determine
            connectivity. If distance between atoms is less than the summ of
            their radii times cutoff they are considered to be coordinated
            
        cutoff_type : str
            Can be 'percent' or 'absolute'. If Absolute then the cutoff is
            considered a distance and replaces the use of Van der Waals radii
            entirely. If multiple atom types are to be used then 'perecent'
            is best.
        	
        Attributes
        ----------
        atoms : ase.Atoms
        	ASE atoms type.
            
        exclude : list of int
        	Indices of atoms to exclude in tabulating the coordination of
            each atom
            
        cutoff : float
            Fudge factor multiplied by Van der Waals radii used to determine
            connectivity. If distance between atoms is less than the summ of
            their radii times cutoff they are considered to be coordinated
            
        cutoff_type : str
            Can be 'percent' or 'absolute'. If Absolute then the cutoff is
            considered a distance and replaces the use of Van der Waals radii
            entirely. If multiple atom types are to be used then 'perecent'
            is best.   
        """
        self.atoms = atoms
        self.exclude = exclude
        self.cutoff = cutoff
        self.cutoff_type = cutoff_type
        self.read()
        
    def read(self):
        """Returns an array of coordination numbers and an array of existing bonds determined by
        distance and covalent radii.  By default a bond is defined as 120% of the combined radii
        or less. This can be changed by setting 'cutoff' to a float representing a 
        factor to multiple by (default = 1.2).
        If 'exclude' is set to an array,  these atomic numbers with be unable to form bonds.
        This only excludes them from being counted from other atoms,  the coordination
        numbers for these atoms will still be calculated,  but will be unable to form
        bonds to other excluded atomic numbers.
        
        Attributes
        ----------
        cn : list of int
        	list of coordination numbers for each atom.
            
        bonded : list of list
        	List of indices of atoms bonded to each atom 
        """
        atoms = self.atoms
        cutoff = self.cutoff
        cutoff_type = self.cutoff_type
        exclude = self.exclude
        # Get all the distances
        distances = np.divide(atoms.get_all_distances(mic=True), cutoff)
        
        # Array of indices of bonded atoms.  len(bonded[x]) == cn[x]
        bonded = []
        indices = list(range(len(atoms)))
        
        # Atomic Numbers
        numbers = atoms.numbers
        # Coordination Numbers for each atom
        cn = []
        
        cr = np.take(CR, numbers)
        if cutoff_type == 'absolute':
            for i in indices:
                cr[i] = cutoff/2.
            distances = atoms.get_all_distances(mic=True)
                
        for i in indices:
            bondedi = []
            for ii in indices:
                # Skip if measuring the same atom
                if i == ii or ii in exclude:
                    continue
                if (cr[i] + cr[ii]) >= distances[i,ii]:
                    bondedi.append(ii)
            # Add this atoms bonds to the bonded list
            bonded.append(bondedi)
        for i in bonded:
            cn.append(len(i))
        self.cn = cn
        self.bonded = bonded
        
    def get_coordination_numbers(self):
        """Implements returns coordination values
        
        Returns
        -------
        cn : list of int
        	list of coordination numbers for each atom.
        """
        return self.cn
    
    def get_bonded(self):
        """returns bonded list
        
        Returns
        -------
        bonded : list of list
        	List of indices of atoms bonded to each atom
        """
        return self.bonded
    
    def get_gcn(self,site=[], surface_type="fcc"):
        """Returns the generalized coordination number of the site given.  To define
        A site,  just give a list containing the indices of the atoms forming the site.
        The calculated coordination numbers and bonding needs to be supplied, and the
        surface type also needs to be changed if the surface is not represented by bulk fcc.
        
        Parameters
        ----------
        
        site : list of int
            indices for which to calculate the GCN values
            
        surface_type : str
            Indicates natural bulk of the material which influences
            the maximum coordination environment of a site
        
        Returns
        -------
        gcn : float
        	Generalized coordination values for the desired site    
        """
        bonded = self.bonded
        cn = self.cn
        # Define the types of bulk accepted
        gcn_bulk = {"fcc": [12., 18., 22., 26., 26.], "bcc": [14., 22., 28., 32.]}
        sum_gcn = 0.
        if len(site) == 0:
            return 0
        counted = []
        for i in site:
            counted.append(i)
        for i in site:
            for ii in bonded[i]:
                if ii in counted:
                    continue
                counted.append(ii)
                sum_gcn += cn[ii]
        gcn = sum_gcn / gcn_bulk[surface_type][len(site) - 1]
        return gcn