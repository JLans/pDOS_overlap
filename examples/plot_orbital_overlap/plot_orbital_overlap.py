"""
==============================================
Calculating orbital overlap using pdos_overlap
==============================================

This example shows how calculate the overlap of gas phase molecular orbitals
with an adsorbate and surface atom.

"""

import os
import numpy as np
from pdos_overlap import VASP_DOS
from pdos_overlap.plotting_tools import set_figure_settings
from pdos_overlap import get_adsorbate_indices
from pdos_overlap import PDOS_OVERLAP
from pdos_overlap.coordination import get_geometric_data

#######################################################################################
# Load DOSCAR file
# ----------------
#
# First we will, get the example data, load a DOSCAR file and use it to
# instantiate a VASP_DOS object.
gas = 'CO'
surface = 'Pt111'
set_figure_settings('paper')
np.set_printoptions(linewidth=100)
#These files are too large to store in the examples directory
example_path = r'C:\Users\lansf\Documents\Data\PROBE_PDOS\vasp_dos_files'
GAS_DOSCAR = os.path.join(example_path, gas + '/DOSCAR')
GAS_CONTCAR = os.path.join(example_path, gas + '/CONTCAR')
ADSORBATE_DOSCAR = os.path.join(example_path, gas + '+' + surface + '/DOSCAR')
ADSORBATE_CONTCAR = os.path.join(example_path, gas + '+' + surface + '/CONTCAR')

#######################################################################################
# Generate VASP_DOS objects
# -------------------------
#
# VASP_DOS objects for both the gas (vacuum) and the adsorbate+surface system
GAS_PDOS = VASP_DOS(GAS_DOSCAR)
REFERENCE_PDOS = VASP_DOS(ADSORBATE_DOSCAR)

#######################################################################################
# Get adsorbate and site indices and initialize PDOS_OVERLAP object
# -----------------------------------------------------------------
#
# This method utilizes two VASP_DOS objects, a gas and an adsorption system.
# It uses the adosorbtion system (REFERENCE_PDOS) to map gas molecular orbitals
# to adsorbate molecular orbitals. It then calculates the adsorption site
# atomic orbital energy overlaps with the adsorbate molecular orbital energies.
adsorbate_indices, site_indices = get_adsorbate_indices(GAS_CONTCAR\
                                                        , ADSORBATE_CONTCAR)
#Initialize Coordination object. Repeat is necessary so it doesn't count itself
CO_overlap = PDOS_OVERLAP(GAS_PDOS, REFERENCE_PDOS, adsorbate_indices\
                          , site_indices, min_occupation=0.9\
                          , upshift=0.5, energy_weight=4)
    
#######################################################################################
# Plot projected density
# ----------------------
#
# We plot the projected density of the gas, adsorbate, and adsorption site.
CO_overlap.plot_projected_density()

#######################################################################################
# Find the optimal upshift factor
# -------------------------------
#
# The optimal upshift factor shifts the molecular orbital energies to
# minimize the sum the orbital scores used in matching gas and adsorbate orbitals.
# This has the effect of increasing certainty and roughly corresponds to the 
# average shift in molecular orbital energies when a gas adsorbs to the surface
# as a fraction of the fermi energy.
optimized_upshift = CO_overlap.optimize_energy_shift(bound=[-0.5,1.5]\
                                                     , reset=True, plot=True)
print(optimized_upshift)
 
#######################################################################################
# Print orbital CO_overlap attributes
# -----------------------------------
#
# Differences in features are used in computing orbital scores. 
# Scores are used to map gas molecular orbitals ot adsorbate molecular orbitals.
print('Print molecular gas and adsorbate orbital features, respectively.')
print(CO_overlap.gas_features)
print(CO_overlap.adsorbate_features)
print('#####################################################################')
print('Orbital matching scores')
print(CO_overlap.orbital_scores)
print('#####################################################################')
print('Gas to adsorbate indices and band centers')
print(CO_overlap.gas_2_adsorbate)

#######################################################################################
# Plot energy overlap
# -------------------
# We select energy overlap histograms with the adsorbate molecular orbitals
# that influence spectra. Gas orbitals 1,2, and 3 interact with the surface.
gas_indices = [i for i in range(5) if CO_overlap.gas_2_adsorbate[i][0] in [1,2,3]]
adsorbate_indices = CO_overlap.gas_2_adsorbate[gas_indices,1].astype('int')
CO_overlap.plot_energy_overlap(adsorbate_indices)

#######################################################################################
# Print orbital interactions
# --------------------------
# Plot orbital interaction of the first gas molecular orbital with a surface
# s, pz, and dz2 orbitals. These are identified from first figure above
nano = 'Pt44'
nano_DOSCAR = os.path.join(example_path, nano + '/DOSCAR')
nano_CONTCAR = os.path.join(example_path, nano + '/CONTCAR')
#obtain atom indices and atom type as 'surface' or 'bulk'
nano_indices, GCNs, atom_types = get_geometric_data(nano_CONTCAR)
#initialize a PDOS object for the nanoparticle
nano_PDOS = VASP_DOS(nano_DOSCAR)
#calculate orbital interactions
BULK_DOSCAR = os.path.join(example_path,'Pt_nano/Pt147/DOSCAR')
# VASP_DOS objects for both the gas (vacuum) and the adsorbate+surface system
GAS_PDOS = VASP_DOS(GAS_DOSCAR)
REFERENCE_PDOS = VASP_DOS(ADSORBATE_DOSCAR)
BULK_PDOS = VASP_DOS(BULK_DOSCAR)
print('Interactions with 4sigma orbital')
orbital_interaction = CO_overlap.calculate_orbital_interaction(gas_indices[0]\
                    , nano_PDOS, nano_indices[atom_types[...] == 'surface'][0]\
                         , ['s','pz','dz2'], BULK_PDOS, bulk_atom=43\
                             , sum_density=False, sum_spin=True)
print(orbital_interaction)
print('Interactions with 1pi orbital')
orbital_interaction = CO_overlap.calculate_orbital_interaction(gas_indices[1]\
                    , nano_PDOS, nano_indices[atom_types[...] == 'surface'][0]\
                         , ['dyz','dxz'], BULK_PDOS, bulk_atom=43\
                             , sum_density=False, sum_spin=True)
print(orbital_interaction)
print('Interactions with 5sigma orbital')
orbital_interaction = CO_overlap.calculate_orbital_interaction(gas_indices[2]\
                    , nano_PDOS, nano_indices[atom_types[...] == 'surface'][0]\
                         , ['s','pz','dz2'], BULK_PDOS, bulk_atom=43\
                             , sum_density=False, sum_spin=True)
print(orbital_interaction)
