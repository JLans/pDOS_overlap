"""
==============================================
Calculating orbital overlap using pdos_overlap
==============================================

This example shows how calculate the overlap of gas phase molecular orbitals \
with an adsorbate and surface atom.

"""

import os
import numpy as np
from pdos_overlap.vasp_dos import VASP_DOS
from pdos_overlap.plotting_tools import set_figure_settings
from pdos_overlap import get_adsorbate_indices
from pdos_overlap import PDOS_OVERLAP
from pdos_overlap.overlap_population import OVERLAP_POPULATION

#######################################################################################
# Load DOSCAR file
# ----------------
#
# First we will, get the example data, load a DOSCAR file and use it to
# instantiate a VASP_DOS object.
gas = 'NO'
adsorbate = 'NO_fcc'
surface = 'Pt111'
set_figure_settings('paper')
np.set_printoptions(linewidth=100)
#These files are too large to store in the examples directory
lobster_path = r'C:\Users\lansf\Documents\Data\PROBE_PDOS\lobster_files'
GAS_DOSCAR = os.path.join(lobster_path, gas + '/DOSCAR.lobster')
GAS_CONTCAR = os.path.join(lobster_path, gas + '/CONTCAR')
ADSORBATE_DOSCAR = os.path.join(lobster_path, 'gas+Pt_G.03_noW/'+surface + '+'\
                          + adsorbate + '/DOSCAR.lobster')
ADSORBATE_CONTCAR = os.path.join(lobster_path, 'gas+Pt_G.03_noW/'+surface + '+'\
                          + adsorbate + '/CONTCAR')

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
reference_indices, site_indices = get_adsorbate_indices(GAS_CONTCAR\
                                                        , ADSORBATE_CONTCAR)
#Initialize Coordination object. Repeat is necessary so it doesn't count itself
NO_overlap = PDOS_OVERLAP(GAS_PDOS, REFERENCE_PDOS, reference_indices\
                          , site_indices, min_occupation=1\
                          , upshift=0.5, energy_weight=3)
    
#######################################################################################
# Plot projected density
# ----------------------
#
# We plot the projected density of the gas, adsorbate, and adsorption site.
NO_overlap.plot_projected_density()

#######################################################################################
# Find the optimal upshift factor
# -------------------------------
#
# The optimal upshift factor shifts the gas molecular orbital energies to
# minimize the sum the orbital scores used in matching gas and adsorbate orbitals.
# This has the effect of increasing certainty and roughly corresponds to the 
# average shift in molecular orbital energies when a gas adsorbs to the surface
optimized_upshift = NO_overlap.optimize_energy_shift(bound=[-10, 10]\
                                                     , reset=True, plot=True)
print(optimized_upshift)
 
#######################################################################################
# Print orbital CO_overlap attributes
# -----------------------------------
#
# Differences in features are used in computing orbital scores. 
# Scores are used to map gas molecular orbitals ot adsorbate molecular orbitals.
print('Print molecular gas and adsorbate orbital features, respectively.')
print(NO_overlap.gas_features)
print(NO_overlap.adsorbate_features)
print('#####################################################################')
print('Orbital matching scores')
print(NO_overlap.orbital_scores)
print('#####################################################################')
print('Gas to adsorbate indices and band centers')
print(NO_overlap.gas_2_adsorbate)

#######################################################################################
# Identify bonding orbitals
# -------------------------
#
# We calcluate the amount of density for each orbital that is in a bonding region
# We can do this both for the gas and for the adsorbate

#gas
COOPCAR_NO = os.path.join(lobster_path, gas + '/COOPCAR.lobster')
POP_NO = OVERLAP_POPULATION(COOPCAR_NO)
bonding_states = POP_NO.get_bonding_states(NO_overlap.gas_orbital_indices\
                                               , NO_overlap.GAS_PDOS.get_energies()\
                                               , set_antibonding_zero=False)
print('Gas bonding states')
print(bonding_states)
    
#adsorbate
COOPCAR_NO = os.path.join(lobster_path, 'gas+Pt_G.03_noW/'+surface + '+'\
                          + adsorbate + '/COOPCAR.lobster')
POP_NO = OVERLAP_POPULATION(COOPCAR_NO)
bonding_states = POP_NO.get_bonding_states(NO_overlap.adsorbate_orbital_indices\
                                               , NO_overlap.REFERENCE_PDOS.get_energies()\
                                               , set_antibonding_zero=True
                                               , emax = NO_overlap.REFERENCE_PDOS.e_fermi)
print('Adsorbate bonding fraction')
print(bonding_states)

bonding_states = POP_NO.get_bonding_states(NO_overlap.adsorbate_orbital_indices
                                               , NO_overlap.REFERENCE_PDOS.get_energies()
                                               , interactions=[6]
                                               , set_antibonding_zero=True
                                               , emax = NO_overlap.REFERENCE_PDOS.e_fermi)
print('N-O bonding fraction')
print(bonding_states)
#######################################################################################
# Plot energy overlap
# -------------------
#
# We select energy overlap histograms with the adsorbate molecular orbitals
# that influence spectra. Gas orbitals 1,2, and 3 interact with the surface.
NO_overlap.plot_energy_overlap([1,2,3,4], ['s','d'])
