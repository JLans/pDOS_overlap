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
gas = 'C2H4'
adsorbate = 'C2H4_bridge'
surface = 'Pt111'
set_figure_settings('paper')
np.set_printoptions(linewidth=100)
#These files are too large to store in the examples directory
lobster_path = r'C:\Users\lansf\Documents\Data\PROBE_PDOS\lobster_files'
GAS_DOSCAR = os.path.join(lobster_path, gas + '/DOSCAR.lobster')
GAS_CONTCAR = os.path.join(lobster_path, gas + '/CONTCAR')
ADSORBATE_DOSCAR = os.path.join(lobster_path, 'gas+Pt_G.03_noW/'+surface + '+'\
                          + adsorbate + '/DOSCAR.lobster')
ADSORBATE_CONTCAR = os.path.join(lobster_path, 'gas+Ptnano/'+surface + '+'\
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
C2H4_overlap = PDOS_OVERLAP(GAS_PDOS, REFERENCE_PDOS, reference_indices\
                          , site_indices, min_occupation=1\
                          , upshift=0.5, energy_weight=3)
    
#######################################################################################
# Plot projected density
# ----------------------
#
# We plot the projected density of the gas, adsorbate, and adsorption site.
C2H4_overlap.plot_projected_density()

#######################################################################################
# Find the optimal upshift factor
# -------------------------------
#
# The optimal upshift factor shifts the gas molecular orbital energies to
# minimize the sum the orbital scores used in matching gas and adsorbate orbitals.
# This has the effect of increasing certainty and roughly corresponds to the 
# average shift in molecular orbital energies when a gas adsorbs to the surface
"""
I need to fix this for NO
"""
optimized_upshift = C2H4_overlap.optimize_energy_shift(bound=[-10, 10]\
                                                     , reset=True, plot=True)
print(optimized_upshift)
 
#######################################################################################
# Print orbital CO_overlap attributes
# -----------------------------------
#
# Differences in features are used in computing orbital scores. 
# Scores are used to map gas molecular orbitals ot adsorbate molecular orbitals.
print('Print molecular gas and adsorbate orbital features, respectively.')
print(C2H4_overlap.gas_features)
print(C2H4_overlap.adsorbate_features)
print('#####################################################################')
print('Orbital matching scores')
print(C2H4_overlap.orbital_scores)
print('#####################################################################')
print('Gas to adsorbate indices and band centers')
print(C2H4_overlap.gas_2_adsorbate)

#######################################################################################
# Identify bonding orbitals
# -------------------------
#
# We calcluate the amount of density for each orbital that is in a bonding region
# We can do this both for the gas and for the adsorbate

#gas
COOPCAR_NO = os.path.join(lobster_path, gas + '/COOPCAR.lobster')
POP_NO = OVERLAP_POPULATION(COOPCAR_NO)
bonding_fraction = POP_NO.get_bonding_fraction(C2H4_overlap.gas_orbital_indices\
                                               , C2H4_overlap.GAS_PDOS.get_energies()\
                                               , set_antibonding_zero=False)
print('Gas bonding fraction')
print(bonding_fraction)
    
#adsorbate
COOPCAR_NO = os.path.join(lobster_path, 'gas+Pt_G.03_noW/'+surface + '+'\
                          + adsorbate + '/COOPCAR.lobster')
POP_NO = OVERLAP_POPULATION(COOPCAR_NO)
bonding_fraction = POP_NO.get_bonding_fraction(C2H4_overlap.adsorbate_orbital_indices\
                                               , C2H4_overlap.REFERENCE_PDOS.get_energies()\
                                               , set_antibonding_zero=True
                                               , emax = C2H4_overlap.REFERENCE_PDOS.e_fermi)
print('Adsorbate bonding fraction')
print(bonding_fraction)
print(C2H4_overlap.adsorbate_occupations)
print(C2H4_overlap.adsorbate_band_centers)

#######################################################################################
# Plot energy overlap
# -------------------
#
# We select energy overlap histograms with the adsorbate molecular orbitals
# that influence spectra. Gas orbitals 1,2, and 3 interact with the surface.
C2H4_overlap.plot_energy_overlap([0,1,2,3,4], ['s','d'])
