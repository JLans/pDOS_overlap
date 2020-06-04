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

#######################################################################################
# Load DOSCAR file
# ----------------
#
# First we will, get the example data, load a DOSCAR file and use it to
# instantiate a VASP_DOS object.

set_figure_settings('paper')
np.set_printoptions(linewidth=100)
#These files are too large to store in the examples directory
example_path = r'C:\Users\lansf\Documents\Data\PROBE_PDOS\vasp_dos_files'
GAS_DOSCAR = os.path.join(example_path, 'CO/DOSCAR')
GAS_CONTCAR = os.path.join(example_path, 'CO/CONTCAR')
ADSORBATE_DOSCAR = os.path.join(example_path, 'CO+Pt111/DOSCAR')
ADSORBATE_CONTCAR = os.path.join(example_path, 'CO+Pt111/CONTCAR')

GAS_PDOS = VASP_DOS(GAS_DOSCAR)
REFERENCE_PDOS = VASP_DOS(ADSORBATE_DOSCAR)

#######################################################################################
# Obtain orbital overlap scores and mapping from gas to adsorbate indices
# -----------------------------------------------------------------------
#
# This method utilizes two VASP_DOS objects, a gas and an adsorption system.
# It uses the adosorbtion system (REFERENCE_PDOS) to map gas molecular orbitals
# to adsorbate molecular orbitals. It then calculates the adsorption site
# atomic orbital energy overlaps with the adsorbate molecular orbital energies.

adsorbate_indices, site_indices = get_adsorbate_indices(GAS_CONTCAR, ADSORBATE_CONTCAR)
#Initialize Coordination object. Repeat is necessary so it doesn't count itself
CO_overlap = PDOS_OVERLAP(GAS_PDOS, REFERENCE_PDOS, adsorbate_indices\
                          , site_indices, min_occupation=0.9)
#print overlap scores
print('Orbital matching scores')
print(CO_overlap.orbital_scores)
print('#####################################################################')
print('Gas to adsorbate indices and band centers')
print(CO_overlap.gas_2_adsorbate)

#######################################################################################
# Plot projected density
# ----------------------
#
# We plot the projected density of the gas, adsorbate, and adsorption site.
CO_overlap.plot_projected_density()

#######################################################################################
# Plot energy overlap
# -------------------
# We select energy overlap histograms with the adsorbate molecular orbitals
# that influence spectra.
indices = [i for i in range(5) if CO_overlap.gas_2_adsorbate[i][0] in [1,2,3]]
adsorbate_indices = CO_overlap.gas_2_adsorbate[indices,1].astype('int')
CO_overlap.plot_energy_overlap(adsorbate_indices)
