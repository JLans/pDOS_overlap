"""
==================================
Calculating band center using vdos
==================================

This example shows how to plot projected density of states

"""

import os
from pdos_overlap import get_example_data
from pdos_overlap import VASP_DOS
from pdos_overlap.plotting_tools import set_figure_settings

#######################################################################################
# Load DOSCAR file
# ----------------
#
# First we will, get the example data, load a DOSCAR file and use it to
# instantiate a VASP_DOS object.

set_figure_settings('paper')
example_path = get_example_data()
DOSCAR = os.path.join(example_path, 'C2H4/DOSCAR')
PDOS = VASP_DOS(DOSCAR)

#######################################################################################
# Calculate and print band centers
# --------------------------------
#
# This method uses the the site and spin orbital projected density. It sums the
# spin orbital densities to get energy sub-level band centers.

orbitals = [key for key in PDOS.orbital_dictionary.keys() if 's' in key or 'p' in key]
    
band_centers = PDOS.get_band_center([0], orbital_list=orbitals\
                                    , max_energy=PDOS.e_fermi)

for count, orbital in enumerate(orbitals):
    print(orbital + ' band center :' + str(band_centers[count]))
