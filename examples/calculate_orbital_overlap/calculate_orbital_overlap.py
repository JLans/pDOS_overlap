# -*- coding: utf-8 -*-
"""
Created on Fri May 29 13:45:17 2020

@author: lansf
"""

"""
==============================================
Calculating orbital overlap using pdos_overlap
==============================================

This example shows how calculate the overlap of gas phase molecular orbitals
with an adsorbate and surface atom.

"""

import os
from vasp_dos import get_example_data
from vasp_dos import VASP_DOS
from vasp_dos.plotting_tools import set_figure_settings
from vasp_dos import get_adsorbate_indices
from vasp_dos import PDOS_OVERLAP

#######################################################################################
# Load DOSCAR file
# ----------------
#
# First we will, get the example data, load a DOSCAR file and use it to
# instantiate a VASP_DOS object.

set_figure_settings('paper')
example_path = get_example_data()
GAS_DOSCAR = os.path.join(example_path, 'CO/DOSCAR')
GAS_CONTCAR = os.path.join(example_path, 'CO/CONTCAR')
ADSORBATE_DOSCAR = os.path.join(example_path, 'CO+Pt111/DOSCAR')
ADSORBATE_CONTCAR = os.path.join(example_path, 'CO+Pt111/CONTCAR')

GAS_PDOS = VASP_DOS(GAS_DOSCAR)
REFERENCE_PDOS = VASP_DOS(ADSORBATE_DOSCAR)

#######################################################################################
# Calculate and print band centers
# --------------------------------
#
# This method uses the the site and spin orbital projected density. It sums the
# spin orbital densities to get energy sub-level band centers.

adsorbate_indices = get_adsorbate_indices(GAS_CONTCAR, ADSORBATE_CONTCAR)
CO_overlap = PDOS_OVERLAP(GAS_PDOS, REFERENCE_PDOS,adsorbate_indices=adsorbate_indices)

