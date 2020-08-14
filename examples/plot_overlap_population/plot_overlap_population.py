"""
=================================================================
Plotting cyrstal orbital overlap population obtained from lobster
=================================================================

This example shows how to plot overlap population data
See http://www.cohp.de/ for details
"""

import os
import numpy as np
import matplotlib.pyplot as plt
from pdos_overlap.overlap_population import get_example_data
from pdos_overlap.overlap_population import OVERLAP_POPULATION
from pdos_overlap.plotting_tools import set_figure_settings

#######################################################################################
# Load COOPCAR file
# -----------------
#
# First we will, get the example data, load a COOPCAR file and use it to
# instantiate an OVERLAP_POPULATION object

set_figure_settings('paper')
example_path = get_example_data()
COOPCAR = os.path.join(example_path, 'C2H4/COOPCAR.lobster')

POP = OVERLAP_POPULATION(COOPCAR)

#######################################################################################
# Identify bonding interactions and check for spin
# ------------------------------------------------
#
print(POP.interactions)
print(POP.is_spin)


#######################################################################################
# Obtain projected overlap
# ------------------------
#
# We projected orbital overlap for the C-C bond and C-H bonds in C2H4
# We group the CH bonds and ensure to sum for spins as all electrons are paired

CC_overlap = POP.get_pcoop(interactions=[0], sum_pcoop=False, sum_spin=True)
CH_overlap = POP.get_pcoop(interactions=[1,2,3,4], sum_pcoop=True, sum_spin=True)

#######################################################################################
# Plot the bonding populaiton with respect to the CC and CH bonds
# ---------------------------------------------------------------
#
# A positive value on the x-axis indicates are greater proportion of states in
# in the bond than outside of the bond

plt.figure(figsize=(3,5))
plt.plot(CC_overlap, POP.get_energies(), zorder=3)
plt.plot(CH_overlap, POP.get_energies(), zorder=2)
plt.plot([np.min([CC_overlap, CH_overlap]), np.max([CC_overlap, CH_overlap])]\
         ,[POP.e_fermi, POP.e_fermi],'k--', zorder=1, linewidth=5)
plt.legend(['C-C overlap population','C-H overlap population','fermi level'],loc='best')
plt.xlabel('Orbital overlap')
plt.ylabel('Energy [eV]')
plt.show()
