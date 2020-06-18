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
from pdos_overlap.overlap_population import OVERLAP_POPULATION
from pdos_overlap.plotting_tools import set_figure_settings

#######################################################################################
# Load COOPCAR file
# -----------------
#
# First we will, get the example data, load a COOPCAR file and use it to
# instantiate an OVERLAP_POPULATION object

set_figure_settings('paper')
data_path = r'C:\Users\lansf\Documents\Data\PROBE_PDOS\lobster_files'
COOPCAR_C2H4 = os.path.join(data_path, 'C2H4/COOPCAR.lobster')
COOPCAR_CO = os.path.join(data_path, 'CO/COOPCAR.lobster')
COOPCAR_NO = os.path.join(data_path, 'NO/COOPCAR.lobster')
POP_C2H4 = OVERLAP_POPULATION(COOPCAR_C2H4)
POP_CO = OVERLAP_POPULATION(COOPCAR_CO)
POP_NO = OVERLAP_POPULATION(COOPCAR_NO)


#######################################################################################
# Obtain projected overlap
# ------------------------
#
# We projected orbital overlap for the C-C bond and C-H bonds in C2H4
# We group the CH bonds and ensure to sum for spins as all electrons are paired

CC_overlap = POP_C2H4.get_pcoop(interactions=[0], sum_pcoop=False, sum_spin=True)
CH_overlap = POP_C2H4.get_pcoop(interactions=[1,2,3,4], sum_pcoop=True, sum_spin=True)
CO_overlap = POP_CO.get_pcoop(sum_spin=True)
NO_overlap = POP_NO.get_pcoop(sum_spin=True)

#######################################################################################
# Plot the bonding populaiton with respect to the CC and CH bonds
# ---------------------------------------------------------------
#
# A positive value on the x-axis indicates are greater proportion of states in
# in the bond than outside of the bond

plt.figure(figsize=(3,5))
plt.plot(CC_overlap, POP_C2H4.get_energies(), zorder=3)
plt.plot(CH_overlap, POP_C2H4.get_energies(), zorder=2)
plt.plot([np.min([CC_overlap, CH_overlap]), np.max([CC_overlap, CH_overlap])]\
         ,[POP_C2H4.e_fermi, POP_C2H4.e_fermi],'k--', zorder=1, linewidth=5)
plt.legend(['C-C overlap population','C-H overlap population','fermi level'],loc='best')
plt.xlabel('Orbital overlap')
plt.ylabel('Energy [eV]')
plt.show()

plt.figure(figsize=(3,5))
plt.plot(CO_overlap, POP_CO.get_energies(), zorder=2)
plt.plot([np.min(CO_overlap), np.max(CO_overlap)]\
         ,[POP_CO.e_fermi, POP_CO.e_fermi],'k--', zorder=1, linewidth=5)
plt.legend(['CO overlap population','fermi level'],loc='best')
plt.xlabel('Orbital overlap')
plt.ylabel('Energy [eV]')
plt.show()

plt.figure(figsize=(3,5))
plt.plot(NO_overlap, POP_NO.get_energies(), zorder=2)
plt.plot([np.min(NO_overlap), np.max(NO_overlap)]\
         ,[POP_NO.e_fermi, POP_NO.e_fermi],'k--', zorder=1, linewidth=5)
plt.legend(['NO overlap population','fermi level'],loc='best')
plt.xlabel('Orbital overlap')
plt.ylabel('Energy [eV]')
plt.show()