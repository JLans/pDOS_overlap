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
data_path = r'C:\Users\lansf\Documents\Data\PROBE_PDOS\lobster_files/surfaces_noW'
COOPCAR_C2H4_atop = os.path.join(data_path, 'Pt111+C2H4_atop/COOPCAR.lobster')
COOPCAR_C2H4_bridge = os.path.join(data_path, 'Pt111+C2H4_bridge/COOPCAR.lobster')
COOPCAR_CO = os.path.join(data_path, 'Pt111+CO/COOPCAR.lobster')
COOPCAR_NO_atop = os.path.join(data_path, 'Pt111+NO_atop/COOPCAR.lobster')
COOPCAR_NO_fcc = os.path.join(data_path, 'Pt111+NO_fcc/COOPCAR.lobster')
POP_C2H4_atop = OVERLAP_POPULATION(COOPCAR_C2H4_atop)
POP_C2H4_bridge = OVERLAP_POPULATION(COOPCAR_C2H4_bridge)
POP_CO = OVERLAP_POPULATION(COOPCAR_CO)
POP_NO_atop = OVERLAP_POPULATION(COOPCAR_NO_atop)
POP_NO_fcc = OVERLAP_POPULATION(COOPCAR_NO_fcc)


#######################################################################################
# Obtain projected overlap
# ------------------------
#
# We projected orbital overlap for the C-C bond and C-H bonds in C2H4
# We group the CH bonds and ensure to sum for spins as all electrons are paired

POP_C2H4_atop.overlap = POP_C2H4_atop.get_pcoop(sum_pcoop=True, sum_spin=True\
                                                , set_antibonding_zero=True)
POP_C2H4_bridge.overlap = POP_C2H4_bridge.get_pcoop(sum_pcoop=True, sum_spin=True\
                                                    , set_antibonding_zero=True)
POP_CO.overlap = POP_CO.get_pcoop(sum_pcoop=True, sum_spin=True\
                                  , set_antibonding_zero=True)
POP_NO_atop.overlap= POP_NO_atop.get_pcoop(sum_pcoop=True, sum_spin=True\
                                           , set_antibonding_zero=True)
POP_NO_fcc.overlap = POP_NO_fcc.get_pcoop(sum_pcoop=True, sum_spin=True\
                                          , set_antibonding_zero=True)

#######################################################################################
# Plot the bonding populaiton with respect to the CC and CH bonds
# ---------------------------------------------------------------
#
# A positive value on the x-axis indicates are greater proportion of states in
# in the bond than outside of the bond

for OVERLAP in [POP_C2H4_atop, POP_C2H4_bridge, POP_CO, POP_NO_atop, POP_NO_fcc]:
    plt.figure(figsize=(3,5))
    plt.plot(OVERLAP.overlap, OVERLAP.get_energies(), zorder=3)
    plt.plot([np.min(OVERLAP.overlap), np.max(OVERLAP.overlap)]\
             ,[OVERLAP.e_fermi, OVERLAP.e_fermi],'k--', zorder=1, linewidth=5)
    plt.legend(['overlap population','fermi level'],loc='best')
    plt.xlabel('Orbital overlap')
    plt.ylabel('Energy [eV]')
    plt.show()