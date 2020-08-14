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
import matplotlib.pyplot as plt
Downloads_folder = os.path.join(os.path.expanduser("~"),'Downloads')
#######################################################################################
# Load DOSCAR file
# ----------------
#
# First we will, get the example data, load a DOSCAR file and use it to
# instantiate a VASP_DOS object.
gas = 'gases/C2H4'
adsorbate = 'C2H4_bridge'
surface = 'Pt111'
set_figure_settings('paper')
np.set_printoptions(linewidth=100)
#These files are too large to store in the examples directory
lobster_path = r'C:\Users\lansf\Documents\Data\PROBE_PDOS\lobster_files_(N+1)bands'
GAS_DOSCAR = os.path.join(lobster_path, gas + '/DOSCAR.lobster')
GAS_CONTCAR = os.path.join(lobster_path, gas + '/CONTCAR')
ADSORBATE_DOSCAR = os.path.join(lobster_path, 'surfaces_noW/'+surface + '+'\
                          + adsorbate + '/DOSCAR.lobster')
ADSORBATE_CONTCAR = os.path.join(lobster_path, 'surfaces_noW/'+surface + '+'\
                          + adsorbate + '/CONTCAR')

#######################################################################################
# Generate VASP_DOS objects
# -------------------------
#
# VASP_DOS objects for both the gas (vacuum) and the adsorbate+surface system
GAS_PDOS = VASP_DOS(GAS_DOSCAR)
REFERENCE_PDOS = VASP_DOS(ADSORBATE_DOSCAR)
#REFERENCE_PDOS.apply_gaussian_filter(10)

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
CO_overlap = PDOS_OVERLAP(GAS_PDOS, REFERENCE_PDOS, reference_indices\
                          , site_indices, min_occupation=0.55\
                          , upshift=0.5, energy_weight=4)
    
#######################################################################################
# Plot projected density
# ----------------------
#
# We plot the projected density of the gas, adsorbate, and adsorption site.
CO_overlap.plot_projected_density(figure_directory=Downloads_folder)

#######################################################################################
# Find the optimal upshift factor
# -------------------------------
#
# The optimal upshift factor shifts the gas molecular orbital energies to
# minimize the sum the orbital scores used in matching gas and adsorbate orbitals.
# This has the effect of increasing certainty and roughly corresponds to the 
# average shift in molecular orbital energies when a gas adsorbs to the surface
optimized_upshift = CO_overlap.optimize_energy_shift(bound=[-10,10]\
                                                     , reset=True, plot=True)
print(optimized_upshift)
 

#######################################################################################
# Identify bonding orbitals
# -------------------------
#
# We calcluate the amount of density for each orbital that is in a bonding region
# We can do this both for the gas and for the adsorbate

#gas
COOPCAR_CO = os.path.join(lobster_path, gas + '/COOPCAR.lobster')
POP_CO_GAS = OVERLAP_POPULATION(COOPCAR_CO)
bonding_states = POP_CO_GAS.get_bonding_states(CO_overlap.gas_orbital_indices\
                                               , CO_overlap.GAS_PDOS.get_energies()\
                                               , set_antibonding_zero=True)
print('Gas bonding states')
print(bonding_states)
    
#adsorbate
COOPCAR_CO = os.path.join(lobster_path, 'surfaces_noW/'+surface + '+'\
                          + adsorbate + '/COOPCAR.lobster')
POP_CO_ADSORBATE = OVERLAP_POPULATION(COOPCAR_CO)
bonding_states = POP_CO_ADSORBATE.get_bonding_states(CO_overlap.adsorbate_orbital_indices\
                                               , CO_overlap.REFERENCE_PDOS.get_energies()\
                                               , set_antibonding_zero=True
                                               , emax = CO_overlap.REFERENCE_PDOS.e_fermi)
print('Adsorbate bonding states')
print(bonding_states)

bonding_states = POP_CO_ADSORBATE.get_bonding_states(CO_overlap.adsorbate_orbital_indices
                                               , CO_overlap.REFERENCE_PDOS.get_energies()
                                               , interactions = [12]
                                               , set_antibonding_zero=False
                                               , emax = CO_overlap.REFERENCE_PDOS.e_fermi)
print('C-O bonding states')
print(bonding_states)
print(CO_overlap.adsorbate_band_centers)
print(CO_overlap.adsorbate_occupations)
#######################################################################################
# Plot energy overlap
# -------------------
#
# We select energy overlap histograms with the adsorbate molecular orbitals
# that influence spectra. Gas orbitals 1,2, and 3 interact with the surface.
CO_overlap.plot_energy_overlap(indices=[0,1,2,3], atomic_orbitals=['s', 'd']
                               , figure_directory=Downloads_folder)

#######################################################################################
# Obtain projected overlap
# ------------------------
#
# We projected orbital overlap for the C-C bond and C-H bonds in C2H4
# We group the CH bonds and ensure to sum for spins as all electrons are paired
GAS_OVERLAP = POP_CO_GAS.get_pcoop([0], sum_pcoop=True)
ADSORBATE_OVERLAP = POP_CO_ADSORBATE.get_pcoop(sum_pcoop=True,set_antibonding_zero=True)
CO_OVERLAP = POP_CO_ADSORBATE.get_pcoop([12],sum_pcoop=True)

#######################################################################################
# Plot the bonding populaiton with respect to the CC and CH bonds
# ---------------------------------------------------------------
#
# A positive value on the x-axis indicates are greater proportion of states in
# in the bond than outside of the bond

fig = plt.figure(figsize=(7.2,4),dpi=400)
abc = ['(a)','(b)','(c)']
axes = fig.subplots(nrows=1, ncols=3)
axes_list = [axes[0], axes[1], axes[2]]
#plotting function
def plot_density(OVERLAP, energies, e_fermi, index):
    axes_list[index].plot(OVERLAP, energies, zorder=2)
    axes_list[index].plot([np.min(OVERLAP), np.max(OVERLAP)]
             ,[e_fermi, e_fermi], 'k--', zorder=1, linewidth=5) 
    axes_list[index].text(0.90,0.96,abc[index],transform=axes_list[index].transAxes)
    
#plot gas density
plot_density(GAS_OVERLAP, POP_CO_GAS.get_energies(), POP_CO_GAS.e_fermi, 0)
#plot adsorbate density
plot_density(CO_OVERLAP, POP_CO_ADSORBATE.get_energies(), POP_CO_ADSORBATE.e_fermi, 1)
#plot adsorption-site density
plot_density(ADSORBATE_OVERLAP, POP_CO_ADSORBATE.get_energies(), POP_CO_ADSORBATE.e_fermi, 2)
fig.text(0.001, 0.5, 'Energy [eV]', va='center', rotation='vertical')
fig.text(0.5, 0.01, 'Overlap density [states/eV]', ha='center')
figure_path = os.path.join(Downloads_folder,'pccop.jpg')
fig.set_tight_layout({'pad':2,'w_pad':1})
plt.savefig(figure_path, format='jpg')
plt.close()
