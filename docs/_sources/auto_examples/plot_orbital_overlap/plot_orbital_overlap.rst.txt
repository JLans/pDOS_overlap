.. only:: html

    .. note::
        :class: sphx-glr-download-link-note

        Click :ref:`here <sphx_glr_download_auto_examples_plot_orbital_overlap_plot_orbital_overlap.py>`     to download the full example code
    .. rst-class:: sphx-glr-example-title

    .. _sphx_glr_auto_examples_plot_orbital_overlap_plot_orbital_overlap.py:


==============================================
Calculating orbital overlap using pdos_overlap
==============================================

This example shows how calculate the overlap of gas phase molecular orbitals
with an adsorbate and surface atom.


.. code-block:: default


    import os
    import numpy as np
    from pdos_overlap import VASP_DOS
    from pdos_overlap.plotting_tools import set_figure_settings
    from pdos_overlap import get_adsorbate_indices
    from pdos_overlap import PDOS_OVERLAP
    from pdos_overlap.coordination import get_geometric_data








Load DOSCAR file
----------------

First we will, get the example data, load a DOSCAR file and use it to
instantiate a VASP_DOS object.


.. code-block:: default

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








Generate VASP_DOS objects
-------------------------

VASP_DOS objects for both the gas (vacuum) and the adsorbate+surface system


.. code-block:: default

    GAS_PDOS = VASP_DOS(GAS_DOSCAR)
    REFERENCE_PDOS = VASP_DOS(ADSORBATE_DOSCAR)








Get adsorbate and site indices and initialize PDOS_OVERLAP object
-----------------------------------------------------------------

This method utilizes two VASP_DOS objects, a gas and an adsorption system.
It uses the adosorbtion system (REFERENCE_PDOS) to map gas molecular orbitals
to adsorbate molecular orbitals. It then calculates the adsorption site
atomic orbital energy overlaps with the adsorbate molecular orbital energies.


.. code-block:: default

    adsorbate_indices, site_indices = get_adsorbate_indices(GAS_CONTCAR\
                                                            , ADSORBATE_CONTCAR)
    #Initialize Coordination object. Repeat is necessary so it doesn't count itself
    CO_overlap = PDOS_OVERLAP(GAS_PDOS, REFERENCE_PDOS, adsorbate_indices\
                              , site_indices, min_occupation=0.9\
                              , upshift=0.5, energy_weight=4)
    







Plot projected density
----------------------

We plot the projected density of the gas, adsorbate, and adsorption site.


.. code-block:: default

    CO_overlap.plot_projected_density()




.. rst-class:: sphx-glr-horizontal


    *

      .. image:: /auto_examples/plot_orbital_overlap/images/sphx_glr_plot_orbital_overlap_001.png
            :class: sphx-glr-multi-img

    *

      .. image:: /auto_examples/plot_orbital_overlap/images/sphx_glr_plot_orbital_overlap_002.png
            :class: sphx-glr-multi-img

    *

      .. image:: /auto_examples/plot_orbital_overlap/images/sphx_glr_plot_orbital_overlap_003.png
            :class: sphx-glr-multi-img


.. rst-class:: sphx-glr-script-out

 Out:

 .. code-block:: none

    C:\Users\lansf\Box Sync\Synced_Files\Coding\Python\Github\pdos_overlap\pdos_overlap\pdos_overlap.py:899: UserWarning: Matplotlib is currently using agg, which is a non-GUI backend, so cannot show the figure.
      plt.show()
    C:\Users\lansf\Box Sync\Synced_Files\Coding\Python\Github\pdos_overlap\pdos_overlap\pdos_overlap.py:899: UserWarning: Matplotlib is currently using agg, which is a non-GUI backend, so cannot show the figure.
      plt.show()
    C:\Users\lansf\Box Sync\Synced_Files\Coding\Python\Github\pdos_overlap\pdos_overlap\pdos_overlap.py:899: UserWarning: Matplotlib is currently using agg, which is a non-GUI backend, so cannot show the figure.
      plt.show()




Find the optimal upshift factor
-------------------------------

The optimal upshift factor shifts the molecular orbital energies to
minimize the sum the orbital scores used in matching gas and adsorbate orbitals.
This has the effect of increasing certainty and roughly corresponds to the 
average shift in molecular orbital energies when a gas adsorbs to the surface
as a fraction of the fermi energy.


.. code-block:: default

    optimized_upshift = CO_overlap.optimize_energy_shift(bound=[-0.5,1.5]\
                                                         , reset=True, plot=True)
    print(optimized_upshift)
 



.. image:: /auto_examples/plot_orbital_overlap/images/sphx_glr_plot_orbital_overlap_004.png
    :class: sphx-glr-single-img


.. rst-class:: sphx-glr-script-out

 Out:

 .. code-block:: none

    C:\Users\lansf\Box Sync\Synced_Files\Coding\Python\Github\pdos_overlap\pdos_overlap\pdos_overlap.py:813: UserWarning: Matplotlib is currently using agg, which is a non-GUI backend, so cannot show the figure.
      plt.show()
    0.36169076485686813




Print orbital CO_overlap attributes
-----------------------------------

Differences in features are used in computing orbital scores. 
Scores are used to map gas molecular orbitals ot adsorbate molecular orbitals.


.. code-block:: default

    print('Print molecular gas and adsorbate orbital features, respectively.')
    print(CO_overlap.gas_features)
    print(CO_overlap.adsorbate_features)
    print('#####################################################################')
    print('Orbital matching scores')
    print(CO_overlap.orbital_scores)
    print('#####################################################################')
    print('Gas to adsorbate indices and band centers')
    print(CO_overlap.gas_2_adsorbate)





.. rst-class:: sphx-glr-script-out

 Out:

 .. code-block:: none

    Print molecular gas and adsorbate orbital features, respectively.
    [[-2.59440624e+01  1.55667344e+00  3.26199490e-13  4.22017815e-01  9.65840299e-15]
     [-1.07725355e+01  6.77718574e-01  9.85077691e-12  8.64998584e-01  2.91694728e-13]
     [-8.37158857e+00  3.93060231e-12  1.38208964e+00  7.22459794e-12  1.38201564e+00]
     [-5.77108906e+00  3.57464960e-01  1.06157881e-12  8.50897818e-01  3.14380969e-14]
     [ 1.18455181e+00  3.22212492e-02  9.60645175e-01  2.46394073e-02  9.60608089e-01]]
    [[-2.41110858e+01  1.73902895e+00  5.09101559e-10  4.56123751e-01  1.26260274e-09]
     [-1.05406293e+01  8.09275578e-01  4.59590029e-07  2.72026686e-01  2.63096436e-07]
     [-7.89123682e+00  1.81585882e-01  3.65520754e-04  1.18150089e+00  2.68352149e-04]
     [-7.00986643e+00  7.00095672e-03  1.16194079e+00  7.04962282e-02  1.15795451e+00]
     [ 7.22337248e-01  1.47708586e-01  1.34778394e+00  2.14257133e-01  1.35510153e+00]]
    #####################################################################
    Orbital matching scores
    [[9.29028568e-01 1.61273572e-04 7.34545868e-08 1.38664310e-07 6.45021613e-09]
     [1.61273572e-04 1.99724666e-01 1.75777184e-02 4.97812655e-03 1.90226298e-05]
     [7.34545868e-08 1.75777184e-02 3.02183417e-04 6.88050593e-01 1.82339250e-03]
     [1.38664310e-07 4.97812655e-03 6.88050593e-01 2.21233455e-03 6.35174227e-03]
     [6.45021613e-09 1.90226298e-05 1.82339250e-03 6.35174227e-03 5.64799688e-01]]
    #####################################################################
    Gas to adsorbate indices and band centers
    [[  0.           0.         -29.20625951 -24.28842369]
     [  1.           1.         -14.03473256 -10.71795641]
     [  2.           3.         -11.63378565  -7.18747574]
     [  3.           2.          -9.03328614  -8.06828387]
     [  4.           4.          -2.07764527   0.54343512]]




Plot energy overlap
-------------------
We select energy overlap histograms with the adsorbate molecular orbitals
that influence spectra. Gas orbitals 1,2, and 3 interact with the surface.


.. code-block:: default

    gas_indices = [i for i in range(5) if CO_overlap.gas_2_adsorbate[i][0] in [1,2,3]]
    adsorbate_indices = CO_overlap.gas_2_adsorbate[gas_indices,1].astype('int')
    CO_overlap.plot_energy_overlap(adsorbate_indices)




.. rst-class:: sphx-glr-horizontal


    *

      .. image:: /auto_examples/plot_orbital_overlap/images/sphx_glr_plot_orbital_overlap_005.png
            :class: sphx-glr-multi-img

    *

      .. image:: /auto_examples/plot_orbital_overlap/images/sphx_glr_plot_orbital_overlap_006.png
            :class: sphx-glr-multi-img

    *

      .. image:: /auto_examples/plot_orbital_overlap/images/sphx_glr_plot_orbital_overlap_007.png
            :class: sphx-glr-multi-img


.. rst-class:: sphx-glr-script-out

 Out:

 .. code-block:: none

    C:\Users\lansf\Box Sync\Synced_Files\Coding\Python\Github\pdos_overlap\pdos_overlap\pdos_overlap.py:844: UserWarning: Matplotlib is currently using agg, which is a non-GUI backend, so cannot show the figure.
      plt.show()




Print orbital interactions
--------------------------
Plot orbital interaction of the first gas molecular orbital with a surface
s, pz, and dz2 orbitals. These are identified from first figure above


.. code-block:: default

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




.. rst-class:: sphx-glr-script-out

 Out:

 .. code-block:: none

    Interactions with 4sigma orbital
    [-0.71071628 -0.31695864 -0.15137432]
    Interactions with 1pi orbital
    [-0.30455614 -0.17911093]
    Interactions with 5sigma orbital
    [-0.41681493 -0.11738777 -0.10148584]





.. rst-class:: sphx-glr-timing

   **Total running time of the script:** ( 0 minutes  27.921 seconds)


.. _sphx_glr_download_auto_examples_plot_orbital_overlap_plot_orbital_overlap.py:


.. only :: html

 .. container:: sphx-glr-footer
    :class: sphx-glr-footer-example



  .. container:: sphx-glr-download sphx-glr-download-python

     :download:`Download Python source code: plot_orbital_overlap.py <plot_orbital_overlap.py>`



  .. container:: sphx-glr-download sphx-glr-download-jupyter

     :download:`Download Jupyter notebook: plot_orbital_overlap.ipynb <plot_orbital_overlap.ipynb>`


.. only:: html

 .. rst-class:: sphx-glr-signature

    `Gallery generated by Sphinx-Gallery <https://sphinx-gallery.github.io>`_
