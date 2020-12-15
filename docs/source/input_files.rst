===========
Input files
===========

.. _label_dft_structure_file:

DFT structure file
******************
Within ``InterPhon``, the interfacial region is supposed to be defined through the statement of constraints on atom movements (selective dynamics).
Phonon evaluation proceeds only in the selected atoms.

See below example of Cu(111) surface where the top three layers are selected as the surface region.
A part of the DFT structure file is shown in VASP, Quantum ESPRESSO, and FHI-aims format.

**VASP**::

    Unknown
    1.00000000000000
        2.5712952614000000    0.0000000000000000    0.0000000000000000
        1.2856476307000000    2.2268070170000001    0.0000000000000000
        0.0000000000000000    0.0000000000000000   27.7901687622000004
    Cu
    7
    Selective dynamics
    Cartesian
        0.0000000000000000    0.0000000000000000    6.2983610849999998   F   F   F  # atom index = 0
        2.5712951080000002    1.4845379230000000    8.3978147799999991   F   F   F  # atom index = 1
        1.2856475540000001    0.7422689609999999   10.5098654109999998   F   F   F  # atom index = 2
        0.0000000000000000    0.0000000000000000   12.6153545979999997   F   F   F  # atom index = 3
        2.5712951080000002    1.4845379230000000   14.7267644020000006   T   T   T  # atom index = 4
        1.2856475540000001    0.7422689609999999   16.8249826169999999   T   T   T  # atom index = 5
        0.0000000000000000    0.0000000000000000   18.9121107990000006   T   T   T  # atom index = 6

**Quantum ESPRESSO**::

    CELL_PARAMETERS angstrom
        2.5712952614000000    0.0000000000000000    0.0000000000000000
        1.2856476307000000    2.2268070170000001    0.0000000000000000
        0.0000000000000000    0.0000000000000000   27.7901687622000004

    ATOMIC_POSITIONS angstrom
    Cu   0.0000000000000000    0.0000000000000000    6.2983610849999998  0 0 0
    Cu   2.5712951080000002    1.4845379230000000    8.3978147799999991  0 0 0
    Cu   1.2856475540000001    0.7422689609999999   10.5098654109999998  0 0 0
    Cu   0.0000000000000000    0.0000000000000000   12.6153545979999997  0 0 0
    Cu   2.5712951080000002    1.4845379230000000   14.7267644020000006
    Cu   1.2856475540000001    0.7422689609999999   16.8249826169999999
    Cu   0.0000000000000000    0.0000000000000000   18.9121107990000006

**FHI-aims**::

    lattice_vector   2.5712952614000000    0.0000000000000000    0.0000000000000000
    lattice_vector   1.2856476307000000    2.2268070170000001    0.0000000000000000
    lattice_vector   0.0000000000000000    0.0000000000000000   27.7901687622000004
    atom   0.0000000000000000    0.0000000000000000    6.2983610849999998 Cu
        constrain_relaxation .true.
    atom   2.5712951080000002    1.4845379230000000    8.3978147799999991 Cu
        constrain_relaxation .true.
    atom   1.2856475540000001    0.7422689609999999   10.5098654109999998 Cu
        constrain_relaxation .true.
    atom   0.0000000000000000    0.0000000000000000   12.6153545979999997 Cu
        constrain_relaxation .true.
    atom   2.5712951080000002    1.4845379230000000   14.7267644020000006 Cu
    atom   1.2856475540000001    0.7422689609999999   16.8249826169999999 Cu
    atom   0.0000000000000000    0.0000000000000000   18.9121107990000006 Cu

.. _label_dft_force_file:

DFT force file
**************
If the DFT calculations (electronic relaxation) are successfully converged,
output files containing atomic forces are created in all DFT programs (*e.g.*, ``vasprun.xml`` and ``OUTCAR`` in VASP).
The force information is used to fill the :ref:`dynamical matrix <label_dynamical_matrix>` elements.

.. _label_kpoint_file:

K-point file
************
The :ref:`label_kpoint_file`, which is supported in VASP KPOINTS_ format, is used for the mesh sampling of k-points.

.. _KPOINTS: https://www.vasp.at/wiki/index.php/KPOINTS

**KPOINTS_dos** (file name is arbitrary)::

    kpoint
    0
    MP  # Monkhorst-Pack grids, use the first character ‘G’ for Gamma-centered grids.
    9 9 1
    0.0 0.0 0.0

**KPOINTS_band** (file name is arbitrary)::

    kpoint
    41
    L  # Line path
    0.00 0.00 0.00  # G
    0.00 0.50 0.00  # M

    0.00 0.50 0.00  # M
    0.333333 0.333333 0.00  # K

    0.333333 0.333333 0.00  # K
    0.00 0.00 0.00  # G
