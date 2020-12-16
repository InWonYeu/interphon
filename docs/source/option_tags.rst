.. _label_option_tags:

===========
Option tags
===========

As listed below, ***InterPhon*** supports a range of option tags to manage phonon computations and plotting styles.
In order to see all of the available options and their default values in command line::

    $ interphon --help

Because it is not convenient to enter many options in command line,
it is supported to pass a bunch of options through a file ``Option_file`` (file name is arbitrary)::

    $ interphon --option_file Option_file

or::

    $ interphon -option Option_file

The general format of ``Option_file`` is::

    [option tag 1] = [option value 1]
    [option tag 2] = [option value 2]
    [option tag 3] = [option value 3]
    ...

.. _label_option_file_example:

and :ref:`example <label_option_file_example>` of ``Option_file`` is::

    DFT = vasp
    DISP = 0.02
    ENLARGE = 4 4 1
    PBC = 1 1 0

    OPTION_DOS = stack  # A plot option of DOS
    NDOS = 1000  # The number of DOS points
    ELIMIT = -1 8  # Energy (THz) limitation of DOS and Band plot
    ATOM_DOS = 6, 5 4  # The Index of atoms to be projected in DOS plot. Among the selected atoms as interface, the atom in 6th, 5th, and 4th lines is located in the topmost layer (1st layer), 2nd layer, and 3rd layer, respectively
    COLOR_DOS = black  # The color of total DOS line
    LEGEND_DOS = 1st layer, 2nd + 3rd layer  # Legends for the projected atoms in DOS

    OPTION_band = projection  # A plot option of band
    ATOM_BAND = 6  # The Index of atoms to be projected in band plot
    K_LABEL_BAND = G M K G  # The label of high-symmetry k-points
    BAR_LABEL_BAND = 1st layer  # The label of colorbar for the projected atoms

.. note::
   Option tags can be either long or short name, irrespective of lowercase and uppercase letters.

.. _label_basic_option_tags:

Basic option tags
*****************

1. ––option_file, –option
-------------------------
::

    help = File of option collections
    value type = File path

    usage:
    $ interphon -option Option_file

2. ––dft_code, –dft
-------------------
::

    help = DFT code name
    value type = str (one of vasp, espresso, aims)
    default = vasp

    usage:
    $ interphon -dft vasp

3. ––displacement, –disp
------------------------
::

    help = Displacement length (unit: Angst)
    value type = float
    default = 0.02

    usage:
    $ interphon -disp 0.02

4. ––enlargement, –enlarge
--------------------------
::

    help = Extension ratio along each a, b, c lattice direction
    value type = str
    default = '2 2 1'

    usage:
    $ interphon -enlarge "2 2 1"

5. ––periodicity, –pbc
----------------------
::

    help = Periodic (True or 1) or not (False or 0) along each a, b, c lattice direction
    value type = str
    default = '1 1 0'

    usage:
    $ interphon -pbc "1 1 0"

6. ––unitcell, –c
-----------------
::

    help = Unit cell file
    value type = File path
    default = POSCAR

    usage:
    $ interphon -c POSCAR

7. ––supercell, –sc
-------------------
::

    help = Supercell file
    value type = File path

    usage:
    $ interphon -sc SUPERCELL

Density of state (DOS) option tags
**********************************

1. ––density_of_state, –dos
---------------------------
::

    help = Flag to DOS
    value type = bool
    default = False (automatically changed to True if the option -kdos is given)

    usage:
    $ interphon -dos

2. ––kpoint_dos, –kdos
----------------------
::

    help = K-point file for DOS
    value type = File path

    usage:
    $ interphon -kdos KPOINTS_dos

3. ––sigma, –sig
----------------
::

    help = Sigma of gaussian smearing (0.0: tetrahedron method)
    value type = float
    default = 0.1

    usage:
    $ interphon -sig 0.1

4. ––number_dos, –ndos
----------------------
::

    help = The number of DOS points
    value type = int
    default = 200

    usage:
    $ interphon -ndos 200

5. ––projection_atom_dos, –atom_dos
-----------------------------------
::

    help = The Index of atoms to be projected in DOS plot
    value type = str

    usage:
    $ interphon -atom_dos "1 2 3, 4 5 6"

6. ––projection_legend_dos, –legend_dos
---------------------------------------
::

    help = Legends for the projected atoms
    value type = str

    usage:
    $ interphon -legend_dos "1st layer, 2nd layer"

7. ––energy_limit, –elimit
--------------------------
::

    help = Energy (THz) limitation of DOS and Band plot
    value type = str

    usage:
    $ interphon -elimit "-1 8"

8. ––tdos_color_dos, –color_dos
-------------------------------
::

    help = The color of total DOS line
    value type = str
    default = tab:orange (should be supported in matplotlib)

    usage:
    $ interphon -color_dos tab:orange

.. _label_dos_option_dos:

9. ––projection_option_dos, –option_dos
---------------------------------------
::

    help = Option for DOS projection plot
    value type = str (one of plain, line, stack)
    default = plain

    usage:
    $ interphon -option_dos plain

10. ––image_orientation_dos, –orientation_dos
---------------------------------------------
::

    help = Orientation of DOS plot
    value type = str (one of horizontal, vertical)
    default = horizontal

    usage:
    $ interphon -orientation_dos horizontal

11. ––legend_location_dos, –legend_loc_dos
------------------------------------------
::

    help = Location of DOS legend
    value type = str (one of best, upper right, upper left, lower left, lower right, right, center left, center right, lower center, upper center, center)
    default = best

    usage:
    $ interphon -legend_loc_dos "upper right"

Thermal property option tags
****************************

1. ––thermal_property, –thermal
-------------------------------
::

    help = Flag to thermal property
    value type = bool
    default = False

    usage:
    $ interphon -thermal

2. ––temperature_minimum, –tmin
-------------------------------
::

    help = Temperature minimum (unit: K)
    value type = int
    default = 0

    usage:
    $ interphon -tmin 0

3. ––temperature_maximum, –tmax
-------------------------------
::

    help = Temperature maximum (unit: K)
    value type = int
    default = 1000

    usage:
    $ interphon -tmax 1000

4. ––temperature_step, –tstep
-----------------------------
::

    help = Temperature step (unit: K)
    value type = int
    default = 10

    usage:
    $ interphon -tstep 10

Band structure option tags
**************************

1. ––phonon_band, –band
-----------------------
::

    help = Flag to phonon band
    value type = bool
    default = False (automatically changed to True if the option -kband is given)

    usage:
    $ interphon -band

2. ––kpoint_band, –kband
------------------------
::

    help = K-point file for Band
    value type = File path

    usage:
    $ interphon -kband KPOINTS_band

3. ––kpoint_label_band, –k_label_band
-------------------------------------
::

    help = The label of high-symmetry k-points
    value type = str

    usage:
    $ interphon -k_label_band "G M K G"

4. ––projection_atom_band, –atom_band
-------------------------------------
::

    help = The Index of atoms to be projected in Band plot
    value type = str

    usage:
    $ interphon -atom_band "1 2 3"

5. ––total_color_band, –color_band
----------------------------------
::

    help = The color of Band line
    value type = str
    default = tab:orange

    usage:
    $ interphon -color_band tab:orange

.. _label_band_option_band:

6. ––projection_option_band, –option_band
-----------------------------------------
::

    help = Option for Band projection plot
    value type = str (one of plain, projection)
    default = plain

    usage:
    $ interphon -option_band plain

7. ––colorbar_label_band, –bar_label_band
-----------------------------------------
::

    help = The label of colorbar for projection plot
    value type = str

    usage:
    $ interphon -bar_label_band "1st layer"

8. ––colorbar_location_band, –bar_loc_band
------------------------------------------
::

    help = Location of colorbar
    value type = str (one of right, bottom)
    default = right

    usage:
    $ interphon -bar_loc_band right

phonon mode option tags
***********************

1. ––phonon_mode, –mode
-----------------------
::

    help = Flag to phonon mode
    value type = bool
    default = False

    usage:
    $ interphon -mode

2. ––index_mode, –ind_mode
--------------------------
::

    help = The index of phonon mode
    value type = int
    default = 0 (0: the lowest band line, 1: second lowest band line, etc.)

    usage:
    $ interphon -ind_mode 0

3. ––k_point_mode, –kpt_mode
----------------------------
::

    help = The K-point of phonon mode
    value type = str
    default = '0.0 0.0 0.0' (corresponding to Gamma point)

    usage:
    $ interphon -kpt_mode "0.0 0.0 0.0"

.. caution::
   The k-point given by the option –kpt_mode should be included in k-points of band line path.
