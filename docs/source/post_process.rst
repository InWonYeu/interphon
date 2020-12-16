.. _label_post_process:

=================================
Phonon calculations: Post-process
=================================

Overview
********

After the DFT force calculations (*e.g.*, IBRION_ = -1 tag in VASP corresponding to electronic relaxation at fixed atomic positions)
for the displaced supercells (``POSCAR-0*``) are finished in each *FORCE-0** folder,
the evaluation and analysis of interfacial phonons can be executed.

.. _IBRION: https://www.vasp.at/wiki/index.php/IBRION

In this post-process execution, the :ref:`dynamical matrix <label_dynamical_matrix>` at each :ref:`K-point <label_kpoint_file>`
is internally constructed by reading the :ref:`label_dft_force_file`.
The internal processes end with printing of the phonon properties,
such as density of states (DOS), thermal properties (vibrational entropy and free energy),
band structures, and phonon modes in the complete forms of both data files and graphics.

Post-process execution
**********************
As an example, let's assume the :ref:`label_dft_force_file` is ``vasprun.xml`` in VASP format
and each force file resides in *Phonon root/FORCE-0** folders.
First, go to the *Phonon root* folder where the ``POSCAR`` file resides.
Then, execute **InterPhon** post-process by one of the following ways: :ref:`label_post_process_command_line` and :ref:`label_post_process_python_interpreter`.
If done successfully, files of :ref:`phonon properties <label_post_process_property_file>` (``band.dat``, ``band.png``, etc.)
and :ref:`Post-process record file <label_post_process_record_file>` (``post_process.yaml``) will be generated.

.. note::
   :ref:`Pre-process record file <label_pre_process_record_file>` (``pre_process.yaml``) and supercell file (``SUPERCELL``),
   which are generated in previous :ref:`Pre-process <label_pre_process>`,
   should reside in the *Phonon root* folder with the :ref:`label_dft_structure_file` of a targeted unit cell.

.. _label_post_process_command_line:

1. Command line
---------------
The following commands is for the simultaneous analysis for all of the supported phonon properties.
If you don't need some analyses, exclude corresponding flag options (:ref:`label_option_tags`).

.. note::

   The :ref:`label_kpoint_file` (``KPOINTS_dos`` and ``KPOINTS_band`` below),
   which is supported in VASP KPOINTS_ format, is required for the mesh sampling of k-points.

.. _KPOINTS: https://www.vasp.at/wiki/index.php/KPOINTS

General long format::

    $ interphon [DFT_force_files] --dft_code ['vasp' or 'espresso' or 'aims'] --kpoint_dos KPOINTS_dos --density_of_state --thermal_property --kpoint_band KPOINTS_band --phonon_band --phonon_mode

General short format::

    $ interphon [DFT_force_files] -dft ['vasp' or 'espresso' or 'aims'] -kdos KPOINTS_dos -dos -thermal -kband KPOINTS_band -band -mode

Simple usage with the help of default setting of :ref:`label_option_tags`::

    $ interphon FORCE-0*/vasprun.xml -kdos KPOINTS_dos -thermal -kband KPOINTS_band -mode

.. _label_post_process_python_interpreter:

2. Python interpreter
---------------------
>>> from InterPhon.core import PostProcess
>>> import glob
>>> user_args = {'dft_code': 'vasp', 'displacement': 0.02, 'enlargement': "4 4 1", 'periodicity': "1 1 0"}
>>> force_path = glob.glob('FORCE-0*/vasprun.xml')
>>> force_path.sort()
>>> post_dos = PostProcess('POSCAR', 'SUPERCELL', code_name=user_args.get('dft_code'))
>>> post_dos.set_user_arg(user_args)
>>> post_dos.set_reciprocal_lattice()
>>> post_dos.set_force_constant(force_path, code_name=user_args.get('dft_code'))
>>> post_dos.set_k_points('KPOINTS_dos')
>>> post_dos.eval_phonon()
>>> from InterPhon.analysis import DOS
>>> post_dos.dos = DOS(post_dos)
>>> post_dos.dos.set()
>>> post_dos.dos.write()
>>> post_dos.dos.plot()
>>> from InterPhon.analysis import ThermalProperty
>>> post_dos.thermal = ThermalProperty(post_dos)
>>> post_dos.thermal.set()
>>> post_dos.thermal.write()
>>> post_dos.thermal.plot()
>>> from copy import deepcopy
>>> post_band = deepcopy(post_dos)
>>> post_band.set_k_points('KPOINTS_band')
>>> post_band.eval_phonon()
>>> from InterPhon.analysis import Band
>>> post_band.band = Band(post_band)
>>> post_band.band.set()
>>> post_band.band.write()
>>> post_band.band.plot()
>>> post_band.band.plot_with_dos(dos_object=post_dos.dos)
>>> from InterPhon.analysis import Mode
>>> post_band.mode = Mode(post_band)
>>> post_band.mode.set()
>>> post_band.mode.write()
>>> post_band.mode.plot(unit_cell='POSCAR', code_name=user_args.get('dft_code'))  # This requires (Optional) ASE
