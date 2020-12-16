============
Introduction
============

**InterPhon** allows efficient extraction of interfacial phonons at arbitrary reciprocal points in 1D or 2D
Brillouin zone, by automatically processing the information obtained by 3D-based DFT codes.

The only requirements for package execution are a :ref:`label_dft_structure_file` representing the structure of interest
and the corresponding :ref:`label_dft_force_file` providing the forces acting on each atom.

Overview
********

The computational cost of *ab initio* phonon calculations is proportional to the number of atoms per unit cell.
In conventional *ab initio* interface calculations, however, a large number of atoms and a low degree of symmetry are
intrinsically inevitable to represent interfacial region in a 3D periodic cell. As a result,
the ability to apply *ab initio* phonon calculations to interfaces has been limited
by excessive computational cost imposed by their large number of atoms and broken symmetry.

Among all of these atoms, however,
the bonding environment of most of the atoms located far from the interface is similar to
that of the bulk atoms, presenting the same vibrations as the ideal bulk.
Deviation of vibrations from bulk occurs only on the atoms in the vicinity of the interface.
Therefore, to alleviate the excessive computational cost,
**InterPhon** focuses on the interfacial atoms by allowing users to easily select atoms to be considered as the interface.

Strategy
********

Phonon evaluation of **InterPhon** is based on direct approach using the finite displacement method (FDM).

In contrast to conventional phonon codes for 3D bulk systems such as Phonopy_ and PHON_, where generate the dynamical matrix of
all constituent atoms in a unit cell by requiring the ‘complete FDM’, **InterPhon**, which requires ‘selective FDM’,
generates the dynamical matrix of a subset of the constituent atoms in a unit cell and does not assume 3D periodicity.

.. _Phonopy: https://phonopy.github.io/phonopy/
.. _PHON: https://www.sciencedirect.com/science/article/pii/S0010465509001064

Architecture
************

The **InterPhon** package comprises following five sub-packages:

- ``error`` package consists of error modules defined by developer to guide users.

- ``inout`` package consists of parser modules to read and write files in different DFT formats.

- ``core`` package consists of core modules responsible for process control.

- ``util`` package consists of util modules to extend the functionality.

- ``analysis`` package consists of analysis modules to characterize phonons.
