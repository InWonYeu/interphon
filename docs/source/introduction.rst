============
Introduction
============

**InterPhon** allows the efficient extraction of interfacial phonons at arbitrary wave vectors in 1D or 2D
reciprocal space, the corresponding real space being a 1D or in-plane 2D periodic system,
by automatically processing the information obtained by 3D-based DFT codes.

The only requirements for the program execution are a :ref:`label_dft_structure_file` representing the atomic structure
of interest and the corresponding :ref:`label_dft_force_file` providing the forces acting on each atom.

Overview
********

The computational cost of *ab initio* phonon calculations is proportional to the number of atoms per unit cell.
In conventional *ab initio* interface calculations, however, a large number of atoms and a low degree of symmetry are
intrinsically inevitable to represent interfacial region in a 3D periodic cell. As a result,
the ability to apply *ab initio* phonon calculations to interfaces has been limited by excessive computational cost.

Among all of these atoms, however,
the bonding environment of most of the atoms located far from the interface is similar to
that of the bulk atoms, presenting the same vibrations as the ideal bulk.
Deviation of vibrations from bulk occurs only on the atoms in the vicinity of the interface.
Therefore, using atomic forces obtained from a 3D periodic DFT framework,
**InterPhon** evaluates the distinct interfacial phonons by allowing users
to define an interfacial region that is smaller than the full simulation box,
such that the phonons are calculated only for the atoms in the interface region.

Strategy
********

**InterPhon** operation is based on direct approach using the finite displacement method (FDM).

In contrast to conventional phonon codes for 3D bulk systems such as Phonopy_ and PHON_, which generate the dynamical matrix of
all constituent atoms in a unit cell by requiring the ‘complete FDM’, **InterPhon**, which requires ‘selective FDM’,
generates the dynamical matrix of a subset of the constituent atoms in a unit cell and does not assume 3D periodicity.

.. _Phonopy: https://phonopy.github.io/phonopy/
.. _PHON: https://www.sciencedirect.com/science/article/pii/S0010465509001064

Architecture
************

The **InterPhon** package consists of following five sub-packages:

- ``error`` sub-package includes error modules defined by developer to guide users.

- ``inout`` sub-package includes inout modules to read and write files in different DFT formats.

- ``core`` sub-package includes core modules responsible for the central processes.

- ``util`` sub-package includes util modules to extend the functionality.

- ``analysis`` sub-package includes analysis modules to characterize phonon properties.
