===========
Preparation
===========

.. _label_preparation_before_phonon_calculations:

Preparation before phonon calculations
**************************************
Based on the derivation of :ref:`label_phonon_formalism` below, the optimization of structure must be done
before the beginning of phonon calculations (:ref:`label_pre_process`).
The structure optimization (*i.e.*, optimizing the positions of constituent atoms into equilibrium) is carried out by relaxation calculation.
It is generally recommended to use a strict force criterion for the relaxation convergence (*e.g.*, corresponding to EDIFFG_ tag in VASP).
This is because highly accurate values of forces acting on each atom are required for phonon evaluation.

.. _EDIFFG: https://www.vasp.at/wiki/index.php/EDIFFG

.. _label_phonon_formalism:

Phonon formalism
****************

In crystals, the equilibrium position of an atom is described by the labels of lattice (:math:`l`) and basis (:math:`s`):

.. math::

   R_{ls} = R_{l}+\tau_{s}

where :math:`R_{l}` indicates a specific unit cell among infinite lattices and :math:`\tau_{s}` is the position of the atom in the unit cell.
In real situation where vibrational motions around their equilibrium positions {:math:`R_{ls}`} with displacements {:math:`{u_{ls}}`} take place,
the potential energy :math:`U` of a crystal can be written by a Taylor expansion around the equilibrium positions up to the second order (harmonic approximation):

.. math::

   U = U_{0}+\frac{1}{2} \sum_{ls\alpha}\sum_{l's'\beta} \frac{\partial^{2}U}{\partial u_{ls\alpha}\partial u_{l's'\beta}} u_{ls\alpha}u_{l's'\beta}

where :math:`U_0` is the zeroth term corresponding to the potential energy at equilibrium; and :math:`\alpha` and :math:`\beta` indicate the Cartesian directions.
Note that the first derivative of the potential energy (corresponding to force) vanishes, since the force on each atom is zero near the equilibrium position.
The second derivative, called the second order force constants (:math:`K_{ls\alpha,l's'\beta}`),
is numerically evaluated using the midpoint method within the FDM scheme:

.. math::

   K_{ls\alpha,l's'\beta} = \frac{\partial^{2}U}{\partial u_{ls\alpha}\partial u_{l's'\beta}} = -\frac{\partial F_{ls\alpha}}{\partial u_{l's'\beta}}
   = -\frac{F_{ls\alpha}(R_{l's'\beta}+u) - (R_{l's'\beta}-u)}{2u}

where :math:`F_{ls\alpha}(R_{l's'\beta}+u)` is the force on the atom at :math:`R_{ls}` along :math:`\alpha` direction
in response to the displacement of the atom at :math:`R_{l's'}` by :math:`u` along :math:`\beta` direction.

.. _label_dynamical_matrix:

The periodic boundary condition (PBC) of crystals imposes translational invariance on :math:`K_{ls\alpha,l's'\beta}`
(*i.e.*, :math:`\sum_{l'} K_{ls\alpha,l's'\beta} = \sum_{l'} K_{0s\alpha,l's'\beta}`) and plae wave nature on :math:`u_{ls\alpha}` (:math:`=u^0_{s\alpha} e^{i(qR_{ls}-wt)}`)
as for a sound wave. Finally, by making several transformations,
the eigenvalue equation of the :ref:`dynamical matrix <label_dynamical_matrix>` with elements of :math:`D_{s\alpha,s'\beta}(q)` is obtained:

.. math::
   w^2 v^0_{s\alpha} = \sum_{s'\beta} D_{s\alpha,s'\beta}(q) v^0_{s'\beta},
.. math::
   D_{s\alpha,s'\beta}(q) = \frac{1}{(M_s M_s')^{1/2}} \sum_{l'} K_{0s\alpha,l's'\beta} e^{iqR_{l'}} e^{iq(\tau_{s'}-\tau_{s})}

where :math:`w` is the square root of the eigenvalue corresponding to the normal frequency of lattice vibrations,
and a vector with components :math:`v^0_{s\alpha}` (:math:`=M^{1/2}_s u^0_{s\alpha}`) is the eigenvector
and is referred to as the mass-weighted normal mode.
