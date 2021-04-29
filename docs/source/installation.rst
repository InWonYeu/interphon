============
Installation
============

Requirements
************

* `Python>=3 <https://www.python.org/>`_
* `Click <https://click.palletsprojects.com/en/7.x/>`_
* `Numpy <https://numpy.org/doc/stable/>`_
* `Matplotlib <https://matplotlib.org/>`_
* `PyYAML <https://pyyaml.org/>`_
* (Optional) ASE_

.. note::
    The :ref:`visualization of phonon mode <label_ase_visualization>` works with ASE_.
    If you want this function, append 'ase' to install_requires in ``setup.py``::

        install_requires.append('ase')

.. _ASE: https://wiki.fysik.dtu.dk/ase/index.html

Installation using pip
**********************

Install stable version of **InterPhon** via pip_, which gets the code stored in PyPI_::

    $ pip install interphon

.. _PyPI: https://pypi.org/project/InterPhon/
.. _PIP: https://pip.pypa.io/en/stable/

Installation from source code
*****************************

.. :Git clone:

Alternatively, install the latest version in development from `InterPhon GitHub <https://github.com/InWonYeu/interphon>`_ by git clone and setup::

    $ git clone https://github.com/inwonyeu/interphon.git
    $ cd interphon/
    $ python setup.py install

Installation test
*****************

Once the **InterPhon** package is successfully installed, you can import **InterPhon** within Python interpreter::

    >>> import InterPhon

Or the following command can be executed (with an error in the unprepared state)::

    $ interphon


