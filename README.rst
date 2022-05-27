=========
QM/MM/PME
=========

:Author: John Pederson
:Author Email: jpederson6@gatech.edu
:Project: QM_MM
:Date Written: March 28, 2022
:Last Date Modified: May 26, 2022

Summary
-------

This package implements single-point energy calculations for the
QM/MM/PME method described by John Pederson and Professor Jesse 
McDaniel:

DOI: `10.1063/5.0087386 <https://aip.scitation.org/doi/10.1063/5.0087386>`_

Installation
------------

This software depends on a `modified fork of openmm
<https://github.com/jmcdaniel43/OpenMM-7.4>`_ and a `modified fork of 
psi4 <https://github.com/jmcdaniel43/psi4>`_.  These modified
repositories must be compiled from source.

The modified openmm requires the following dependencies:

- cython
- doxygen
- swig

The modified psi4 requires the following dependencies:

- gcc>=4.9
- gau2grid
- pint
- pydantic
- libxc
- numpy>=1.19.2

The QM/MM repository requires the following additional dependencies:

- lxml
- scipy

Once the modified psi4 and openmm repositories are built, the QM_MM
repository may be cloned.  No source requires compilation.

Usage
-----

Any python run file should include this directory in its path in order
to use this software.  An anaconda environment should point towards the
modified openmm and psi4 repositories.  Required input files include PDB,
topology XML, and forcefield XML for the OpenMM interface.  All other
options may be passed as options in the instantiation of the
MMSubsystem, QMSubsystem, and QMMMSystem objects.

Documentation
-------------

The documentation for this project can be found `here
<http://johnppederson.com/qm_mm/html/index.html>`_!

Authors
-------

Jesse McDaniel

John Pederson
