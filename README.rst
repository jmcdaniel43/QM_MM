=========
QM/MM/PME
=========

:author: John Pederson
:author_email: jpederson6@gatech.edu
:project: QM_MM
:date_written: March 28, 2022
:last_date_modified: April 11, 2022

Summary
-------

This package implements single-point energy calculations for the
QM/MM/PME method described by John Pederson and Professor Jesse 
McDaniel:

DOI:10.1063/5.00xxxxx

Installation
------------

This software depends on a `modified fork of openmm
<https://github.com/jmcdaniel43/OpenMM-7.4>`_ and a `modified fork of 
psi4 <https://github.com/jmcdaniel43/psi4>`_.  These modified
repositories must be compiled from source.

The modified openmm requires the following dependencies:

- cython=0.29.17
- doxygen=1.8.18
- swig=3.0.12

The modified psi4 requires the following dependencies:

- gcc>=4.9
- gau2grid=1.3.1
- pint=0.9 or 0.11
- pydantic=1.5.1
- libxc=4.3.4
- numpy>=1.19.2

The QM/MM repository requires the following additional dependencies:

- lxml>=4.5.0
- scipy>=1.4.1

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

Authors
-------

Jesse McDaniel
John Pederson
