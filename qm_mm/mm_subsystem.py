#! /usr/bin/env python3
# -*- coding: utf-8 -*-
"""
OpenMM interface to model the MM subsystem of the QM/MM system.
"""
import sys

import numpy as np
import simtk.openmm
import simtk.unit

from .shared import *
from .system import System
from .utils import *


class MMSubsystem(MM_base, System):
    """
    Child class of the MM base and System classes for QM/MM simulations.

    Must use compiled version of custom OpenMM to allow for QM/MM/PME.

    Parameters
    ----------
    pdb_list: list of str
        The directories containing the PDB files which define the system
        geometry.
    residue_xml_list: list of str
        The directories containing the XML files which define the system
        topology.
    ff_xml_list: list of str
        The directories containing the XML files which define the system
        interactions in the MM Hamiltonian.
    platform: str
        The platform for OpenMM.  These include "Reference", "CPU",
        "CUDA", and "OpenCL".
    qmmm_pme: bool, Optional, default=False
        Determine whether or not to include extended electrostatics
        using the QM/MM/PME method.
    qmmm_pme_gridnumber: int, Optional, default=100
        The number of grid points to use along each of the box vectors
        of the principal cell during the PME procedure.
    qmmm_pme_alpha: float, Optional, default=5.0
        The Gaussian width of the smeared point charges in the Ewald
        summation scheme, in inverse nanometers.  See OpenMM's
        documentation for further discussion.
    kwargs: dict
        Additional arguments to send to the MM base class.
    """

    def __init__(
            self, 
            pdb_list,
            residue_xml_list,
            ff_xml_list,
            platform,
            qmmm_pme = False,
            qmmm_pme_gridnumber = 100,
            qmmm_pme_alpha = 5.0,
            **kwargs,
        ):
        MM_base.__init__(
            self,
            pdb_list,
            residue_xml_list,
            ff_xml_list,
            **kwargs,
        )
        System.__init__(self)
        self.args = [pdb_list, residue_xml_list, ff_xml_list]
        self.kwargs = kwargs
        self.qmmm_pme = qmmm_pme
        self.qmmm_pme_gridnumber = qmmm_pme_gridnumber
        self.qmmm_pme_alpha = qmmm_pme_alpha
        # Sets the PME parameters in OpenMM.  The grid size is important
        # for the accuracy of the external potental in the DFT 
        # quadrature, since this is interpolated from the PME grid.
        if self.qmmm_pme:
            self.nbondedForce.setPMEParameters(
                self.qmmm_pme_alpha,
                self.qmmm_pme_gridnumber,
                self.qmmm_pme_gridnumber,
                self.qmmm_pme_gridnumber,
            )
        properties = {}
        if self.qmmm_pme:
            properties["ReferenceVextGrid"] = "true"
        if platform == "Reference":
            self.platform = openmm.Platform.getPlatformByName("Reference")
        elif platform == "CPU":
            self.platform = openmm.Platform.getPlatformByName("CPU")
        elif platform == "OpenCL":
            self.platform = openmm.Platform.getPlatformByName("OpenCL")
            if self.qmmm_pme:
                print(
                    """Only Reference and CPU OpenMM platforms are
                    currently supported for the QM/MM/PME
                    implementation."""
                )
                sys.exit()
        elif platform == "CUDA":
            self.platform = openmm.Platform.getPlatformByName("CUDA")
            properties["Precision"] = "mixed"
            if self.qmmm_pme:
                print(
                    """Only Reference and CPU OpenMM platforms are
                    currently supported for the QM/MM/PME
                    implementation."""
                )
                sys.exit()
        else:
            print("Platform '{}' is unrecognized.".format(platform))
            sys.exit()
        self.simulation = simtk.openmm.app.Simulation(
            self.modeller.topology,
            self.system,
            self.integrator,
            self.platform,
            properties,
        )
        self.simulation.context.setPositions(self.modeller.positions)
        residue_part_list = []
        element_symbols = []
        charges = []
        for residue in self.pdb.topology.residues():
            part_list = []
            for part in residue._atoms:
                part_list.append(part.index)
                element_symbols.append(part.element.symbol)
                (q, sig, eps) = self.nbondedForce.getParticleParameters(
                    part.index
                )
                charges.append(q._value)
            residue_part_list.append(part_list)
        self.residue_part_list = residue_part_list
        self.element_symbols = element_symbols
        self.charges = charges
        if self.qmmm_pme:
            state = self.simulation.context.getState(
                getEnergy=True,
                getForces=True,
                getPositions=True,
                getVext_grids=True,
                getPME_grid_positions=True
            )
        else:
            state = self.simulation.context.getState(
                getEnergy=False,
                getForces=False,
                getVelocities=False,
                getPositions=True,
            )
        self.box = state.getPeriodicBoxVectors(True) / simtk.unit.angstrom
        self.positions = state.getPositions(True) / simtk.unit.angstrom

    def build_potential_grid(self):
        """
        Return the PME external potential grid.

        Returns
        -------
        potential_grid: List of float
            The PME external potential grid which OpenMM uses to
            calculate electrostatic interactions, in 
            kJ/mol/proton charge.
        """
        state = self.simulation.context.getState(
            getEnergy=True,
            getForces=True,
            getVelocities=True,
            getPositions=True,
            getVext_grids=True,
            getPME_grid_positions=True,
        )
        potential_grid = state.getVext_grid()
        return potential_grid

    def write_pdb(self, name):
        """
        Write the current state of the subsystem to a PDB file.

        Parameters
        ----------
        name: str
            The name of output file.
        """
        state = self.simulation.context.getState(
            getEnergy=False,
            getForces=False,
            getVelocities=False,
            getPositions=True,
            enforcePeriodicBox=True,
        )
        positions = state.getPositions()
        self.simulation.topology.setPeriodicBoxVectors(
            state.getPeriodicBoxVectors(),
        )
        simtk.openmm.app.PDBFile.writeFile(
            self.simulation.topology,
            positions,
            open(name + '.pdb','w'),
        )
