#! /usr/bin/env python3
# -*- coding: utf-8 -*-
"""
QM/MM system to handle interactions between the QM and MM subsystems.
"""
import sys

import numpy as np

from .mm_subsystem import MMSubsystem
from .pbc_subsystem import PBCSubsystem
from .system import System
from .utils import *


class QMMMSystem(System):
    """
    QM/MM system object to handle interactions between QM and MM regions.

    Parameters
    ----------
    qm_subsystem: QMSubsystem object
        Contains all relevant information for calculations involving the
        QM subsystem.
    mm_subsystem: MMSubsystem object
        Contains all relevant information for calculations involving the
        MM subsystem.
    group_part_dict: dict of list of list of int
        The indices of particles in the system grouped by a given key.
        Example keys include "qm_atom", "qm_drude", and "analytic".
    embedding_cutoff: float or list, Optional, default=0
        The short-range electrostatic embedding cutoff.  Must be a float,
        unless the embedding_method is specified to be "hybrid", in
        which case the cutoff must be a two-element list of floats.  The
        lesser of the two floats is taken to be the analytic embedding
        cutoff.  The greater of the two floats is taken to be a 
        real-space embedding cutoff.  Real-space embedding is performed
        on the particles between the analytic embedding and real-space
        embedding cutoffs, and analytic embedding is performed on the
        particles within the analytic embedding cutoff.
    embedding_method: str, Optional, default="none"
        The method to use for short-range electrostatic embedding.
        Options include "analytic", "realspace", and "hybrid".  Analytic
        embedding incorporates short-range partial point charges as
        one-electron operators within the QM Hamiltonian.  Real-space
        embedding incorporates short-range partial point charges as 
        non-truncated Coulomb interactions between smeared point charges
        within the PME potential grid, which is interpolated into the
        exchange correlation functional quadrature for the QM/MM/PME
        method.  Hybrid embedding employs both analytic and real-space
        embedding at different cutoffs.
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
    """
    # Enforce inheritance of properties from System class so that setter
    # functions may be overridden.
    for attr in dir(System):
        if type(eval("System." + attr)) == property:
            exec(attr + " = System." + attr)

    def __init__(
            self,
            qm_subsystem,
            mm_subsystem,
            group_part_dict,
            embedding_cutoff=0,
            embedding_method="none",
            qmmm_pme=False,
            qmmm_pme_gridnumber=100,
            qmmm_pme_alpha=5.0,
        ):
        System.__init__(self)
        self.qm_subsystem = qm_subsystem
        self.mm_subsystem = mm_subsystem
        self.subsystems = [self.qm_subsystem, self.mm_subsystem]
        # Collect residue, element, and charge data from OpenMM.
        self.residue_part_list = self.mm_subsystem.residue_part_list
        self.group_part_dict = group_part_dict
        self.element_symbols = self.mm_subsystem.element_symbols
        self.charges = self.mm_subsystem.charges
        self.positions = self.mm_subsystem.positions
        self.box = self.mm_subsystem.box
        for residue in self._residue_part_list:
            if all(atom in self._group_part_dict["qm_atom"]
                   for atom in residue):
                self._group_part_dict["qm_drude"] = list(
                    set(residue)-set(self._group_part_dict["qm_atom"]),
                )
        # Default embedding settings.
        self.embedding_cutoff = embedding_cutoff
        self.embedding_method = embedding_method
        self.qmmm_pme = qmmm_pme
        self.qmmm_pme_gridnumber = qmmm_pme_gridnumber
        self.qmmm_pme_alpha = qmmm_pme_alpha
        if self.qmmm_pme:
            supported_platforms = ["Reference", "CPU"]
            # Ensure that there is an MM subsystem which can provide the
            # external potential grid.
            if self.mm_subsystem.platform.getName() not in supported_platforms:
                self.rs_subsystem = MMSubsystem(
                    *self.mm_subsystem.args,
                    "CPU",
                    qmmm_pme=True,
                    qmmm_pme_gridnumber=self.qmmm_pme_gridnumber,
                    qmmm_pme_alpha=self.qmmm_pme_alpha,
                    **self.mm_subsystem.kwargs,
                )
            else:
                self.rs_subsystem = self.mm_subsystem
            self.pbc_subsystem = PBCSubsystem(
                qmmm_pme_gridnumber = self.qmmm_pme_gridnumber,
                qmmm_pme_alpha = self.qmmm_pme_alpha,
                group_part_dict = self._group_part_dict,
                charges = self._charges,
                positions = self._positions,
                box = self._box,
                n_threads = self.qm_subsystem.n_threads,
            )
            self.subsystems.append(self.rs_subsystem)
            self.subsystems.append(self.pbc_subsystem)

    def single_point_calculation(self, tare=False, forces=False):
        """
        Perform a single-point energy/force calculation for the system.

        The energy calculated is effectively the partially
        self-consistent electrostatic interaction energy of the QM
        subsystem and MM subsystem (and extended subsystem, if the
        QM/MM/PME method is being used).

        Parameters
        ----------
        tare: bool, Optional, default=False
            Determine whether or not to tare the gas-phase energy from
            the calculation.
        forces: bool, Optional, default=False
            Determine whether or not to perform Psi4 gradient
            calculation for forces.

        Returns
        -------
        energy: float
            The electronic energy of the QM subsystem.
        forces: NumPy Array object, Optional
            The forces acting on the QM subsystem and electrostatic
            forces acting on the analytically embedded point charges if
            using analytic embedding.
        """
        if self.embedding_method == "hybrid":
            if len(self.embedding_cutoff) != 2:
                print(
                    """Hybrid embedding requires two embedding cutoffs!
                    Please pass a tuple or list of two embedding cutoffs
                    when instantiating the QMMMSystem object with hybrid
                    embedding."""
                )
                sys.exit()
            embedding_list_low = self.generate_embedding_list(
                threshold=min(self.embedding_cutoff),
            )
            self.group_part_dict["analytic"] = embedding_list_low
            embedding_list_high = self.generate_embedding_list(
                threshold=max(self.embedding_cutoff),
            )
            embedding_list_real = list(
                set(embedding_list_high)-set(embedding_list_low)
            )
            self.group_part_dict["realspace"] = embedding_list_real
        elif self.embedding_method != "none":
            embedding_list = self.generate_embedding_list()
            self.group_part_dict[self.embedding_method] = embedding_list
        self.qm_subsystem.generate_geometry()
        arguments = {}
        if self.qmmm_pme:
            potential_grid = self.rs_subsystem.build_potential_grid()
            ref_quadrature = self.qm_subsystem.build_ref_quadrature()
            if forces:
                (quad_extd_pot, nuc_extd_pot, nuc_extd_grad) = self.pbc_subsystem.build_extd_pot(
                    ref_quadrature,
                    potential_grid,
                    return_grad=True,
                )
            else:
                (quad_extd_pot, nuc_extd_pot) = self.pbc_subsystem.build_extd_pot(
                    ref_quadrature,
                    potential_grid,
                )
            arguments["quad_extd_pot"] = quad_extd_pot
            arguments["nuc_extd_pot"] = nuc_extd_pot
        self.qm_subsystem.compute_energy(**arguments)
        energy = self.qm_subsystem.energy
        if tare:
            self.qm_subsystem.compute_gas_phase_energy()
            energy -= self.qm_subsystem.gas_phase_energy
        if forces:
            arguments["nuc_extd_grad"] = nuc_extd_grad
            self.qm_subsystem.compute_forces(**arguments)
            forces = self.qm_subsystem.forces
            return (energy, forces)
        else:
            return energy

    def generate_embedding_list(self, threshold=None):
        """
        Create the embedding list for the current state of the system.

        The distances from the QM atoms are computed using the centroid
        of the non-QM molecule from the centroid of the QM atoms.  The
        legacy method involves computing distances using the first atom
        position from the non-QM molecule instead.

        Parameters
        ----------
        threshold: float, Optional, default=None
            The threshold distance within which to include particles in
            the embedding list.  If no threshold is specified, the
            embedding_cutoff attribute is taken to be the threshold.

        Returns
        -------
        embedding_list: list of tuple of int
            The list of embedded particles, arranged by residue in
            tuples.
        """
        if not threshold:
            threshold = self.embedding_cutoff
        positions = self._positions
        qm_atom_list = self._group_part_dict["qm_atom"]
        qm_centroid = [sum([positions[atom][k] for atom in qm_atom_list])
                       / len(qm_atom_list) for k in range(3)]
        embedding_list = []
        for residue in self._residue_part_list:
            nth_centroid = [sum([positions[atom][k] for atom in residue])
                            / len(residue) for k in range(3)]
            # Legacy embedding
            #nth_centroid = positions[residue[0]]
            r_vector = least_mirror_vector(
                nth_centroid,
                qm_centroid,
                self._box,
            )
            distance = sum([x**2 for x in r_vector])**0.5
            is_qm = any(atom in residue for atom in qm_atom_list)
            if distance < threshold and not is_qm:
                    embedding_list.append(tuple(residue))
        return embedding_list
    
    @residue_part_list.setter
    def residue_part_list(self, residue_part_list):
        self._residue_part_list = residue_part_list
        for subsystem in self.subsystems:
            subsystem.residue_part_list = residue_part_list

    @group_part_dict.setter
    def group_part_dict(self, group_part_dict):
        self._group_part_dict = group_part_dict
        for subsystem in self.subsystems:
            subsystem.group_part_dict = group_part_dict

    @element_symbols.setter
    def element_symbols(self, element_symbols):
        self._element_symbols = element_symbols
        for subsystem in self.subsystems:
            subsystem.element_symbols = element_symbols
    
    @charges.setter
    def charges(self, charges):
        self._charges = charges
        for subsystem in self.subsystems:
            subsystem.charges = charges

    @positions.setter
    def positions(self, positions):
        self._positions = positions
        for subsystem in self.subsystems:
            subsystem.positions = positions
    
    @box.setter
    def box(self, box):
        self._box = box
        self._bohr_box = [[k*self.angstrom_to_bohr for k in vector]
                          for vector in box]
        for subsystem in self.subsystems:
            subsystem.box = box
