#! /usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Psi4 interface to model the QM subsystem of the QM/MM system.
"""
import copy
import sys

import psi4
import psi4.core

from .system import System
from .utils import *


class QMSubsystem(System):
    """
    Psi4 interface for the QM subsystem.

    Parameters
    ----------
    basis_set: str
        Name of desired Psi4 basis set.
    functional: str
        Name of desired Psi4 density functional.
    quadrature_spherical: int
        Number of spherical (angular and azimuthal) points for the 
        exchange-correlation functional quadrature.
    quadrature_radial: int
        Number of radial points for the exchange-correlation functional
        quadrature.
    qm_charge: int
        Charge of the QM system in proton charge.
    qm_spin: int
        Spin of the QM system.
    scf_type: str, Optional, default="df"
        The SCF type for the calculation, as per Psi4 documentation.
    qmmm_pme: bool, Optional, default=False
        Determine whether or not to use the QM/MM/PME method.
    n_threads: int, Optional, default=1
        Number of threads across which to parallelize the QM calculation.
    read_guess: bool, Optional, default=False
        Determine whether to base calculations on previous wavefunction
        objects (default) or on the Psi4 default SAD.
    group_part_dict: dict of list of list of int, Optional, default=None
        The indices of particles in the system grouped by a given key.
        Example keys include "qm_atom", "qm_drude", and "analytic".
    element_symbols: list of str, Optional, default=None
        The element symbols of all particles in the system.
    charges: NumPy Array object, Optional, default=None
        The charges of all particles in the system, in proton charge
        units.
    positions: NumPy Array object, Optional, default=None
        The positions of all particles in the system, in Angstroms.
    box: NumPy Array object, Optional, default=None
        The box vectors defining the periodic system, in Angstroms.
    """

    def __init__(
            self,
            basis_set,
            functional,
            quadrature_spherical,
            quadrature_radial,
            qm_charge,
            qm_spin,
            scf_type="df",
            qmmm_pme=False,
            n_threads=1,
            read_guess=False,
            group_part_dict=None,
            element_symbols=None,
            charges=None,
            positions=None,
            box=None,
        ):
        System.__init__(self)
        self.basis_set = basis_set
        self.functional = functional
        self.quadrature_spherical = quadrature_spherical
        self.quadrature_radial = quadrature_radial
        self.qm_charge = qm_charge
        self.qm_spin = qm_spin
        self.scf_type = scf_type
        self.qmmm_pme = qmmm_pme
        self.n_threads = n_threads
        psi4.set_num_threads(n_threads)
        options = {
            "basis": self.basis_set,
            "dft_spherical_points": self.quadrature_spherical,
            "dft_radial_points": self.quadrature_radial,
            "scf_type": self.scf_type,
            "pme": str(self.qmmm_pme).lower(),
        }
        psi4.set_options(options)
        self.read_guess = read_guess
        self.wfn = None
        self.charge_field = None
        self.group_part_dict = group_part_dict
        self.element_symbols = element_symbols
        self.charges = charges
        self.positions = positions
        if box:
            self.box = box

    def generate_geometry(self):
        """
        Create the geometry string to feed into the Psi4 calculation.
        """
        # Use unrestricted Kohn-Sham if molecule is not in the singlet
        # state.
        qm_atom_list = self._group_part_dict["qm_atom"]
        qm_centroid = np.array([sum([self._positions[atom][k] 
                                for atom in qm_atom_list])
                                / len(qm_atom_list) for k in range(3)])
        if self.qm_spin > 1:
            psi4.core.set_local_option("SCF", "REFERENCE", "UKS")
        # Add MM charges in QMregion for analytic embedding.
        if "analytic" in self._group_part_dict:
            embedding_list = self._group_part_dict["analytic"]
            charge_field = []
            for residue in embedding_list:
                nth_centroid = [sum([self._positions[atom][k] for atom in residue])
                                / len(residue) for k in range(3)]
                new_centroid = least_mirror_vector(
                    nth_centroid,
                    qm_centroid,
                    self._box,
                )
                displacement = np.array(new_centroid) - np.array(nth_centroid)
                for atom in residue:
                    position = (self._positions[atom]
                                + displacement
                                + qm_centroid) * self.angstrom_to_bohr
                    charge_field.append(
                        [
                            self._charges[atom],
                            [position[0],position[1],position[2]],
                        ]
                    )
            self.charge_field = charge_field
        # Construct geometry string.
        geometrystring = (' \n '
                          + str(self.qm_charge) + " " + str(self.qm_spin) + " \n"
                          + " noreorient  \n  " + " nocom  \n  ")
        for atom in qm_atom_list:
            position = self._positions[atom]
            geometrystring = (geometrystring + " " 
                              + str(self._element_symbols[atom]) + " " 
                              + str(position[0]) + " " 
                              + str(position[1]) + " " 
                              + str(position[2]) + " \n")
        geometrystring = geometrystring + ' symmetry c1 \n '
        # now create Psi4 geometry object.
        self.geometry = psi4.geometry(geometrystring)

    def build_ref_quadrature(self):
        """
        Build a reference quadrature for the external potential grid.

        Returns
        -------
        ref_quadrature: NumPy Array object
            A reference quadrature constructed from the current geometry
            of the QM subsystem.
        """
        sup_func = psi4.driver.dft.build_superfunctional(
            self.functional,
            True,
        )[0]
        basis = psi4.core.BasisSet.build(
            self.geometry,
            "ORBITAL",
            self.basis_set,
        )
        ref_v = psi4.core.VBase.build(basis, sup_func, "RV")
        ref_v.initialize()
        ref_quadrature = ref_v.get_np_xyzw()
        return ref_quadrature

    def compute_energy(self, quad_extd_pot=None, nuc_extd_pot=None):
        """
        Perform the Psi4 energy calculation.

        Parameters
        ----------
        quad_extd_pot: NumPy Array object, Optional, default=None
            The extended potential evaluated on the Psi4 quadrature grid.
        nuc_extd_pot: NumPy Array object, Optional, default=None
            The extended potential evaluated at the nuclear coordinates.
        """
        options = {}
        if self.qmmm_pme:
            psi4.core.set_local_option("SCF","PME", self.qmmm_pme)
            options["quad_extd_pot"] = quad_extd_pot
            options["nuc_extd_pot"] = nuc_extd_pot
        if self.charge_field:
            options["external_potentials"] = self.charge_field
        (self.energy, self.wfn) = psi4.energy(
                self.functional,
                return_wfn=True,
                **options,
        )

    def compute_forces(self, quad_extd_pot=None, nuc_extd_pot=None, nuc_extd_grad=None):
        """
        Perform the Psi4 force calculation.

        Parameters
        ----------
        quad_extd_pot: NumPy Array object, Optional, default=None
            The extended potential evaluated on the Psi4 quadrature grid.
        nuc_extd_pot: NumPy Array object, Optional, default=None
            The extended potential evaluated at the nuclear coordinates.
        nuc_extd_grad: NumPy Array object, Optional, default=None
            The gradient of the extended potential evaluated at the
            nuclear coordinates.
        """
        options = {}
        if self.qmmm_pme:
            psi4.core.set_local_option("SCF","PME", self.qmmm_pme)
            options["quad_extd_pot"] = quad_extd_pot
            options["nuc_extd_pot"] = nuc_extd_pot
            options["nuc_extd_grad"] = nuc_extd_grad
        if self.charge_field:
            options["external_potentials"] = self.charge_field
        self.forces = -np.asarray(psi4.gradient(
            self.functional,
            return_wfn=False,
            **options,
        ))

    @property
    def box(self):
        """
        The box vectors defining the periodic system.
        """
        return self._box
    
    @box.setter
    def box(self, box):
        self._box = box
        self._bohr_box = [[k*self.angstrom_to_bohr for k in vector]
                          for vector in box]
