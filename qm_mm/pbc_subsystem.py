#! /usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Extended subsystem for QM/MM/PME.

Imports
-------
sys: Standard
numpy: Third Party
scipy: Third Party
system: Local
utils: Local
"""
import sys

import numpy as np
import scipy

from .system import System
from .utils import *


class PBCSubsystem(System):
    """
    The extended environment subsystem for the QM/MM/PME method.

    Parameters
    ----------
    qmmm_pme_gridnumber: int
        The number of grid points to use along each of the box vectors
        of the principal cell during the PME procedure.
    qmmm_pme_alpha: float
        The Gaussian width of the smeared point charges in the Ewald
        summation scheme, in inverse nanometers.  See OpenMM's
        documentation for further discussion.
    group_part_dict: dict of list of list of int
        The indices of particles in the system grouped by a given key.
        Example keys include "qm_atom", "qm_drude", and "analytic".
    charges: NumPy Array object
        The charges of all particles in the system, in proton charge
        units.
    positions: NumPy Array object
        The positions of all particles in the system, in Angstroms.
    box: NumPy Array object
        The box vectors defining the periodic system, in Angstroms.
    """
    
    def __init__(self, 
            qmmm_pme_gridnumber,
            qmmm_pme_alpha,
            group_part_dict,
            charges,
            positions,
            box,
        ):
        System.__init__(self)
        self.qmmm_pme_gridnumber = qmmm_pme_gridnumber
        self.qmmm_pme_alpha = qmmm_pme_alpha
        self.group_part_dict = group_part_dict
        self.charges = charges
        self.positions = positions
        self.box = box

    def build_pme_exclusions(self, ref_quadrature):
        """
        Collect the PME points to which exclusions will be applied.

        The points include the region containing the quadrature grid.

        Parameters
        ----------
        ref_quadrature: NumPy Array object
            A reference quadrature constructed from the geometry of the
            QM subsystem.
        """
        # Create real-space coordinates of the PME grid in Bohr.
        x = np.linspace(0, sum([self._bohr_box[0][i]**2 for i in range(3)])**0.5,
                        self.qmmm_pme_gridnumber, endpoint=False)
        y = np.linspace(0, sum([self._bohr_box[1][i]**2 for i in range(3)])**0.5,
                        self.qmmm_pme_gridnumber, endpoint=False)
        z = np.linspace(0, sum([self._bohr_box[2][i]**2 for i in range(3)])**0.5,
                        self.qmmm_pme_gridnumber, endpoint=False)
        X, Y, Z = np.meshgrid(x, y, z, indexing='ij')
        x = X.flatten()[:,np.newaxis]
        y = Y.flatten()[:,np.newaxis]
        z = Z.flatten()[:,np.newaxis]
        self.pme_xyz = np.concatenate((x,y,z), axis=1)
        # Extract real-space coordinates from the reference quadrature.
        x, y, z, w = ref_quadrature
        quadrature_grid = []
        for i in range(len(x)):
              quadrature_grid.append([x[i], y[i], z[i]])
        quadrature_grid = np.array(quadrature_grid)
        # Project box to reciprocal space.
        inverse_box = np.linalg.inv(self._bohr_box)
        # Project quadrature grid to reciprocal space.
        quadrature_grid_project = self.project_to_pme_grid(
            quadrature_grid,
            inverse_box,
            self.qmmm_pme_gridnumber,
        )
        xi = quadrature_grid_project
        # Perform floor and ceil operations for each dimension.  Each
        # quadrature grid point in reciprocal space should be within a
        # cube, or voxel, of PME grid points.  All 8 grid points
        # defining that voxel should be included in the points to which
        # exclusions are applied because these grid points will be used
        # to interpolate the potential onto the quadrature.
        xx = [np.floor(xi[:,0].reshape((-1,1))),
              np.ceil(xi[:,0].reshape((-1,1)))]
        xy = [np.floor(xi[:,1].reshape((-1,1))),
              np.ceil(xi[:,1].reshape((-1,1)))]
        xz = [np.floor(xi[:,2].reshape((-1,1))),
              np.ceil(xi[:,2].reshape((-1,1)))]
        # Populate all 8 permutations of the above operations.
        xf = []
        for i in range(2):
            for j in range(2):
                for k in range(2):
                    xf.append(np.concatenate((xx[i],xy[j],xz[k]),axis=1))
        # Remove non-unique points.
        xf = np.unique(np.concatenate(tuple(xf),axis=0),axis=0)
        xf[xf==self.qmmm_pme_gridnumber] = 0
        self.pme_exclusion_list = xf

    def apply_pme_exclusions(self, external_grid):
        """
        Apply exlcusions to relevant external potential grid points.

        Parameters
        ----------
        external_grid: list of float
            The external potential grid calculated by OpenMM.

        Returns
        -------
        external_grid: NumPy Array object
            The external potential grid with exclusions applied for the
            potential associated with the QM atoms and analytic
            embedding particles in the principal box, as well as
            real-space embedding corrections if selected.
        """
        external_grid = np.array(external_grid) / self.hartree_to_kjmol
        # Prepare an array of indices which will correspond to the
        # corrected v_exteranl and cast to type int.
        indices = (self.pme_exclusion_list[:,0]*self.qmmm_pme_gridnumber**2
                   + self.pme_exclusion_list[:,1]*self.qmmm_pme_gridnumber
                   + self.pme_exclusion_list[:,2]).astype(np.int)
        positions = []
        charges = []
        for atom in self._group_part_dict["qm_atom"]:
            positions.append(self._positions[atom])
            charges.append(self._charges[atom])
        if "qm_drude" in self._group_part_dict:
            for atom in self._group_part_dict["qm_drude"]:
                positions.append(self._positions[atom])
                charges.append(self._charges[atom])
        if "analytic" in self._group_part_dict:
            for residue in self._group_part_dict["analytic"]:
                for atom in residue:
                    positions.append(self._positions[atom])
                    charges.append(self._charges[atom])
        print("PME Exclusions: {}".format(len(positions)))
        # Cast positions to n_atoms x 1 x 3 array, which allows the
        # pme_xyz coordinates in real space to be broadcast onto the
        # QM region positions later.
        positions = np.array(positions)[:,np.newaxis,:] * self.angstrom_to_bohr
        charges = np.array(charges)
        # Get real-space PME grid coordinates for the excluded indices
        pme_exclusions = self.pme_xyz[indices,:]
        # Get the least-mirror postions between PME positions and the
        # QM region atom positions, which will broadcast to produce an
        # n_atom x n_gridpoints x 3 array.
        dr = least_mirror_array(pme_exclusions, positions, self._bohr_box)
        # Get inverse distance, which will produce an
        # n_atom x n_gridpoint array.
        inv_dr = 1 / np.linalg.norm(dr, axis=2)
        alpha_bohr = self.qmmm_pme_alpha / self.nm_to_bohr
        alpha_dr = alpha_bohr / inv_dr
        grid_temp = charges[np.newaxis,:]*inv_dr.T*scipy.special.erf(alpha_dr).T
        # Get indices for atoms that are too close to grid points and
        # substitute the proper, analytical convergent potential.
        mask = np.where(scipy.special.erf(alpha_dr) <= 1*10**(-6))
        grid_temp[mask[1],mask[0]] = alpha_bohr*charges[mask[0]]*2*np.pi**(-0.5)
        # Subtract exlcuded contributions from the v_external as an
        # n_gridpoints x 1 array.
        external_grid[indices] -= np.sum(grid_temp, axis=1) 
        if "realspace" in self._group_part_dict:
            positions = []
            charges = []
            for residue in self._group_part_dict["realspace"]:
                for atom in residue:
                    positions.append(self._positions[atom])
                    charges.append(self._charges[atom])
            # Cast positions to n_atoms x 1 x 3 array, which will allow
            # the PME exclusion positions in real space to be broadcast
            # onto the embedding region positions.
            positions = np.array(positions)[:,np.newaxis,:]*self.angstrom_to_bohr
            charges = np.array(charges)
            # Get least mirror postions between pmegrid positions and 
            # the embedding atom positions, which will broadcast to 
            # produce an n_atom x n_gridpoints x 3 array.
            dr = least_mirror_array(
                pme_exclusions, 
                positions, 
                self._bohr_box,
            )
            # Get inverse distance, which will produce an
            # n_atom x n_gridpoint array.
            inv_dr = 1 / np.linalg.norm(dr, axis=2)
            alpha_dr = alpha_bohr / inv_dr
            external_grid[indices] += np.sum(
                (charges[np.newaxis,:]*inv_dr.T*scipy.special.erfc(alpha_dr).T),
                axis=1,
            )
        return external_grid

    @staticmethod
    def project_to_pme_grid(points, inverse_box, pme_gridnumber):
        """
        Project points onto a PME grid in reciprocal space.

        This algorithm is identical to that used in method
        'pme_update_grid_index_and_fraction' in OpenMM source code,
        ReferencePME.cpp.

        Parameters
        ----------
        points: NumPy Array object
            The real-space points, in Bohr, to project onto the
            reciprocal-space PME grid.
        inverse_box: NumPy Array object
            The inverse of the box vector defining the principal cell of
            the system, in inverse Bohr.
        pme_gridnumber: int
            The number of grid points along each box vector within the
            system.

        Return
        ------
        scaled_grid_points: NumPy Array object
            The projected points in reciprocal-space, in inverse Bohr.
        """
        scaled_grid_points = np.matmul(points, inverse_box)
        scaled_grid_points = (scaled_grid_points - np.floor(scaled_grid_points))*pme_gridnumber
        scaled_grid_points = np.mod(scaled_grid_points.astype(int), pme_gridnumber) + (scaled_grid_points - scaled_grid_points.astype(int))
        return scaled_grid_points
    
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
