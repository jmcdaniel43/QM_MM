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
import itertools
import sys

import dask.array as da
import dask.array.linalg as la
import numpy as np
import scipy.special
import scipy.interpolate

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
    n_threads: int, Optional, default=1
        Number of threads across which to parallelize the QM calculation.
    """
    
    def __init__(self, 
            qmmm_pme_gridnumber,
            qmmm_pme_alpha,
            group_part_dict,
            charges,
            positions,
            box,
            n_threads=1,
        ):
        System.__init__(self)
        self.qmmm_pme_gridnumber = qmmm_pme_gridnumber
        self.qmmm_pme_alpha = qmmm_pme_alpha
        self.group_part_dict = group_part_dict
        self.charges = charges
        self.positions = positions
        self.box = box
        self.n_threads = n_threads
        
    def build_extd_pot(self, ref_quadrature, extd_pot, return_grad=False):
        """
        Build extended potential grids to pass to Psi4.

        Parameters
        ----------
        ref_quadrature: NumPy Array object
            A reference quadrature grid from Psi4 with the appropriate
            geometry of the system.
        extd_pot: NumPy Array object
            The extended potential provided by OpenMM.
        return_grad: bool, Optional, default=False
            Determine whether or not to return the gradient of the
            extended potential at the nuclear coordinates.

        Returns
        -------
        quad_extd_pot: NumPy Array object
            The extended potential evaluated on the Psi4 quadrature grid.
        nuc_extd_pot: NumPy Array object
            The extended potential evaluated at the nuclear coordinates.
        nuc_extd_grad: NumPy Array object, Optional
            The gradient of the extended potential evaluated at the
            nuclear coordinates.
        """
        inverse_box = np.linalg.inv(self._bohr_box)
        xdim = np.array(
            [i for i in range(-1,self.qmmm_pme_gridnumber+1)]
        )
        grid = (xdim, xdim, xdim)
        quadrature_grid = np.transpose(np.array(ref_quadrature))[:,0:3]
        positions = []
        for atom in self._group_part_dict["qm_atom"]:
            positions.append(self._positions[atom])
        nuclei_grid = np.array(positions) * self.angstrom_to_bohr
        # Determining exlcusions
        self.build_pme_exclusions(quadrature_grid, inverse_box)
        # Applying exclusions to extd_pot
        extd_pot = self.apply_pme_exclusions(extd_pot)
        # Preparing arrays for interpolation
        extd_pot_3d = np.reshape(
            extd_pot,
            (self.qmmm_pme_gridnumber, self.qmmm_pme_gridnumber, self.qmmm_pme_gridnumber),
        )
        extd_pot_3d = np.pad(extd_pot_3d, 1, mode="wrap")
        # performing interpolations
        quad_extd_pot = self.interp_extd_pot(
            quadrature_grid,
            grid,
            extd_pot_3d,
            inverse_box,
        )
        nuc_extd_pot = self.interp_extd_pot(
            nuclei_grid,
            grid,
            extd_pot_3d,
            inverse_box,
        )
        if return_grad:
            nuc_extd_grad = self.interp_extd_grad(
                nuclei_grid,
                grid,
                extd_pot_3d,
                inverse_box,
            )
            return quad_extd_pot, nuc_extd_pot, nuc_extd_grad
        else:
            return quad_extd_pot, nuc_extd_pot

    def build_pme_exclusions(self, points, inverse_box):
        """
        Collect the PME points to which exclusions will be applied.

        The points include the region containing the quadrature grid.

        Parameters
        ----------
        points: NumPy Array object
            The grid points to be excluded.
        inverse_box: NumPy Array object
            The inverse box vectors, in inverse Bohr.
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
        # Project quadrature grid to reciprocal space.
        points_project = self.project_to_pme_grid(
            points,
            inverse_box,
        )
        xi = points_project
        indices = np.floor(xi).T
        edges = list(itertools.product(*[[i, i + 1] for i in indices]))
        edges = [np.stack(x,axis=-1) for x in edges]
        xf = np.unique(np.concatenate(tuple(edges),axis=0),axis=0)
        xf[xf==self.qmmm_pme_gridnumber] = 0
        self.pme_exclusion_list = xf

    def apply_pme_exclusions(self, external_grid):
        """
        Apply exlcusions to relevant external potential grid points.

        Parameters
        ----------
        external_grid: NumPy Array object
            The external potential grid calculated by OpenMM.

        Returns
        -------
        external_grid: NumPy Array object
            The external potential grid with exclusions applied for the
            potential associated with the QM atoms and analytic
            embedding particles in the principal box, as well as
            real-space embedding corrections if selected.
        """
        external_grid = da.from_array(external_grid, chunks=(len(external_grid)//self.n_threads,)) / self.hartree_to_kjmol
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
        # Cast positions to n_atoms x 1 x 3 array, which allows the
        # pme_xyz coordinates in real space to be broadcast onto the
        # QM region positions later.
        positions = da.from_array(positions, chunks=(len(positions)//self.n_threads, 3))[:,np.newaxis,:] * self.angstrom_to_bohr
        charges = da.from_array(charges, chunks=(len(charges)//self.n_threads,))
        # Get real-space PME grid coordinates for the excluded indices
        pme_exclusions = da.from_array(self.pme_xyz[indices,:])
        # Get the least-mirror postions between PME positions and the
        # QM region atom positions, which will broadcast to produce an
        # n_atom x n_gridpoints x 3 array.
        dr = least_mirror_array(pme_exclusions, positions, self._bohr_box)
        # Get inverse distance, which will produce an
        # n_atom x n_gridpoint array.
        inv_dr = 1 / la.norm(dr, axis=2)
        alpha_bohr = self.qmmm_pme_alpha / self.nm_to_bohr
        alpha_dr = alpha_bohr / inv_dr
        grid_temp = charges[np.newaxis,:] * da.transpose(inv_dr) * da.transpose(scipy.special.erf(alpha_dr))
        # Get indices for atoms that are too close to grid points and
        # substitute the proper, analytical convergent potential.
        mask = da.where(scipy.special.erf(alpha_dr) <= 1*10**(-6), 1, 0)
        #print(mask.shape,flush=True)
        #grid_temp[mask] = alpha_bohr*charges[mask[0,:]]*2*np.pi**(-0.5)
        # Subtract exlcuded contributions from the v_external as an
        # n_gridpoints x 1 array.
        external_grid[indices] -= grid_temp.sum(axis=1)
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
            positions = da.from_array(positions, chunks=(len(positions)//self.n_threads, 3))[:,np.newaxis,:] * self.angstrom_to_bohr
            charges = da.from_array(charges, chunks=(len(charges)//self.n_threads,))
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
            inv_dr = 1 / la.norm(dr, axis=2)
            alpha_dr = alpha_bohr / inv_dr
            grid_temp = charges[np.newaxis,:] * da.transpose(inv_dr) * da.transpose(scipy.special.erfc(alpha_dr))
            external_grid[indices] += grid_temp.sum(axis=1)
        print(f"external_grid chunk: {external_grid.chunks}")
        external_grid = external_grid.compute(scheduler="threads", num_workers=self.n_threads)
        return external_grid

    def interp_extd_pot(self, points, grid, extd_pot_3d, inverse_box):
        """
        Builds the extended potential to pass to the quadrature.
        
        Parameters
        ----------
        points: NumPy Array object
            The grid points to be interpolated into the external
            potential grid.
        grid: tuple of NumPy Array object
            Tuple of grid axis values.
        extd_pot_3d: NumPy Array object
            The 3d extended potential grid.
        inverse_box: NumPy Array object
            The inverse of the box vectors, in inverse Bohr.

        Returns
        -------
        interp_points: NumPy Array object
            The interpolated extended potential values at the points.
        """
        points_project = self.project_to_pme_grid(points, inverse_box)
        interp_points = scipy.interpolate.interpn(
            grid,
            extd_pot_3d,
            points_project,
            method='linear',
        )
        return interp_points

    def interp_extd_grad(self, points, grid, extd_pot_3d, inverse_box):
        """
        Create the chain rule for the extended potential on the nuclei.

        Parameters
        ----------
        points: NumPy Array object
            The grid points to be interpolated into the external
            potential grid.
        grid: tuple of NumPy Array object
            Tuple of grid axis values.
        extd_pot_3d: NumPy Array object
            The 3d extended potential grid.
        inverse_box: NumPy Array object
            The inverse of the box vectors, in inverse Bohr.

        Returns
        -------
        extd_grad: NumPy Array object
            The gradient of the interpolated extended potential values
            at the points.
        """
        points_project = self.project_to_pme_grid(
            points,
            inverse_box,
        )
        # This code is largely based on
        # scipy.interpolate.RegularGridInterpolator._evaluate linear.
        interp_function = scipy.interpolate.RegularGridInterpolator(
            grid,
            extd_pot_3d,
        )
        indices, norm_dist, out_of_bounds = interp_function._find_indices(
            points_project.T,
        )
        edges = list(itertools.product(*[[i, i + 1] for i in indices]))
        extd_x=0
        extd_y=0
        extd_z=0
        for edge_indices in edges:
            weight_x = 1
            weight_y = 1
            weight_z = 1
            for j, (ei, i, yi) in enumerate(zip(edge_indices, indices, norm_dist)):
                if j == 0:
                    weight_x *= np.where(ei == i, -1.0, 1.0)
                    weight_z *= np.where(ei == i, 1 - yi, yi)
                    weight_y *= np.where(ei == i, 1 - yi, yi)
                if j == 1:
                    weight_y *= np.where(ei == i, -1.0, 1.0)
                    weight_x *= np.where(ei == i, 1 - yi, yi)
                    weight_z *= np.where(ei == i, 1 - yi, yi)
                if j == 2:
                    weight_z *= np.where(ei == i, -1.0, 1.0)
                    weight_y *= np.where(ei == i, 1 - yi, yi)
                    weight_x *= np.where(ei == i, 1 - yi, yi)
            extd_x += np.array(interp_function.values[edge_indices]) * weight_x
            extd_y += np.array(interp_function.values[edge_indices]) * weight_y
            extd_z += np.array(interp_function.values[edge_indices]) * weight_z
        extd_du = np.concatenate(
            (extd_x.reshape((-1,1)), extd_y.reshape((-1,1)), extd_z.reshape((-1,1))),
            axis=1,
        )
        extd_grad = self.qmmm_pme_gridnumber * (extd_du @ inverse_box)
        return extd_grad

    def project_to_pme_grid(self, points, inverse_box):
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

        Returns
        -------
        scaled_grid_points: NumPy Array object
            The projected points in reciprocal-space, in inverse Bohr.
        """
        scaled_grid_points = np.matmul(points, inverse_box)
        scaled_grid_points = ((scaled_grid_points
                               - np.floor(scaled_grid_points))
                              * self.qmmm_pme_gridnumber)
        scaled_grid_points = np.mod(
            scaled_grid_points.astype(int),
            self.qmmm_pme_gridnumber,
        ) + (scaled_grid_points - scaled_grid_points.astype(int))
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
