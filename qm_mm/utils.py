#! /usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Helper functions accessed by multiple classes.
"""
import dask.array as da
import numpy as np


def least_mirror_array(i_array, j_array, box):
    """
    Returns the least mirror coordinates of i_array with respect to
    j_array given a set of box vectors.

    Parameters
    ----------
    i_array: Dask Array object
        Position array.
    j_array: Dask Array object
        Reference array.
    box: list of list of float
        Box vectors from OpenMM.

    Returns
    -------
    r_array: Dask Array object
        Least mirror array of the position array with respect to the
        reference array.
    """
    r_array = i_array - j_array
    A = da.from_array(box[2])[:,np.newaxis,np.newaxis] 
    B = da.transpose(da.floor(r_array[:,:,2]/box[2][2]+0.5))
    C = da.transpose(A * B)
    r_array = r_array - C
    A = da.from_array(box[1])[:,np.newaxis,np.newaxis] 
    B = da.transpose(da.floor(r_array[:,:,1]/box[1][1]+0.5))
    C = da.transpose(A * B)
    r_array = r_array - C
    A = da.from_array(box[0])[:,np.newaxis,np.newaxis] 
    B = da.transpose(da.floor(r_array[:,:,0]/box[0][0]+0.5))
    C = da.transpose(A * B)
    r_array = r_array - C
    return r_array


def least_mirror_vector(i_vector, j_vector, box):
    """
    Returns the least mirror coordinates of i_vector with respect to
    j_vector given a set of box vectors.

    Parameters
    ----------
    i_vector: NumPy array
        Position vector.
    j_vector: NumPy array
        Reference vector.
    box: list of list of float
        Box vectors from OpenMM.

    Returns
    -------
    r_vector: NumPy array
        Least mirror vector of the position vector with respect to the
        reference vector.
    """
    r_vector = [i_vector[k] - j_vector[k] for k in range(3)]
    r_vector -= box[2] * np.floor(r_vector[2]/box[2][2] + 0.5)
    r_vector -= box[1] * np.floor(r_vector[1]/box[1][1] + 0.5)
    r_vector -= box[0] * np.floor(r_vector[0]/box[0][0] + 0.5)
    return r_vector
