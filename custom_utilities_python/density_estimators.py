#! /usr/bin/env python

# Victor Calderon
# February 16, 2016
# Vanderbilt University

"""
Set of density estimators for k-nearest searches.
"""

import numpy as num
from scipy.spatial import KDTree

def Nth_Nearest_Neighbor_search(data, nth=1, var='dist'):
    """
    Generates array of nearest neighbor densities based on distance.

    Parameters
    ----------
    data: array_like, Shape (N,3), three-dimensional array.
        N-dimensional array of positions (x,y,z) in phase-space for `N` items.

    nth: int
        Number of nth nearest neighbor, up to which to calculate distances.

    var: string, optional (default='dist')
        variable to evaluate
        Options: 'dist' or 'dens'

    Returns
    -------
    nth_var_arr: array_like (default=1, but for coding reasons)
        Numpy array of distances/densities to nth nearest neighbor.
        - if 'var=`dist`' -> Numpy Array of distance to nth nearest neighbor.
        - if 'var=`dens`' -> Array of densities to nth nearest neighbor.
    """
    assert(data.ndim==2 and data.shape[-1]==3)
    assert(type(nth)==float or type(nth)==int)
    assert(nth>=0)
    assert(var=='dist' or var=='dens')
    try: 
        nth=int(int) 
    except: 
        raise ValueError ("'{0}' is not a number".format(nth))
    # KDTree calcs
    data_tree = KDTree(data)
    KD_query  = data_tree.query(data, nth+1)
    dist_arr  = KD_query[0]
    if var=='dist':
        nth_dist_arr = dist_arr.copy() if nth ==0 else dist_arr[:,nth]
        nth_var_arr  = nth_dist_arr.copy()
    elif var=='dens':
        nth_dist_arr = dist_arr.copy() if nth ==0 else dist_arr[:,nth]
        nth_rho_arr  = dist_arr.copy() if nth ==0 else nth/(4.*num.pi*(nth_dist_arr**3)/3.)
        nth_var_arr  = nth_rho_arr.copy()

    return nth_var_arr
