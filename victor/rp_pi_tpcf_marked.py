# -*- coding: utf-8 -*-
from __future__ import (absolute_import, division, print_function,
                        unicode_literals)
"""
Calculate the redshift space two point correltion function
"""
####import modules########################################################################
import sys
import numpy as np
from math import pi, gamma
from halotools.mock_observables.clustering_helpers import *
from halotools.mock_observables.tpcf_estimators import *
from halotools.mock_observables.pair_counters import xy_z_marked_npairs, xy_z_npairs
from halotools.mock_observables.pair_counters.marked_double_tree_helpers import *
# from halotools.mock_observables import rp_pi_tpcf
##########################################################################################


__all__=['rp_pi_tpcf_marked', 'wp_marked']
__author__ = ['Victor Calderon']
__description__ = ['Used Duncan Campbell code and added the marked functionality']


np.seterr(divide='ignore', invalid='ignore') #ignore divide by zero in e.g. DD/RR


def rp_pi_tpcf_marked(sample1, rp_bins, pi_bins, sample2=None, randoms=None, 
    period=None, do_auto=True, do_cross=True, estimator='Natural', num_threads=1, 
    max_sample_size=int(1e6), approx_cell1_size = None, 
    approx_cell2_size = None, approx_cellran_size = None, wfunc=1, 
    weights1=None, weights2=None, return_rand=True, use_rand=None):
    """ 
    Calculate the redshift space correlation function, :math:`\\xi(r_{p}, \\pi)`
    
    Calclulate the correlation function as a function of pair seperation perpendicular to 
    the line-of-sight(LOS) and parallel to the LOS.
    
    The first two dimensions (x, y) define the plane for perpendicular distances. 
    The third dimension (z) is used for parallel distances.  i.e. x,y positions are on 
    the plane of the sky, and z is the radial distance coordinate.  This is the 'distant 
    observer' approximation.
    
    Example calls to this function appear in the documentation below. 
    See the :ref:`mock_obs_pos_formatting` documentation page for 
    instructions on how to transform your coordinate position arrays into the 
    format accepted by the ``sample1`` and ``sample2`` arguments.   
    For thorough documentation of all features, see :ref:`rp_pi_tpcf_usage_tutorial`. 
    
    Parameters 
    ----------
    sample1 : array_like
        Npts x 3 numpy array containing 3-D positions of points.
    
    rp_bins : array_like
        array of boundaries defining the radial bins perpendicular to the LOS in which 
        pairs are counted.
    
    pi_bins : array_like
        array of boundaries defining the p radial bins parallel to the LOS in which 
        pairs are counted.
    
    sample2 : array_like, optional
        Npts x 3 numpy array containing 3-D positions of points.
    
    randoms : array_like, optional
        Nran x 3 numpy array containing 3-D positions of points.  If no ``randoms`` are 
        provided analytic randoms are used (only valid for periodic boundary conditions).
    
    period : array_like, optional
        Length-3 array defining axis-aligned periodic boundary conditions. If only 
        one number, Lbox, is specified, period is assumed to be [Lbox]*3.
    
    estimator : string, optional
        options: 'Natural', 'Davis-Peebles', 'Hewett' , 'Hamilton', 'Landy-Szalay'
    
    do_auto : boolean, optional
        do auto-correlation?
    
    do_cross : boolean, optional
        do cross-correlation?
    
    num_threads : int, optional
        number of threads to use in calculation. Default is 1. A string 'max' may be used
        to indicate that the pair counters should use all available cores on the machine.
    
    max_sample_size : int, optional
        Defines maximum size of the sample that will be passed to the pair counter. 
        If sample size exceeds ``max_sample_size``, the sample will be randomly down-sampled 
        such that the subsample is equal to ``max_sample_size``.
    
    approx_cell1_size : array_like, optional 
        Length-3 array serving as a guess for the optimal manner by which 
        the `~halotools.mock_observables.pair_counters.FlatRectanguloidDoubleTree` 
        will apportion the sample1 points into subvolumes of the simulation box. 
        The optimum choice unavoidably depends on the specs of your machine. 
        Default choice is to use [max(rp_bins),max(rp_bins),max(pi_bins)] in each 
        dimension, which will return reasonable result performance for most use-cases. 
        Performance can vary sensitively with this parameter, so it is highly 
        recommended that you experiment with this parameter when carrying out  
        performance-critical calculations. 

    approx_cell2_size : array_like, optional 
        Analogous to ``approx_cell1_size``, but for ``sample2``.  See comments for 
        ``approx_cell1_size`` for details. 
    
    approx_cellran_size : array_like, optional 
        Analogous to ``approx_cell1_size``, but for ``randoms``.  See comments for 
        ``approx_cell1_size`` for details. 

    Returns 
    -------
    correlation_function(s) : numpy.ndarray
        *len(rp_bins)-1* by *len(pi_bins)-1* ndarray containing the correlation function 
        :math:`\\xi(r_p, \\pi)` computed in each of the bins defined by input ``rp_bins``
        and ``pi_bins``.
        
        .. math::
            1 + \\xi(r_{p},\\pi) = \\mathrm{DD}r_{p},\\pi) / \\mathrm{RR}r_{p},\\pi)
            
        if ``estimator`` is set to 'Natural', where  :math:`\\mathrm{DD}(r_{p},\\pi)`
        is calculated by the pair counter, and :math:`\\mathrm{RR}(r_{p},\\pi)` is counted 
        internally using "analytic randoms" if ``randoms`` is set to None 
        (see notes for further details).
        
        If ``sample2`` is passed as input (and not exactly the same as ``sample1``), 
        three arrays of shape *len(rp_bins)-1* by *len(pi_bins)-1* are returned:
        
        .. math::
            \\xi_{11}(r_{p},\\pi), \\xi_{12}(r_{p},\\pi), \\xi_{22}(r_{p},\\pi),
        
        the autocorrelation of ``sample1``, the cross-correlation between ``sample1`` and 
        ``sample2``, and the autocorrelation of ``sample2``, respectively. If 
        ``do_auto`` or ``do_cross`` is set to False, the appropriate result(s) are 
        returned.
    
    Notes
    -----
    Pairs are counted using 
    `~halotools.mock_observables.pair_counters.xy_z_npairs`.  This pair 
    counter is optimized to work on points distributed in a rectangular cuboid volume, 
    e.g. a simulation box.  This optimization restricts this function to work on 3-D 
    point distributions.
    
    If the points are distributed in a continuous "periodic box", then ``randoms`` are not 
    necessary, as the geometry is very simple, and the monte carlo integration that 
    randoms are used for in complex geometries can be done analytically.
    
    If the ``period`` argument is passed in, all points' ith coordinate 
    must be between 0 and period[i].
    
    Examples
    --------
    For demonstration purposes we create a randomly distributed set of points within a 
    periodic unit cube. 
    
    >>> Npts = 1000
    >>> Lbox = 1.0
    >>> period = np.array([Lbox,Lbox,Lbox])
    
    >>> x = np.random.random(Npts)
    >>> y = np.random.random(Npts)
    >>> z = np.random.random(Npts)
    
    We transform our *x, y, z* points into the array shape used by the pair-counter by 
    taking the transpose of the result of `numpy.vstack`. This boilerplate transformation 
    is used throughout the `~halotools.mock_observables` sub-package:
    
    >>> coords = np.vstack((x,y,z)).T
    
    >>> rp_bins = np.logspace(-2,-1,10)
    >>> pi_bins = np.logspace(-2,-1,10)
    >>> xi = rp_pi_tpcf(coords, rp_bins, pi_bins, period=period)
    
    """
    
    function_args = [sample1, rp_bins, pi_bins, sample2, randoms, period, do_auto,\
                     do_cross, estimator, num_threads, max_sample_size,\
                     approx_cell1_size, approx_cell2_size, approx_cellran_size]
    
    sample1, rp_bins, pi_bins, sample2, randoms, period, do_auto, do_cross, num_threads,\
        _sample1_is_sample2, PBCs = _rp_pi_tpcf_process_args(*function_args)
    weights1, weights2 = (
        _marked_npairs_process_weights(sample1, sample2, 
            weights1, weights2, wfunc))
    
    def random_counts(sample1, sample2, randoms, rp_bins, pi_bins, period,\
                      PBCs, num_threads, do_RR, do_DR, _sample1_is_sample2,\
                      approx_cell1_size, approx_cell2_size, approx_cellran_size):
        """
        Count random pairs.  There are two high level branches:
            1. w/ or wo/ PBCs and randoms.
            2. PBCs and analytical randoms
        There are also logical bits to do RR and DR pair counts, as not all estimators
        need one or the other, and not doing these can save a lot of calculation.
        
        Analytical counts are N**2*dv*rho, where dv can is the volume of the spherical 
        shells, which is the correct volume to use for a continious cubic volume with PBCs
        """
        
        def cylinder_volume(R,h):
            """
            Calculate the volume of a cylinder(s), used for the analytical randoms.
            """
            return pi*np.outer(R**2.0,h)
        
        #No PBCs, randomms must have been provided.
        if randoms is not None:
            if do_RR==True:
                RR = xy_z_npairs(randoms, randoms, rp_bins, pi_bins, period=period,\
                                 num_threads=num_threads,\
                                 approx_cell1_size=approx_cellran_size,\
                                 approx_cell2_size=approx_cellran_size)
                RR = np.diff(np.diff(RR,axis=0),axis=1)
            else: RR=None
            if do_DR==True:
                D1R = xy_z_npairs(sample1, randoms, rp_bins, pi_bins, period=period,\
                                  num_threads=num_threads,\
                                  approx_cell1_size=approx_cell1_size,\
                                  approx_cell2_size=approx_cellran_size)
                D1R = np.diff(np.diff(D1R,axis=0),axis=1)
            else: D1R=None
            if _sample1_is_sample2: #calculating the cross-correlation
                D2R = None
            else:
                if do_DR==True:
                    D2R = xy_z_npairs(sample2, randoms, rp_bins, pi_bins, period=period,\
                                      num_threads=num_threads,\
                                      approx_cell1_size=approx_cell2_size,\
                                      approx_cell2_size=approx_cellran_size)
                    D2R = np.diff(np.diff(D2R,axis=0),axis=1)
                else: D2R=None
            
            return D1R, D2R, RR
        #PBCs and no randoms--calculate randoms analytically.
        elif randoms is None:
            
            #set the number of randoms equal to the number of points in sample1
            NR = len(sample1)
            
            #do volume calculations
            v = cylinder_volume(rp_bins,2.0*pi_bins) #volume of spheres
            dv = np.diff(np.diff(v, axis=0),axis=1) #volume of annuli
            global_volume = period.prod()
            
            #calculate randoms for sample1
            N1 = np.shape(sample1)[0]
            rho1 = N1/global_volume
            D1R = (N1)*(dv*rho1) #read note about pair counter
            
            #calculate randoms for sample2
            N2 = np.shape(sample2)[0]
            rho2 = N2/global_volume
            D2R = N2*(dv*rho2) #read note about pair counter
                
            #calculate the random-random pairs.
            rhor = NR**2/global_volume
            RR = (dv*rhor) #RR is only the RR for the cross-correlation.
            
            return D1R, D2R, RR
    
    def pair_counts(sample1, sample2, rp_bins, pi_bins, period,\
                    N_thread, do_auto, do_cross, _sample1_is_sample2,\
                    approx_cell1_size, approx_cell2_size, marks1, marks2,
                    wfunc):
        """
        Count data pairs.
        """
        D1D1 = xy_z_marked_npairs(sample1, sample1, rp_bins, pi_bins, period=period,\
                           num_threads=num_threads,\
                           approx_cell1_size=approx_cell1_size,\
                           approx_cell2_size=approx_cell1_size,
                           weights1=marks1, weights2=marks2,wfunc=wfunc)
        D1D1 = np.diff(np.diff(D1D1,axis=0),axis=1)
        if _sample1_is_sample2:
            D1D2 = D1D1
            D2D2 = D1D1
        else:
            if do_cross==True:
                D1D2 = xy_z_marked_npairs(sample1, sample2, rp_bins, pi_bins, period=period,\
                                  num_threads=num_threads,\
                                  approx_cell1_size=approx_cell1_size,\
                                  approx_cell2_size=approx_cell2_size,\
                                  weights1=marks1, weights2=marks2,wfunc=wfunc)
                D1D2 = np.diff(np.diff(D1D2,axis=0),axis=1)
            else: D1D2=None
            if do_auto==True:
                D2D2 = xy_z_marked_npairs(sample2, sample2, rp_bins, pi_bins, period=period,\
                                   num_threads=num_threads,\
                                   approx_cell1_size=approx_cell2_size,\
                                   approx_cell2_size=approx_cell2_size,\
                                   weights1=marks2, weights2=marks2,wfunc=wfunc)
                D2D2 = np.diff(np.diff(D2D2,axis=0),axis=1)
            else: D2D2=None
        
        return D1D1, D1D2, D2D2
    
    do_DD, do_DR, do_RR = _TP_estimator_requirements(estimator)
              
    # How many points are there (for normalization purposes)?
    N1 = len(sample1)
    N2 = len(sample2)
    if randoms is not None:
        NR = len(randoms)
    else:
        #set the number of randoms equal to the number of points in sample1
        #this is arbitrarily set, but must remain consistent!
        NR = N1
    
    #count pairs
    D1D1,D1D2,D2D2 = pair_counts(sample1, sample2, rp_bins, pi_bins, period,\
                                 num_threads, do_auto, do_cross, _sample1_is_sample2,\
                                 approx_cell1_size, approx_cell2_size, weights1, weights2, wfunc)
    if return_rand:
        D1R, D2R, RR = random_counts(sample1, sample2, randoms, rp_bins, pi_bins, period,\
                                     PBCs, num_threads, do_RR, do_DR, _sample1_is_sample2,\
                                     approx_cell1_size, approx_cell2_size, approx_cellran_size)
    
    if _sample1_is_sample2:
        xi_11 = _TP_estimator(D1D1,D1R,RR,N1,N1,NR,NR,estimator)
        if return_rand:
            return xi_11, [D1R, D2R, RR], [D1D1,D1D2,D2D2]
        else:
            return xi_11
    else:
        if (do_auto==True) & (do_cross==True): 
            xi_11 = _TP_estimator(D1D1,D1R,RR,N1,N1,NR,NR,estimator)
            xi_12 = _TP_estimator(D1D2,D1R,RR,N1,N2,NR,NR,estimator)
            xi_22 = _TP_estimator(D2D2,D2R,RR,N2,N2,NR,NR,estimator)
            return xi_11, xi_12, xi_22
        elif (do_cross==True):
            xi_12 = _TP_estimator(D1D2,D1R,RR,N1,N2,NR,NR,estimator)
            return xi_12
        elif (do_auto==True):
            xi_11 = _TP_estimator(D1D1,D1R,D1R,N1,N1,NR,NR,estimator)
            xi_22 = _TP_estimator(D2D2,D2R,D2R,N2,N2,NR,NR,estimator)
            return xi_11

# -*- coding: utf-8 -*-

"""
functions to calculate clustering statistics, e.g. two point correlation functions.
"""

def wp_marked(sample1, rp_bins, pi_bins, sample2=None, randoms=None, period=None,\
       do_auto=True, do_cross=True, estimator='Natural', num_threads=1,\
       max_sample_size=int(1e6), approx_cell1_size=None, approx_cell2_size=None,\
       approx_cellran_size=None, marks1=None, marks2=None, wfunc=1, \
       return_rand=True):
    """ 
    Calculate the projected two point correlation function, :math:`w_{p}(r_p)`,
    where :math:`r_p` is the seperation perpendicular to the line-of-sight (LOS).
    
    Calculation of :math:`w_{p}(r_p)` requires the user to supply bins in :math:`\\pi`,
    the seperation parallel to the line of sight, and the result will in general depend 
    on both the binning and the maximum :math:`\\pi` seperation integrated over.  See 
    notes for further details.
    
    The first two dimensions define the plane for perpendicular distances.  The third 
    dimension is used for parallel distances.  i.e. x,y positions are on the plane of the
    sky, and z is the redshift coordinate. This is the 'distant observer' approximation.
    
    Example calls to this function appear in the documentation below. 
    See the :ref:`mock_obs_pos_formatting` documentation page for 
    instructions on how to transform your coordinate position arrays into the 
    format accepted by the ``sample1`` and ``sample2`` arguments. 
      
    See also :ref:`galaxy_catalog_analysis_tutorial4`. 

    Parameters 
    ----------
    sample1 : array_like
        Npts x 3 numpy array containing 3-D positions of points. 
    
    rp_bins : array_like
        array of boundaries defining the bins perpendicular to the LOS in which 
        pairs are counted.
    
    pi_bins : array_like
        array of boundaries defining the bins parallel to the LOS in which 
        pairs are counted.
    
    sample2 : array_like, optional
        Npts x 3 numpy array containing 3-D positions of points.
    
    randoms : array_like, optional
        Nran x 3 numpy array containing 3-D positions of points.
    
    period : array_like, optional
        Length-k array defining axis-aligned periodic boundary conditions. If only 
        one number, Lbox, is specified, period is assumed to be [Lbox]*3.
        If none, PBCs are set to infinity.
    
    do_auto : boolean, optional
        do auto-correlation?  Default is True.
    
    do_cross : boolean, optional
        do cross-correlation?  Default is True.
    
    estimator : string, optional
        options: 'Natural', 'Davis-Peebles', 'Hewett' , 'Hamilton', 'Landy-Szalay'
    
    num_threads : int, optional
        number of threads to use in calculation. Default is 1. A string 'max' may be used
        to indicate that the pair counters should use all available cores on the machine.
    
    max_sample_size : int, optional
        Defines maximum size of the sample that will be passed to the pair counter. 
        
        If sample size exceeds max_sample_size, the sample will be randomly down-sampled 
        such that the subsample is equal to max_sample_size.
    
    approx_cell1_size : array_like, optional 
        Length-3 array serving as a guess for the optimal manner by which 
        the `~halotools.mock_observables.pair_counters.FlatRectanguloidDoubleTree` 
        will apportion the sample1 points into subvolumes of the simulation box. 
        The optimum choice unavoidably depends on the specs of your machine. 
        Default choice is to use [max(rp_bins),max(rp_bins),max(pi_bins)] in each 
        dimension, which will return reasonable result performance for most use-cases. 
        Performance can vary sensitively with this parameter, so it is highly 
        recommended that you experiment with this parameter when carrying out  
        performance-critical calculations. 

    approx_cell2_size : array_like, optional 
        Analogous to ``approx_cell1_size``, but for sample2.  See comments for 
        ``approx_cell1_size`` for details. 
    
    approx_cellran_size : array_like, optional 
        Analogous to ``approx_cell1_size``, but for randoms.  See comments for 
        ``approx_cell1_size`` for details. 

    Returns 
    -------
    correlation_function(s) : numpy.array
        *len(rp_bins)-1* length array containing the correlation function :math:`w_p(r_p)` 
        computed in each of the bins defined by input ``rp_bins``.
        
        
        If ``sample2`` is not None (and not exactly the same as ``sample1``), 
        three arrays of length *len(rp_bins)-1* are returned: 
        
        .. math::
            w_{p11}(r_p), \\ w_{p12}(r_p), \\ w_{p22}(r_p),

        the autocorrelation of ``sample1``, the cross-correlation between ``sample1`` 
        and ``sample2``, and the autocorrelation of ``sample2``.  If ``do_auto`` or ``do_cross`` 
        is set to False, the appropriate result(s) is not returned.
    
    Notes
    -----
    The projected correlation function is calculated by integrating the 
    redshift space two point correlation function using 
    `~halotools.mock_observables.rp_pi_tpcf`:
    
    .. math::
        w_p(r_p) = \\int_0^{\\pi_{\\rm max}}\\xi(r_p,\\pi)\\mathrm{d}\\pi
    
    where :math:`\\pi_{\\rm max}` is maximum(``pi_bins``) and :math:`\\xi(r_p,\\pi)` 
    is the redshift space correlation function.
    
    Notice that the results will generally be sensitive to the choice of ``pi_bins``, as
    they indicate where to evalulate :math:`\\xi(r_p,\\pi)` and :math:`\\mathrm{d}\\pi`.
    
    Examples
    --------
    For demonstration purposes we create a randomly distributed set of points within a 
    periodic unit cube. 
    
    >>> Npts = 1000
    >>> Lbox = 1.0
    >>> period = np.array([Lbox,Lbox,Lbox])
    
    >>> x = np.random.random(Npts)
    >>> y = np.random.random(Npts)
    >>> z = np.random.random(Npts)
    
    We transform our *x, y, z* points into the array shape used by the pair-counter by 
    taking the transpose of the result of `numpy.vstack`. This boilerplate transformation 
    is used throughout the `~halotools.mock_observables` sub-package:
    
    >>> coords = np.vstack((x,y,z)).T
    
    >>> rp_bins = np.logspace(-2,-1,10)
    >>> pi_bins = np.logspace(-2,-1,10)
    >>> xi = wp(coords, rp_bins, pi_bins, period=period)
    
    See also 
    --------
    :ref:`galaxy_catalog_analysis_tutorial4`

    """
    
    #process input parameters
    function_args = [sample1, rp_bins, pi_bins, sample2, randoms, period, do_auto,\
                     do_cross, estimator, num_threads, max_sample_size,\
                     approx_cell1_size, approx_cell2_size, approx_cellran_size]
    sample1, rp_bins, pi_bins, sample2, randoms, period, do_auto, do_cross, num_threads,\
        _sample1_is_sample2, PBCs = _rp_pi_tpcf_process_args(*function_args)
    
    if _sample1_is_sample2:
        sample2=None
    
    #pass the arguments into the redshift space TPCF function
    if return_rand:
        result, DD_arr, RR_arr = rp_pi_tpcf_marked(sample1, rp_bins=rp_bins, pi_bins=pi_bins,\
                                     sample2 = sample2, randoms=randoms,\
                                     period = period, do_auto=do_auto, do_cross=do_cross,\
                                     estimator=estimator, num_threads=num_threads,\
                                     max_sample_size=max_sample_size,\
                                     approx_cell1_size=approx_cell1_size,\
                                     approx_cell2_size=approx_cell2_size,\
                                     approx_cellran_size=approx_cellran_size, 
                                     weights1=marks1, weights2=marks2, wfunc=wfunc,\
                                     return_rand=True)
    else:
    
    #integrate the redshift space TPCF to get w_p
    def integrate_2D_xi(x,pi_bins):
        return 2.0*np.sum(x*np.diff(pi_bins), axis=1)

    #return the results.
    if _sample1_is_sample2:
        wp_D1D1 = integrate_2D_xi(result,pi_bins)
        return wp_D1D1
    else:
        if (do_auto==True) & (do_cross==True):
            wp_D1D1 = integrate_2D_xi(result[0],pi_bins)
            wp_D1D2 = integrate_2D_xi(result[1],pi_bins)
            wp_D2D2 = integrate_2D_xi(result[2],pi_bins)
            return wp_D1D1, wp_D1D2, wp_D2D2
        elif (do_auto==True) & (do_cross==False):
            wp_D1D1 = integrate_2D_xi(result[0],pi_bins)  
            wp_D2D2 = integrate_2D_xi(result[1],pi_bins)
            return wp_D1D1, wp_D2D2
        elif (do_auto==False) & (do_cross==True):
            wp_D1D2 = integrate_2D_xi(result,pi_bins)
            return wp_D1D2




