#! /usr/bin/env python

# Victor Calderon
# February 15, 2016
# Vanderbilt University

"""
Set of plotting functions commonly used in my codes.
"""

import matplotlib
matplotlib.use( 'Agg' )
import matplotlib.pyplot as plt
import os
import numpy as num

def Rotating_GIF(output_dir, ax, files, output_gif, prefix='tem'):
    """
    Produces GIF for every rotation of figure in three axes.
    And deleted temporary figures.

    Parameters
    ----------
    output_dir: str
            Location of output images

    ax: axis_like object, matplotlib_like object
            Figure axis to be rotated

    files: array_like
            List of absolute paths to output figures.

    output_gif: str
            Path to output gif figure.
            Creates a gif file with output figures.

    prefix: str
            Prefix of output images
    """
    ## Prefix of output images
    output_files_pref = '{0}/{1}'.format(output_dir,prefix)
    ## Rotations
    for jj in range(0, 361, 15 ):
        ''' Rotation about z-axis'''
        elevation = 0
        azimuth   = jj
        ax.view_init(elev = elevation, azim = azimuth )
        fname     = '{0}_0_{1}.png'.format( output_files_pref, jj)
        plt.savefig(  fname )
        files.append( fname )
    for jj in range(0, 91,10 ):
        ''' Rotation about y-axis '''
        elevation = jj
        azimuth   = 0
        ax.view_init(elev = elevation, azim = azimuth )
        fname     = '{0}_1_{1}.png'.format( output_files_pref, jj)
        plt.savefig(  fname )
        files.append( fname )
    for jj in range(0, 91, 10 ):
        ''' Rotation about z-axis '''
        elevation = 90
        azimuth   = jj
        ax.view_init(elev = elevation, azim = azimuth )
        fname     = '{0}_2_{1}.png'.format( output_files_pref, jj)
        plt.savefig(  fname )
        files.append( fname )
    for jj in range(90, 0, -10 ):
        ''' Rotation about x-axis '''
        elevation = jj
        azimuth   = 90
        ax.view_init(elev = elevation, azim = azimuth )
        fname     = '{0}_3_{1}.png'.format( output_files_pref, jj)
        plt.savefig(  fname )
        files.append( fname )
    for jj in range(90, 181, 10 ):
        ''' Rotation about y-axis '''
        elevation = 0
        azimuth   = jj
        ax.view_init(elev = elevation, azim = azimuth )
        fname     = '{0}_4_{1}.png'.format( output_files_pref, jj)
        plt.savefig(  fname )
        files.append( fname )
    for jj in range(0, 91, 10 ):
        elevation = jj
        azimuth   = 180
        ax.view_init(elev = elevation, azim = azimuth )
        fname     = '{0}_5_{1}.png'.format( output_files_pref, jj)
        plt.savefig(  fname )
        files.append( fname )
    for jj in range(180, 271, 10 ):
        elevation = 90
        azimuth   = jj
        ax.view_init(elev = elevation, azim = azimuth )
        fname     = '{0}_6_{1}.png'.format( output_files_pref, jj)
        plt.savefig(  fname )
        files.append( fname )
    for jj in range(90, 0, -10 ):
        elevation = jj
        azimuth   = 270
        ax.view_init(elev = elevation, azim = azimuth )
        fname     = '{0}_7_{1}.png'.format( output_files_pref, jj)
        plt.savefig(  fname )
        files.append( fname )
    ## Creating GIF file
    delay  = 32
    repeat = True
    loop   = -1 if repeat else 0
    os.system('convert -delay %d -loop %d %s %s' %( delay,loop," ".join(files), \
        output_gif) )
    print('\t\t Output GIF file: {0}'.format(output_gif))
    ## Removing temp files
    for fname in files: os.remove(fname)

def Rotating_GIF_2axis(output_dir, ax1, ax2, files, output_gif, prefix='tem'):
    """
    Produces GIF for every rotation of figure in three axes.
    And deleted temporary figures.

    Parameters
    ----------
    output_dir: str
            Location of output images

    ax1: axis_like object, matplotlib_like object
            Primary Figure axis to be rotated

    ax2: axis_like object, matplotlib_like object
            Secondary Figure axis to be rotated

    files: array_like
            List of absolute paths to output figures.

    output_gif: str
            Path to output gif figure.
            Creates a gif file with output figures.

    prefix: str
            Prefix of output images
    """
    ## Prefix of output images
    output_files_pref = '{0}/{1}'.format(output_dir,prefix)
    ## Rotations
    for jj in range(0, 361, 15 ):
        ''' Rotation about z-axis'''
        elevation = 0
        azimuth   = jj
        ax1.view_init(elev = elevation, azim = azimuth )
        ax2.view_init(elev = elevation, azim = azimuth )
        fname     = '{0}_0_{1}.png'.format( output_files_pref, jj)
        plt.savefig(  fname )
        files.append( fname )
    for jj in range(0, 91,10 ):
        ''' Rotation about y-axis '''
        elevation = jj
        azimuth   = 0
        ax1.view_init(elev = elevation, azim = azimuth )
        ax2.view_init(elev = elevation, azim = azimuth )
        fname     = '{0}_1_{1}.png'.format( output_files_pref, jj)
        plt.savefig(  fname )
        files.append( fname )
    for jj in range(0, 91, 10 ):
        ''' Rotation about z-axis '''
        elevation = 90
        azimuth   = jj
        ax1.view_init(elev = elevation, azim = azimuth )
        ax2.view_init(elev = elevation, azim = azimuth )
        fname     = '{0}_2_{1}.png'.format( output_files_pref, jj)
        plt.savefig(  fname )
        files.append( fname )
    for jj in range(90, 0, -10 ):
        ''' Rotation about x-axis '''
        elevation = jj
        azimuth   = 90
        ax1.view_init(elev = elevation, azim = azimuth )
        ax2.view_init(elev = elevation, azim = azimuth )
        fname     = '{0}_3_{1}.png'.format( output_files_pref, jj)
        plt.savefig(  fname )
        files.append( fname )
    for jj in range(90, 181, 10 ):
        ''' Rotation about y-axis '''
        elevation = 0
        azimuth   = jj
        ax1.view_init(elev = elevation, azim = azimuth )
        ax2.view_init(elev = elevation, azim = azimuth )
        fname     = '{0}_4_{1}.png'.format( output_files_pref, jj)
        plt.savefig(  fname )
        files.append( fname )
    for jj in range(0, 91, 10 ):
        elevation = jj
        azimuth   = 180
        ax1.view_init(elev = elevation, azim = azimuth )
        ax2.view_init(elev = elevation, azim = azimuth )
        fname     = '{0}_5_{1}.png'.format( output_files_pref, jj)
        plt.savefig(  fname )
        files.append( fname )
    for jj in range(180, 271, 10 ):
        elevation = 90
        azimuth   = jj
        ax1.view_init(elev = elevation, azim = azimuth )
        ax2.view_init(elev = elevation, azim = azimuth )
        fname     = '{0}_6_{1}.png'.format( output_files_pref, jj)
        plt.savefig(  fname )
        files.append( fname )
    for jj in range(90, 0, -10 ):
        elevation = jj
        azimuth   = 270
        ax1.view_init(elev = elevation, azim = azimuth )
        ax2.view_init(elev = elevation, azim = azimuth )
        fname     = '{0}_7_{1}.png'.format( output_files_pref, jj)
        plt.savefig(  fname )
        files.append( fname )
    ## Creating GIF file
    delay  = 32
    repeat = True
    loop   = -1 if repeat else 0
    os.system('convert -delay %d -loop %d %s %s' %( delay,loop," ".join(files), \
        output_gif) )
    print('\t\t Output GIF file: {0}'.format(output_gif))
    ## Removing temp files
    for fname in files: os.remove(fname)

def GIF_MOVIE(files, output_gif, delay=60, repeat=True, removef=False):
    """
    Given a list if 'files', it creates a gif file, and deletes temp files.

    Parameters
    ----------
    files: array_like
            List of abs. paths to temporary figures

    output_gif: str
            Absolute path to output gif file.
    """
    loop = -1 if repeat else 0
    os.system('convert -delay %d -loop %d %s %s' %( delay,loop," ".join(files), \
        output_gif) )

    if removef:
        for fname in files: os.remove(fname)

def Med_scatter_plot(x_data, y_data, ax, statfunc=num.median, alpha=.34, mode='perc',\
    perc_opt='data', fill=True, *args):
    """
    Calculates median, scatter or percentiles.
    Optionally: plots the median relation and scatter/percentiles of the data.

    x_data: array_like, Shape (N, ...), one-dimensional array of x and y values.
        Array of x- and y-values for `N` number of data points.

    y_data: array_like, Shape (N,2), two-dimensional array of x and y values.
        Array of x- and y-values for `N` number of data points.

    ax: axis_like object, matplotlib object
        axis to plot the results.

    statfunc: numpy statistical function (default=numpy.median)
        Statistical function to evaluate the data

    alpha: float (default = .34)
        Percentile value to be estimated
        Value between (0,100)

    mode: string
        Type of plot to produce.
        - 'perc': Plots the values between the 50 pm `alpha` about the med, 
                    plus line for statfunc.
        - 'scat': Plots the mean of each bin, along with the scatter of the data.

    perc_opt: string, optional
        Option for how to calculate the percentiles in each bin of data.
        - perc_opt='data': It uses the actual data enclosed by `alpha` about med.
        - perc_opt='numpy': Use the numpy.percentile function to estimate the 
                            percentiles.

    fill: boolean, optional (default=True)
        Option for filling the area of interest.

    args: array_like, optional
        Array of arguments to be passed to matplotlib functions.
    """
    assert(mode=='perc' or mode=='scat')
    assert(alpha>=0. and alpha<=1.)
    assert()
    n_bins = y_data.shape[0]
    stat_sort_ydat=num.sort(y_data)
    if mode=='perc':
        stat_sort_dat_val = statfunc(stat_sort_ydat, axis=1)
        if perc_opt=='data':
            stat_sort_dat_per=num.array([
                [stat_sort_ydat[xx][int(((alpha/2.))*len(stat_sort_ydat[xx]))],
                stat_sort_ydat[xx][int((1.-(alpha/2.))*len(stat_sort_ydat[xx]))]]
                for xx in range(len(stat_sort_ydat))])
        elif perc_opt=='numpy':
            stat_sort_dat_per=num.array([
                [num.percentile(stat_sort_ydat[xx],100.*(alpha/2.)),
                num.percentile(stat_sort_ydat[xx],100.*(1.-(alpha/2.)))]
                for xx in range(len(stat_sort_ydat))])
        if fill:
            ax.plot(x_data, stat_sort_dat_val, *args)
            ax.fill_between(x_data, y1=stat_sort_dat_per.T[0], 
                y2=stat_sort_dat_per.T[1], *args)
        else:
            ax.errorbar(x_data, stat_sort_dat_val, yerr=stat_sort_dat_per, *args)































