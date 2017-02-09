#!/usr/bin/env python
'''
Description:
    Setup file paths for systems that I usually work on.
'''

#Duncan Campbell
#September 2, 2012
#Yale University
#Setup file paths for common systems I work on.


__all__        =['get_system','get_base_path','get_scripts_path','get_results_path','get_rawdata_path','get_cleandata_path']
__copyright__  =["Copyright 2016 Victor Calderon, common utilities in python; modified by Adam Szewciw"]
__author__     =['Adam Szewciw']
__email__      =['adam.o.szewciw@vanderbilt.edu']
__maintainer__ =['Adam Szewciw']


def known_systems():
    return ["bender", "Adams-MacBook-Pro-2"]


def get_system():
    """
    get the name of the system
    """

    import os, sys
    path_to_home = os.getenv("HOME")

    host = os.popen('echo $HOSTNAME').read()
    host = host.split('.')[0]
    host = host.split('\n')[0]

    if host in known_systems(): return host
    else:
        # host = None
        # if os.path.isfile(os.getenv("HOME")+'/.bender'):
        #     host = 'bender'
        # elif os.path.isfile(os.getenv("HOME")+'/.victor-mac-pro'):
        #     host = "Adam's-MacBook-Pro-2"
        # else: raise ValueError('unknown system.')
        # return host
        raise ValueError('unknown system.')



def get_base_path(node=None):
    """
    get the base path to project folder for the system
    """

    if node==None: node = get_system()

    if node=='bender':
        path = '/fs1/szewciw/galaxy_clustering/'
    elif node=='Adams-MacBook-Pro-2':
        path = '/Users/Adam/Codes/SEGUE/galaxy_clustering/'
    else:
        raise ValueError('error: unknown data directory for this enviornment!')

    return path


def get_scripts_path(node=None):
    """
    get the base path to various scripts
    """
    if node==None: node = get_system()

    if node!='bender' and node!='Adams-MacBook-Pro-2':
        raise ValueError('error: unknown code directory for this environment')

    path = get_base_path(node)+'codes/'

    return path


def get_results_path(node=None):
    """
    get the base path to results folder
    """
    if node==None: node = get_system()

    if node!='bender' and node!='Adams-MacBook-Pro-2':
        raise ValueError('error: unknown code directory for this environment')

    path = get_base_path(node)+'results/'

    return path


def get_utils_path(node=None):
    """
    get the base path to utilities
    """
    if node==None: node = get_system()

    if node!='bender' and node!='Adams-MacBook-Pro-2':
        raise ValueError('error: unknown code directory for this environment')

    path = get_base_path(node)+'utils/'

    return path


def get_carmen_path(node=None):
    """
    get the base path to folder containing results of carmen simulations
    """
    if node==None: node = get_system()

    if node!='bender':
        raise ValueError('error: data only exists on bender')

    path = '/fs2/andreas/Carmen/Planck/'

    return path


def get_carmenHD_path(node=None):
    """
    get the base path to folder containing results of carmenHD simulations
    """
    if node==None: node = get_system()

    if node!='bender':
        raise ValueError('error: data only exists on bender')

    path = '/fs2/andreas/Carmelota/Planck/'

    return path


def get_corrfunc_path(node=None):
    """
    get the base path to folder containing results of carmenHD simulations
    """
    if node==None: node = get_system()

    if node!='bender':
        raise ValueError('error: just run corrfunc on bender')

    path = '/fs1/szewciw/Corrfunc/'

    return path


def get_pix_file(node=None):
    """
    get the base path to the pixfile
    """
    if node==None: node = get_system()

    if node!='bender':
        raise ValueError('error: just run corrfunc on bender')

    file = '/home/piscioja/SDSSPix/Maps/lss_geometry_north.dr72.pix'

    return file