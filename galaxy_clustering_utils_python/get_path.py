#!/usr/bin/env python
'''
Description:
    Setup file paths for systems that I usually work on.
'''

#Duncan Campbell
#September 2, 2012
#Yale University
#Setup file paths for common systems I work on.


__all__        =['get_system','get_base_path','get_scripts_path','get_adam_path',
                 'get_gillian_path','get_rawdata_path','get_cleandata_path'
                 'get_gsl_lib', 'get_gsl_inc']
__copyright__  =["Copyright 2016 Victor Calderon, common utilities in python; modified by Adam Szewciw"]
__author__     =['Adam Szewciw']
__email__      =['adam.o.szewciw@vanderbilt.edu']
__maintainer__ =['Adam Szewciw']


def known_systems():
    return ["bender", "Adams-MacBook-Pro-2", "stampede", "stampede2"]

def known_login_nodes():
    return ["login1", "login2", "login3", "login4"]

def get_system():
    """
    get the name of the system
    """

    import os, sys
    # path_to_home = os.getenv("HOME")

    host = os.popen('echo $HOSTNAME').read()
    tmp = host.split('.')[0]
    tmp = tmp.split('\n')[0]

    # Check if first term in hostname is a login node
    if tmp in known_login_nodes():
        host = host.split('.')[1]
    else:
        host = tmp

    if host in known_systems(): return host
    else:
        raise ValueError('unknown system.')


def get_base_path(node=None):
    """
    get the base path to project folder for the system
    """
    import os

    if node==None: node = get_system()

    if node=='bender':
        path = '/fs1/szewciw/galaxy_clustering'
    elif node=='Adams-MacBook-Pro-2':
        path = '/Users/Adam/Codes/galaxy_clustering'
    elif node=='stampede' or node=='stampede2':
        path = os.getenv('HOME') + '/galaxy_clustering'
    else:
        raise ValueError('error: unknown data directory for this enviornment!')

    return path


def get_scripts_path(node=None):
    """
    get the base path to various scripts
    """
    if node==None: node = get_system()

    if not (node in known_systems()):
        raise ValueError('error: unknown code directory for this environment')

    path = get_base_path(node)+'/codes'

    return path


def get_adam_path(node=None):
    """
    get the base path to Adam's results folder
    """
    if node==None: node = get_system()

    if not (node in known_systems()):
        raise ValueError('error: unknown code directory for this environment')

    path = get_base_path(node)+'/adam'

    return path

def get_gillian_path(node=None):
    """
    get the base path to Gillian's results folder
    """
    if node==None: node = get_system()

    if not (node in known_systems()):
        raise ValueError('error: unknown code directory for this environment')

    path = get_base_path(node)+'/gillian'

    return path

def get_utils_path(node=None):
    """
    get the base path to utilities
    """
    if node==None: node = get_system()

    if not (node in known_systems()):
        raise ValueError('error: unknown code directory for this environment')

    path = get_base_path(node)+'/utils'

    return path


def get_carmen_path(node=None):
    """
    get the base path to folder containing results of carmen simulations
    """
    if node==None: node = get_system()

    if node!='bender':
        raise ValueError('error: data only exists on bender')

    path = '/fs2/andreas/Carmen/Planck'

    return path


def get_carmenHD_path(node=None):
    """
    get the base path to folder containing results of carmenHD simulations
    """
    if node==None: node = get_system()

    if node!='bender':
        raise ValueError('error: data only exists on bender')

    path = '/fs2/andreas/Carmelota/Planck'

    return path


def get_corrfunc_path(node=None):
    """
    get the base path to folder containing manodeep's corrfunc
    """
    if node==None: node = get_system()

    if node!='bender':
        raise ValueError('error: just run corrfunc on bender')

    path = '/fs1/szewciw/Corrfunc'

    return path


def get_pix_file(node=None):
    """
    get the base path to the pixfile
    """
    if node==None: node = get_system()

    if node!='bender':
        raise ValueError('Pixfile on bender')

    file = '/home/piscioja/SDSSPix/Maps/lss_geometry_north.dr72.pix'

    return file

def get_gsl_lib(node=None):
    """
    get path to this system's gsl lib directory
    """
    import os
    if node==None: node = get_system()

    if not (node in known_systems()):
        raise ValueError('error: unknown code directory for this environment')

    if node=='Adams-MacBook-Pro-2':
        raise ValueError("Don't run this on your mac, dummy")

    if node=='bender':
        path='/usr/local/gsl/latest/nehalem/intel12/nonet/lib'
    elif node=='stampede' or node=='stampede2':
        path=os.environ['TACC_GSL_LIB']

    return path

def get_gsl_inc(node=None):
    """
    get path to this system's gsl include directory
    """
    import os
    if node==None: node = get_system()

    if not (node in known_systems()):
        raise ValueError('error: unknown code directory for this environment')

    if node=='Adams-MacBook-Pro-2':
        raise ValueError("Don't run this on your mac, dummy")

    if node=='bender':
        path='/usr/local/gsl/latest/nehalem/intel12/nonet/include'
    elif node=='stampede' or node=='stampede2':
        path=os.environ['TACC_GSL_INC']

    return path
