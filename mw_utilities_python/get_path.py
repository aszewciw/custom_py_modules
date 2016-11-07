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
        path = '/fs1/szewciw/MW_structure/'
    elif node=='Adams-MacBook-Pro-2':
        path = '/Users/Adam/Codes/SEGUE/MW_structure/'
    else:
        raise ValueError('error: unknown data directory for this enviorment!')

    return path


def get_scripts_path(node=None):
    """
    get the base path to various scripts
    """
    if node==None: node = get_system()

    if node!='bender' and node!='Adams-MacBook-Pro-2':
        raise ValueError('error: unknown code directory for this environment')

    path = get_base_path(node)+'scripts/'

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


def get_rawdata_path(node=None):
    """
    get the base path to raw data
    """
    if node==None: node = get_system()

    if node!='bender' and node!='Adams-MacBook-Pro-2':
        raise ValueError('error: unknown code directory for this environment')

    path = get_base_path(node)+'data_segue_raw/'

    return path


def get_cleandata_path(node=None):
    """
    get the base path to cleaned data folder
    """
    if node==None: node = get_system()

    if node!='bender' and node!='Adams-MacBook-Pro-2':
        raise ValueError('error: unknown code directory for this environment')

    path = get_base_path(node)+'data_segue_gdwarfs_cln/'

    return path