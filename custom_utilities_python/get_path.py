#!/usr/bin/env python
'''
Description:
    Setup file paths for systems that I usually work on.
'''

#Duncan Campbell
#September 2, 2012
#Yale University
#Setup file paths for common systems I work on.

__all__        =['get_system','get_base_path','get_data_path','get_output_path','get_plot_path']
__copyright__  =["Copyright 2016 Victor Calderon, common utilities in python"]
__author__     =['Victor Calderon']
__email__      =['victor.calderon@vanderbilt.edu']
__maintainer__ =['Victor Calderon']


def known_systems():
    return ['bender', 'Victors-MacBook-Pro-2']


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
        host = None
        if os.path.isfile(os.getenv("HOME")+'/.bender'): host = 'bender'
        elif os.path.isfile(os.getenv("HOME")+'/.victor-mac-pro'): 
            host = 'Victors-MacBook-Pro-2'
        else: raise ValueError('unknown system.')
        return host

def get_base_path(node=None):
    """
    get the base path for the system
    """
    
    if node==None: node = get_system()

    if node=='bender':
        path = '/fs1/caldervf/CODES/vandy_group_statistics2/'
    elif node=='Victors-MacBook-Pro-2':
        path = '/Users/victor2/Documents/REPOSITORIES/vandy_group_statistics2/'
    else:
        return 'error: unknown data directory for this enviorment!' 

    return path

def get_code_c(node=None):
    """
    get the base path to codes in c for the system
    """
    if node==None: node = get_system()
    if node=='bender': 
        path='/home/caldervf/Codes2/custom_utilities_c/'
    elif node=='Victors-MacBook-Pro-2':
        path='/Users/victor2/Codes/custom_utilities_c/'
    else:
        raise ValueError ('error: unknown code directory for this environment')

    return path


def get_data_path(node=None):
    """
    get the base path to data storage for the system
    """
    
    if node==None: node = get_system()

    if node=='bender':
        path = get_base_path()+'datafiles/'
    elif node=='Victors-MacBook-Pro-2':
        path = get_base_path()+'datafiles/'
    else:
        return 'error: unknown data directory for this enviorment!'

    return path

def get_output_path(node=None):
    """
    get the base path to get_output_path storage for the system
    """
    
    if node==None: node = get_system()

    if node=='bender':
        path = get_base_path()+'processed_data/'
    elif node=='Victors-MacBook-Pro-2':
        path = get_base_path()+'processed_data/'
    else:
        return 'error: unknown data directory for this enviorment!'

    return path


def get_plot_path(node=None):
    """
    get the base path to plot storage for the system
    """
    
    if node==None: node = get_system()

    if node=='bender':
        path = get_base_path()+'plots/'
    elif node=='Victors-MacBook-Pro-2':
        path = get_base_path()+'plots/'
    else:
        return 'error: unknown data directory for this environment!'

    return path


