#!/usr/bin/env python
'''
Description:
    Checks for existence of files or directories.
'''


__all__        =['check_file', 'check_dir']
__copyright__  =['Copyright 2017 Adam Szewciw, common utilities in python']
__author__     =['Adam Szewciw']
__email__      =['adam.o.szewciw@vanderbilt.edu']
__maintainer__ =['Adam Szewciw']


def check_file(filepath, kill_program=True):
    '''
    Check for the existence of a file.
        Default     -- Kill external program.
        Alternate   -- Continue running external progam.
    '''
    import sys, os

    if not os.path.isfile(filepath):
        sys.stderr.write('Error: file {} does not exist.\n'.format(filepath))
        if kill_program:
            sys.stderr.write('Exiting...\n')
            sys.exit()


def check_dir(dirpath, mkdir=True):
    '''
    Check if a directory exists.
        Default     -- Make directory if it does not exist.
        Alternate   -- Kill external program.
    '''
    import sys, os

    if not os.path.isdir(dirpath):
        sys.stderr.write('Error: directory {} does not exist.\n'.format(dirpath))
        if mkdir:
            sys.stderr.write('Making directory and continuing...\n')
            cmd = 'mkdir ' + dirpath
            os.system(cmd)
        else:
            sys.stderr.write('Exiting...\n')
            sys.exit()