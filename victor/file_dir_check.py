#! /usr/bin/env python

# Victor Calderon
# February 14, 2016
# Vanderbilt University

"""
Checks existence of files, directories, etc.
"""

import os
import sys
import subprocess
import traceback
import numpy as num

def Program_Msg(filename):
	"""
	Program message for 'filename'

	Parameters
	----------
	filename: filetype object
	           It will use the path of the file running 
	           to produce the program message.

	Returns
	-------
	Prg_msg: str
	          String message for given 'filename'
	"""
	Prg_msg = '\n\t\t==> {0} >>> '.format(os.path.basename(filename))
	return Prg_msg

def Index(directory, datatype):
	"""
	Indexes the files in a directory 'directory' with a 
	specific data type.

	Parameters
	----------
	directory: str
	            Absolute path to the folder that is indexed.

	datatype: str
	          Data type of the files to be indexed in the folder.

	Returns
	-------
	file_array: array_like 
	             num.array of indexed files in the folder 'directory' 
	             with specific datatype.

	Examples
	--------
	>>> Index('~/data', '.txt')
	>>> array(['A.txt', 'Z.txt', ...])
	"""
	assert(os.path.exists(directory))
	stack = [directory]
	files = []
	while stack:
		directory = stack.pop()
		for file in os.listdir(directory):
			if file.endswith(datatype):
				fullname = os.path.join(directory, file)
				files.append(fullname)
				if os.path.isdir(fullname) and not os.path.islink(fullname):
					stack.append(fullname)
	files = num.array(files)

	return files

def get_immediate_subdirectories(directory):
	"""
	Immediate subdirectories to this folder

	Parameters
	----------
	directory: str
	            absolute path of the directory to be indexed.

	Returns
	-------
	subdir_arr: array_like
	             Numpy array of subdirectories of parent folder 'directory'.
	"""
	assert(os.path.exists(directory))
	subdir_arr = num.array([name for name in os.listdir(directory) if os.path.isdir(os.path.join(directory, name))])

	return subdir_arr

def Path_Folder(directory):
	"""
	Determines if folder exists. It creates it otherwise.

	Paremeters
	----------
	directory: str
	            Absolute path of the expected directory.
	"""
	if not os.path.exists(directory):
		os.makedirs(directory)

def File_Exists(filename):
	"""
	Determines if file exists. Raises an error otherwise.

	Parameters
	----------
	filename: str
	           Absolute path + filename of the file.

	Exceptions
	----------
	   It traces back error if filename is not found.
	"""
	filedef_msg = Program_Msg(__file__)
	try:
		assert(os.path.isfile(filename))
	except AssertionError:
		print("\n {0} File '{1}' not found!".format(filedef_msg, filename))
		_, _, tb = sys.exc_info()
		traceback.print_tb(tb) # Fixed format
		tb_info = traceback.extract_tb(tb)
		filename, line, func, text = tb_info[-1]
		raise AssertionError ('Error in line{0} in statement {1}'.format(line, text))
		sys.exit(1)

def File_Download_needed(filename_local, filename_web):
    """
    Given a `filename_local`, determine if it exists.
    If not, download the file to specified directory as 
    `filename_local` from `filename_web`.

    Parameters
    ----------
    filename_local: string
        Absolute path + name of the file that resides `locally`.

    filename_web: string
        Web path + name of file that resides `online`.
    """
    message_def = Program_Msg(__file__)
    filename_web_ext = os.path.splitext(filename_web)[1]
    if os.path.isfile(filename_local):
        pass
    else:
        if filename_web_ext=='.gz':
            cmd = 'wget '+ filename_web +' -O '+ filename_local+filename_web_ext
            print('{0} {1}'.format(message_def, cmd))
            subprocess.call(cmd, shell=True)
            cmd = 'gzip -d {0}'.format(filename_local + filename_web_ext)
            print( '{0} {1}'.format(message_def, cmd))
            subprocess.call(cmd, shell=True)
        elif filename_web_ext=='.tar':
            cmd = 'wget '+ filename_web +' -O '+ filename_local+filename_web_ext
            print('{0} {1}'.format(message_def, cmd))
            subprocess.call(cmd, shell=True)
            filename_local_dir = os.path.dirname(filename_local)
            cmd = 'tar zxf '+filename_local+filename_web_ext+' -C '+filename_local_dir
            print( '{0} {1}'.format(message_def, cmd))
            subprocess.call(cmd, shell=True)
        else:
            cmd = 'wget '+filename_web+' -O '+filename_local
            print( '{0} {1}'.format(message_def, cmd))
            subprocess.call(cmd, shell=True)
    File_Exists(filename_local)



