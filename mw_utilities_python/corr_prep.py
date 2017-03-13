#!/usr/bin/env python
'''
Description:
    Preparing for correlation function calculation.
'''

__all__        =['set_rbins']
__copyright__  =["Copyright 2016 Adam Szewciw, "]
__author__     =['Adam Szewciw']
__email__      =['adam.o.szewciw@vanderbilt.edu']
__maintainer__ =['Adam Szewciw']

class R_Bin:
    pass

#------------------------------------------------------------------------------#

def set_rbins(r_min=0.005, r_max=2.0, nbins=12, log=True, filepath=None):
    """
    Make radial bins files to be passed into correlation function.
    Units are in parsecs. Defaults used by Mao et al. are listed above.
    """
    import math, pickle, sys

    if filepath==None: raise ValueError('Unspecified path.')

    if log==True:
        r_min_log = math.log10(r_min)
        r_max_log = math.log10(r_max)
        dr_log = (r_max_log - r_min_log) / nbins
    else:
        dr = (r_max - r_min) / nbins

    bins_list = []

    for i in range(nbins):
        b = R_Bin()

        if log==True:
            b.r_lower_log  = r_min_log + dr_log * i
            b.r_upper_log  = b.r_lower_log + dr_log
            b.r_middle_log = b.r_lower_log + 0.5 * dr_log

            b.r_lower  = 10.0**b.r_lower_log
            b.r_upper  = 10.0**b.r_upper_log
            b.r_middle = 10.0**b.r_middle_log
            b.dr       = b.r_upper - b.r_lower
        else:
            b.r_lower  = r_min + dr * i
            b.r_upper  = b.r_lower + dr
            b.r_middle = b.r_lower + 0.5 * dr
            b.dr       = dr

        bins_list.append(b)

    # pickle output
    output_filename = filepath + 'rbins.dat'
    output_file     = open(output_filename, 'wb')
    pickle.dump(bins_list, output_file)
    output_file.close()

    # ascii output
    output_filename = filepath + 'rbins.ascii.dat'
    output_file = open(output_filename, 'w')
    # output number of bins first
    output_file.write('{}\n'.format(nbins))
    for b in bins_list:
        output_file.write('{0:.14f}\t{1:.14f}\t{2:.14f}\t{3:.14f}\n'
                          .format(b.r_lower, b.r_upper, b.r_middle, b.dr))
    output_file.close()

    sys.stderr.write('Bins list output to {}\n\n'.format(filepath))


#------------------------------------------------------------------------------#
