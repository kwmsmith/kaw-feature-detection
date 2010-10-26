#!/usr/bin/env python

import numpy as np

def eigman(arr, CCW=1.0, CW=-1.0, STR=0.0):
    '''
    does thresholding, yielding regions dominated by isotropic scaling,
    differential rotation (vorticity) and anisotropic stretching
    '''
    tmp1 = arr.copy('F')
    tmp2 = arr.copy('F')
    tmp3 = arr.copy('F')
    from fwderiv import gamma_scal, gamma_rotation, gamma_stretch, vel_xy
    u, v = vel_xy(arr, tmp1, tmp2)
    gscal = gamma_scal(u, v, tmp3)
    grot = gscal.copy('F')
    grot = gamma_rotation(u, v, grot)
    gstr = grot.copy('F')
    gstr = gamma_stretch(u, v, gstr)
    norm = gscal**2 + grot**2 + gstr**2
    gscal /= norm
    grot  /= norm
    gstr  /= norm
    thresh = grot.copy('F')
    thresh.fill(-2.0)
    thresh[grot > gstr] = CCW
    thresh[-grot > gstr] = CW
    thresh[gstr > np.abs(grot)] = STR
    return thresh

def eigman_h5(h5file, fieldnames):
    import tables
    dta = tables.openFile(h5file, mode='r+')
    try:
        for fieldname in fieldnames:
            gp = dta.getNode('/%s' % fieldname)
            try:
                em_gp = dta.createGroup('/', '%s_em' % fieldname, '%s eigenvalue manifold' % fieldname)
            except tables.exceptions.NodeError:
                em_gp = dta.getNode('/%s_em' % fieldname)
            for arr in gp:
                arr_dta = arr.read()
                dta.createArray(em_gp, arr.name, eigman(arr_dta), 'eigenvalue manifold')
    finally:
        dta.close()

if __name__ == '__main__':
    import sys
    from optparse import OptionParser
    usage = "%prog [options] field_name1 [field_name2 ...]"
    parser = OptionParser(usage=usage)
    parser.add_option("-f", "--file", dest="fname",
                        help="hdf5 data file")
    options, args = parser.parse_args()
    if not options.fname:
        parser.print_help()
        sys.exit(1)

    eigman_h5(options.fname, args)
