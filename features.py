#!/usr/bin/env python

import numpy as np
from fwderiv import gamma_scal, gamma_rotation, gamma_stretch, vel_xy, deriv, make_gauss, convolve_rr

last_sigma = None
last_gauss_arr = None

def convolve_gauss(rarr, sigma):
    global last_sigma, last_gauss_arr
    if sigma == last_sigma:
        gauss_arr = last_gauss_arr
    else:
        ga = np.zeros(rarr.shape, dtype=rarr.dtype, order='f')
        gauss_arr = last_gauss_arr = make_gauss(ga, sigma)
    rout = np.empty(rarr.shape, dtype=np.double, order='f')
    return convolve_rr(rarr, gauss_arr, rout)

def hessian(arr):

    Dxx = arr.copy('F')
    Dxy = arr.copy('F')

    Dxx = deriv(arr, Dxx, axis=1, order=1)
    Dxy = deriv(Dxx, Dxy, axis=2, order=1)

    Dxx = deriv(arr, Dxx, axis=1, order=2)

    Dyy = arr.copy('F')
    Dyy = deriv(arr, Dyy, axis=2, order=2)

    hess = Dxx * Dyy - Dxy**2

    return hess

def eigman(arr, CCW=1.0, CW=-1.0, STR=0.0):
    '''
    does thresholding, yielding regions dominated by isotropic scaling,
    differential rotation (vorticity) and anisotropic stretching
    '''
    tmp1 = arr.copy('F')
    tmp2 = arr.copy('F')
    tmp3 = arr.copy('F')
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

def apply_h5(h5file, fieldnames, fcn, suffix, txt, force=False):
    import tables
    dta = tables.openFile(h5file, mode='r+')
    try:
        for fieldname in fieldnames:
            gp = dta.getNode('/%s' % fieldname)
            new_gp_name = '%s_%s' % (fieldname, suffix)
            try:
                new_gp = dta.createGroup('/', new_gp_name, '%s %s' % (fieldname, txt))
            except tables.exceptions.NodeError:
                new_gp = dta.getNode('/%s' % new_gp_name)
            for arr in gp:
                arr_dta = arr.read()
                res = fcn(arr_dta)
                try:
                    dta.removeNode(new_gp, arr.name)
                except tables.exceptions.NoSuchNodeError:
                    pass
                finally:
                    dta.createArray(new_gp, arr.name, res, txt)
    finally:
        dta.close()

# def eigman_h5(h5file, fieldnames):
    # import tables
    # dta = tables.openFile(h5file, mode='r+')
    # try:
        # for fieldname in fieldnames:
            # gp = dta.getNode('/%s' % fieldname)
            # try:
                # em_gp = dta.createGroup('/', '%s_em' % fieldname, '%s eigenvalue manifold' % fieldname)
            # except tables.exceptions.NodeError:
                # em_gp = dta.getNode('/%s_em' % fieldname)
            # for arr in gp:
                # arr_dta = arr.read()
                # dta.createArray(em_gp, arr.name, eigman(arr_dta), 'eigenvalue manifold')
    # finally:
        # dta.close()

if __name__ == '__main__':
    import sys
    from optparse import OptionParser
    usage = "%prog [options] field_name1 [field_name2 ...]"
    parser = OptionParser(usage=usage)
    parser.add_option("-f", "--file", dest="fname",
                        help="hdf5 data file")
    parser.add_option('-t', '--type', dest='type',
                        help='type of feature ext. to use, either "eigman" or "hessian"')
    parser.add_option('--sigma', dest='sigma', type=int, default=-1,
                        help='sigma of convolution gaussian in pixels')
    options, args = parser.parse_args()
    if not options.fname:
        parser.print_help()
        sys.exit(1)
    if options.type not in ('eigman', 'hessian'):
        parser.print_help()
        sys.exit(1)

    tp_dict = {
            'hessian' : dict(fcn=hessian, suffix='hs', txt='hessian'),
            'eigman' : dict(fcn=eigman, suffix='em', txt='eigenvalue manifold'),
            }

    def smoothing(fcn, sigma):
        def _inner(arr):
            # we assume L == 2 pi.
            sig = float(sigma) / arr.shape[0] * 2.0 * np.pi
            return fcn(convolve_gauss(arr, sig))
        return _inner

    if options.sigma > 0:
        for x in tp_dict:
            tp_dict[x]['fcn'] = smoothing(tp_dict[x]['fcn'], options.sigma)

    apply_h5(options.fname, args, **tp_dict[options.type])
