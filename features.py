

def eigman(arr, CCW=1.0, CW=-1.0, STR=0.0):
    '''
    does thresholding, yielding regions dominated by isotropic scaling,
    differential rotation (vorticity) and anisotropic stretching
    '''
    tmp1 = arr.copy('F')
    tmp2 = arr.copy('F')
    tmp3 = arr.copy('F')
    from fwderiv.fwderiv import gamma_scal, gamma_rotation, gamma_stretch, vel_xy
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
