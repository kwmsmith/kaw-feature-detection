import numpy as np
from numpy.linalg import norm

import fwderiv

from nose.tools import set_trace, ok_, eq_

fw_float = fwderiv.fwr_dbl
fw_cpx = fwderiv.fwc_complex_x16


def test_roundtrip():

    def _inner(rarr):
        N = rarr.shape[0]
        rarr_cpy = rarr.copy('F')
        carr = np.empty((N/2+1, N), dtype=fw_cpx, order='F')
        carr = fwderiv.fft_r2c(rarr, carr, N, N)
        carr, rarr = fwderiv.fft_c2r(carr, rarr, N, N)
        rarr /= (N*N)
        ok_(np.allclose(rarr, rarr_cpy, rtol=1e-7, atol=1e-6))

    N = 512
    for i in range(10):
        rarr = np.asfortranarray(np.random.rand(N,N).astype(fw_float))
        yield _inner, rarr

def trig_arr(N, m, n, func):
    X = np.linspace(0.0, m*2*np.pi, N, endpoint=False).reshape(1, N)
    Y = np.linspace(0.0, n*2*np.pi, N, endpoint=False).reshape(N, 1)
    rarr = np.asfortranarray(func(X+Y).astype(fw_float))
    rmax = np.max(np.abs(rarr))
    if rmax: rarr /= rmax
    return X, Y, rarr

class test_fft(object):

    def setup(self):
        self.N = 16

    def test_const(self):
        rarr = np.zeros((self.N,self.N), dtype=fw_float, order='F')
        rarr.fill(10.0)
        rarr_out = rarr.copy('F')
        ok_(np.allclose(0.0, fwderiv.deriv(rarr, rarr_out, axis=1, order=1)))
        ok_(np.allclose(0.0, fwderiv.deriv(rarr, rarr_out, axis=2, order=1)))

    def test_sine(self):
        def _inner(N, m, n):
            X, Y, rarr = trig_arr(N, m, n, np.cos)
            carr = np.empty((N/2+1, N), dtype=fw_cpx, order='F')
            carr = fwderiv.fft_r2c(rarr, carr, N, N)
            carr /= (N * N)
            cabs = np.abs(carr)/np.max(np.abs(carr))
            ok_(np.allclose(cabs[n,m], 1.0))
            carr, rarr = fwderiv.fft_c2r(carr, rarr, N, N)

        self.N = 8
        for m in range(0, self.N/2):
            for n in range(0, self.N/2):
                yield _inner, self.N, m, n

    def test_deriv(self):

        def _inner(N, m, n, axis, func, deriv):
            if axis == 1:
                mult = n
            else:
                mult = m
            X, Y, rarr = trig_arr(N, m, n, func)
            rout = np.empty(rarr.shape, dtype=rarr.dtype, order='F')
            rout = fwderiv.deriv(rarr, rout, axis=axis, order=1)/(N*N)
            rderiv = np.asfortranarray((deriv(X+Y).astype(fw_float))*mult)
            if not (np.allclose(rderiv, rout, rtol=1e-7, atol=1e-6)):
                print np.max(np.abs(rderiv-rout))

        N = 16
        for m in range(0, N/4):
            for n in range(0, N/4):
                yield _inner, N, m, n, 1, np.sin, np.cos
                yield _inner, N, m, n, 2, np.sin, np.cos


class test_gammas(object):

    def setup(self):
        self.N = 128

    def test_gam_scl(self):
        def _inner(N, m, n):
            X, Y, rarr = trig_arr(N, m, n, np.cos)
            u = np.empty(rarr.shape, dtype=rarr.dtype, order='F')
            v = np.empty(rarr.shape, dtype=rarr.dtype, order='F')
            gam_scal = np.empty(rarr.shape, dtype=rarr.dtype, order='F')
            u, v = fwderiv.vel_xy(rarr, u, v)
            u /= (N*N)
            v /= (N*N)
            gam_scal = fwderiv.gamma_scal(u, v, gam_scal)
            gam_scal /= (N*N)
            ok_(np.allclose(gam_scal, 0.0, atol=1e-7, rtol=1e-7), `(N, m, n)`)

        for i in range(self.N/2):
            for j in range(self.N/2):
                _inner(self.N, i, j)
