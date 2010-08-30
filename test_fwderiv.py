import numpy as np
from numpy.linalg import norm

from fwderiv import fwderiv

from nose.tools import set_trace, ok_

def test_roundtrip():

    def _inner(rarr):
        N = rarr.shape[0]
        rarr_cpy = rarr.copy('F')
        carr = np.empty((N/2+1, N), dtype=fwderiv.fwc_complex, order='F')
        carr = fwderiv.fft_r2c(rarr, carr, N, N)
        carr, rarr = fwderiv.fft_c2r(carr, rarr, N, N)
        rarr /= (N*N)
        ok_(np.allclose(rarr, rarr_cpy, rtol=1e-5, atol=1e-4))

    N = 512
    for i in range(100):
        rarr = np.asfortranarray(np.random.rand(N,N).astype(fwderiv.fwr_real))
        rarr *= 100.0
        yield _inner, rarr
