import matplotlib
matplotlib.use('Agg')
import numpy as np

import pylab as pl

from nose.tools import set_trace

def trig_arr(N, m, n):
    from test_fwderiv import trig_arr as ta
    X,Y,cos_arr = ta(N, m, 0, np.cos)
    X,Y,sin_arr = ta(N, 0, n, np.sin)
    return cos_arr + sin_arr

def test_threshold():
    from features import eigman
    N = 512
    m, n = 10,15 
    arr = trig_arr(N, m, n)
    thresh = eigman(arr)
