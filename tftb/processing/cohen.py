import numpy as np
from utils import tftb_window

def wigner_ville(x, t=None, N=None, trace=False):
    if x.ndim == 2:
        xrow, xcol = x.shape
        assert xcol <= 2
    else:
        x = x.reshape((len(x),1))
        xrow, xcol = x.shape
    
    if t is None:
        t = np.arange(xrow)
    elif t.ndim !=1:
        t = t.ravel()
    if N is None:
        N = xrow

    tcol = t.shape[0]
    tfr = np.zeros((N, tcol),dtype=complex)
    
    for icol in xrange(tcol):
        ti = t[icol]
        a,b,c = ti, xrow-ti-1,np.round(N/2)-1
        taumax = np.min((a,b,c))
        tau = np.arange(-taumax, taumax+1)
        indices = np.remainder(N+tau, N)
        tfr[indices, icol] = x[ti+tau, 0] * np.conj(x[ti-tau, xcol-1])
        tau = np.round(N/2)
        if (ti <= xrow - tau) and (ti >= tau + 1):
            tfr[tau, icol] = 0.5 * (x[ti+tau,0] * np.conj(x[ti-tau, xcol-1])) + \
                                   (x[ti-tau,0] * np.conj(x[ti+tau, xcol-1]))

    tfr = np.fft.fft(tfr,axis=0)
    if xcol == 1:
        tfr = np.real(tfr)
    return tfr

def pseudo_wigner_ville(x, t=None, N=None, h=None, trace=False):
    if x.ndim == 2:
        x = x.ravel()
    xrow = len(x)
    
    if N is None:
        N = xrow   
    
    if t is None:
        t = np.arange(xrow)
    hlength = np.floor(N/4.0) + 1 - np.remainder(np.floor(N/4.0), 2)
    
    if h is None:
        h = tftb_window(hlength)
    if len(h)%2 == 0:
        raise TypeError("The smoothing window should be odd!")
    from IPython.core.debugger import Tracer
    Tracer()()
    
    tcol = len(t)
    hrow = len(h)
    Lh = (hrow-1)/2
    h = h/h[Lh]
    
    tfr = np.zeros((N, tcol))
    for icol in xrange(tcol):
        ti = t[icol]
        taumax = np.min([ti,xrow-ti-1,np.round(N/2)-1,Lh])
        tau = np.arange(taumax,taumax+1)
        indices = np.remainder(N+tau, N)
        tfr[indices, icol] = h[Lh+tau] * x[ti+tau-1] * \
                             np.conjugate(x[ti-tau-1])
        tau = np.round(N/2)
        if (ti<=xrow-tau) and (ti>=tau+1) and (tau<=Lh):
            tfr[tau,icol] = 0.5 * (h[Lh+1+tau] * x[ti+tau] * np.conj(x[ti-tau])) + \
                                  (h[Lh+1-tau] * x[ti-tau] * np.conj(x[ti+tau]))
    tfr = np.fft.fft(tfr,axis=0)
    return np.real(tfr)

    