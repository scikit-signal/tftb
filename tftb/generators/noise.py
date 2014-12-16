import numpy as np
from utils import nextpow2
from scipy.signal import hilbert

def noisecu(N):
    if N <= 2:
        noise = (np.random.rand(N,1) - 0.5 + 1j*(np.random.rand(N,1)-0.5)) * \
                np.sqrt(6)
    else:
        noise = np.random.rand(2**nextpow2(N),1) - 0.5
        noise = hilbert(noise)/noise.std()/np.sqrt(2)
        inds = len(noise) - np.arange(N-1,0,step=-1)
        noise = noise[inds]
    return noise


def noisecg(N, a1=None, a2=None):
    """
    Generate analytic complex gaussian noise with mean 0.0 and variance 1.0.

    Parameters
    ----------

    N : int
        Length of the desired output signal.

    a1, a2 : float, optional
        Coefficients of the filter through which the noise is passed.
        If only a1 is provided, the filter is a first order filter:

    Returns
    -------
    noise : 1-D ndarray, shape (N,)

    """
    assert N > 0
    if N <= 2:
        noise = (np.random.randn(N, 1)+1j*np.random.randn(N,1))/np.sqrt(2)
    else:
        noise = np.random.randn(2**nextpow2(N),1)
        noise = hilbert(noise)/noise.std()/np.sqrt(2)
        noise = noise[len(noise) - np.arange(N-1,-1,-1) -1]
    return noise.ravel()

def dopnoise():
    pass


if __name__ == "__main__":
    n = noisecg(128)
