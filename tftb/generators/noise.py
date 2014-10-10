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
