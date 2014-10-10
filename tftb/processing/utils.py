import numpy as np

def tftb_window(N, name='hamming', param=None, param2=None):
    name = name.lower()
    if name in ('rectangle','rect','rectang'):
        h = np.ones((N,))
    elif name in ('hamming','ham','hamm'):
        h = np.hamming(N)
    elif name in ('hanning','han','hann'):
        h = np.hanning(N)
    elif name == 'kaiser':
        if param is None:
            param = 3*np.pi
        h = np.kaiser(N, beta=param)
    elif name == 'nuttall':
        from scipy.signal import nuttall
        h = nuttall(N)
    elif name == 'blackman':
        h = np.blackman(N)
    elif name == 'harris':
        from scipy.signal import blackmanharris
        h = blackmanharris(N)
    elif name in ('bartlett','bart', #eat my shorts,
                  'triang','triangle'):
        h = np.bartlett(N)
    elif name == 'barthann':
        from scipy.signal import barthann
        h = barthann(N)
    elif name == 'papoulis':
        inds = np.arange(1,N+1)*np.pi/(N+1)
        h = np.sin(inds)
    elif name == 'gauss':
        if param is None:
            param = 0.005
        h = np.exp(np.log(param)*np.linspace(-1.0,1.0,N)**2)
    elif name in ('dolph','dolf','dolph-chebyshev'):
        if param is None:
            param = 100
        from scipy.signal import chebwin
        h = chebwin(N, param)
    return h

