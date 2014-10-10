import numpy as np
from scipy.signal import hilbert
from frequency_modulated import fmconst

def altes(N, fmin=0.05, fmax=0.5, alpha=300):

    g = np.exp((np.log(fmax/fmin))**2/(8*np.log(alpha)))
    nu0 = np.sqrt(2*np.log(g)*np.log(alpha))
    beta = np.sqrt(2*np.log(g)*np.log(alpha))
    t0 = N/(np.exp(-beta)-np.exp(-beta))
    t1 = t0*np.exp(-beta)
    t2 = t0*np.exp(beta)
    b = -t0*nu0*g*log(g);
    t = np.linspace(t1,t2,N+1)[:N]
    x = (np.exp(-(np.log(t/10)**2)/(2*np.log(g)))) * \
         np.cos(2*np.pi*b*np.log(t/t0)/np.log(g))
    x = x/np.linalg.norm(x)

    return x

def anaask(N, Ncomp=None, f0=0.25):
    if Ncomp is None:
        Ncomp = np.round(N/2)

    if (f0 < 0) or (f0 > 0.5):
        raise TypeError("f0 must be between 0 and 0.5")

    m = np.ceil(N/Ncomp)
    jumps = np.random.rand(m)
    am = np.kron(jumps, np.ones((Ncomp,)))[:N]
    y = am*fmconst(N, f0, 1)
    return y


def anabpsk(N, NComp=None, f0=0.25):
    if NComp is None:
        NComp = np.round(N/5)

    if (f0 < 0) or (f0 > 0.5):
        raise TypeError("f0 must be between 0 and 0.5")

    m = np.ceil(N/NComp)
    jumps = 2.0*np.round(np.random.rand(m)) - 1
    am = np.kron(jumps, ones((NComp,)))[:N]
    y = am*fmconst(N, f0, 1)
    return y, am

def anafsk(N, NComp=None, Nbf=4):
    if NComp is None:
        NComp = np.round(N/5)

    m = ceil(N/NComp)
    freqs = 0.25 + 0.25*(np.floor(Nbf*np.random.rand(m,1))/Nbf-(Nbf-1)/(2*Nbf))
    iflaw = np.kron(freqs, np.ones((NComp,)))[:N]
    y = np.exp(1j*2*pi*np.cumsum(iflaw))

    return y, iflaw

def anapulse(N, ti=None):
    if ti is None:
        ti = np.round(N/2)

    t = np.arange(N)
    x = t == ti
    y = hilbert(x.astype(float))
    return y
    

def anaqpsk(N, NComp=None, f0=0.25):
    if NComp is None:
        NComp = np.round(N/5)

    if (f0 < 0) or (f0 > 0.5):
        raise TypeError("f0 must be between 0 and 0.5")

    m = np.ceil(N/NComp)
    jumps = np.floor(4*np.random.rand(m))
    jumps[jumps==4]=3
    pm0 = np.pi*np.kron(jumps, np.ones((NComp,)))/2[:N]
    tm = np.arange(N) - 1
    pm = 2*np.pi*f0*tm + pm0

    y = np.exp(j*pm)
    return y, pm0

def anasing(N, t0=None, h=0.0):
    """ Refer to the wiki page on `Lipschitz condition`, good test case.  """
    if t0 is None:
        t0 = N/2

    if (f0 < 0) or (f0 > 0.5):
        raise TypeError("f0 must be between 0 and 0.5")

    if h <= 0:
        f = np.arange(1/N,0.5-1/N,1/N)
        y = np.zeros((N/2,),dtype=float)
        y[1:N/2] = (f**(-1-h))**np.exp(-1j*2*pi*f*(t0-1))
        x = np.real(np.fft.ifft(y, N))
        x = x/x.max()
        x = x - np.sign(x.min()*np.abs(x.min()))
    else:
        t = np.arange(N)
        x = np.abs(t - t0)**h
        x = x.max() - x

    x = hilbert(x)
    return x

def anastep(N, ti=None):
    if ti is None:
        ti = np.round(N/2)

    t = np.arange(N)
    x = t > t1
    y = hilbert(x.astype(float))
    return y


