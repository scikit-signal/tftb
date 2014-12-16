import numpy as np


def atoms(N, coord, display):
    pass

def doppler(N, Fs, f0, d, v, t0=None, c=340.0):
    if t0 is None:
        t0 = N/2

    if d <= 0.0:
        raise TypeError("d must be strictly positive.")
    elif Fs < 0.0:
        raise TypeError("Sampling frequency must be positive.")
    elif (t0 < 1) or (t0 > N):
        raise TypeError("T0 must be between 1 and N")
    elif (f0 < 0) or (f0 > Fs/2):
        raise TypeError("F0 must be between 0 and Fs/2")
    elif v < 0:
        raise TypeError("v must be positive")

    tmt0 = (np.arange(N) - t0)/Fs
    dist = np.sqrt(d**2 + (v*tmt0)**2)
    fm = np.exp(1j*2*pi*f0*(tmt0-dist/c))
    if np.abs(f0) < np.spacing(1):
        am = 0
    else:
        am = 1./np.sqrt(dist)
    iflaw = (1 - v**2*tmt0/dist/c)*f0*Fs
    return fm, am, iflaw


def klauder(N, lam=10.0, f0=0.2):
    
    assert N > 0
    assert ((f0 < 0.5) and (f0 > 0))

    f = np.linspace(0,0.5,N/2+1)
    mod = np.exp(-2*np.pi*lam*f)*f**(2*np.pi*lam*f0 - 0.5)
    wave = mod
    wave[0] = 0
    a, b = wave[:N/2], wave[1:N/2+1][::-1]
    wave = np.hstack((a,b))
    wavet = np.fft.ifft(wave)
    wavet = np.fft.fftshift(wavet)
    x = np.real(wavet)/np.linalg.norm(wavet)
    return x


def mexhat(nu=0.05):
    assert (nu <= 0.5) and (nu >= 0)
    N = 1.5
    alpha = np.pi ** 2 * nu ** 2
    n = np.ceil(N/nu)
    t = np.arange(-n,n+1,step=1)
    h = nu*np.sqrt(np.pi)/2*np.exp(-alpha*t**2)*(1-2*alpha*t**2)
    return h


def gdpower(N, k=0, c=1):
    t0 = 0
    lnu = np.round(N/2)
    nu = np.linspace(0,lnu+1,0.5)
    nu = nu[1:lnu]
    am = nu**((k-2)/6)

    if c == 0:
        raise TypeError("c must be non-zero")

    t = np.arange(N)
    tfx = np.zeros((N,),dtype=float)

    if (k<1) and (k!=0):
        d = N**(k*c)
        t0 = N/10
        tfx[:lnu] = np.exp(-j*2*pi*(to*(nu)+d*(nu)**k/k))*am
        x = np.fft.ifft(tfx)
    elif k == 1:
        from analytic_signals import anapulse
        t0 = N
        x = anapulse(N, t0)
    elif k > 1:
        d = N*2**(k-1)*c
        tfx[:lnu] = np.exp(-1j*2*pi*(t0*(nu)+d*(nu)**k/k))*am
        x = np.fft.ifft(tfx)
    else:
        t0 = N/10
        tfx[:lnu] = np.exp(-1j*2*np.pi*(t0*(nu)+d*log(nu)))*am
        x = np.fft.ifft(tfx)

    if k != 1:
        gpsd = t0+np.abs(np.sign(c)-1)/2*(N+1) + c*nu**(k-1)
    else:
        gpd = t0*np.ones((N/2,))

    x = x - x.mean()
    x = x/np.linalg.norm(x)

    return x

if __name__ == "__main__":
    from matplotlib.pyplot import plot, show
    plot(mexhat())
    show()
