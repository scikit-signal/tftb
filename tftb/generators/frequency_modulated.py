import numpy as np

def fmconst(N, fnorm=0.25, t0=None):
    """
    Generate a signal with constant frequency modulation.

    Parameters
    ----------
    N : int
        Length of the desired signal

    fnorm: float, optional
        Normalized frequency of the output. Default is 0.25

    t0 : int, optional
        Time center. Default is N/2

    Returns
    -------
    y : 1-D ndarrays, shape (N,)
        The output signal.

    iflaw : 1-D ndarray, shape(N,)
        The instantaneous frequency law.

    """
    if t0 is None:
        t0 = np.round(N/2)

    if N <= 0:
        raise TypeError
    elif abs(fnorm) > 0.5:
        raise TypeError
    else:
        tmt0 = np.arange(N) - t0
        y = np.exp(1j*2.0*np.pi*fnorm*tmt0)
        y = y/y[t0]
        iflaw = fnorm*np.ones((N,1))
        return y, iflaw


def fmhyp(N, p1, p2):
    if (len(p1)!=2) or (len(p2)!=2):
        raise TypeError
    elif N <= 0:
        raise TypeError

    if (p1[0] > N) or (p1[0] < 1):
        raise TypeError
    elif (p2[0]>N) or (p2[0]<1):
        raise TypeError
    elif (p1[1]<0) or (p2[1]<0):
        raise TypeError

    c = (p2[1] - p1[1])/(1.0/p2[0] - 1/p1[0])
    f0 = p1[1] - c/p1[0]

    t = np.arange(N)
    phi = 2 * np.pi * (f0*t + c*np.log(np.abs(t)))
    iflaw = (f0 + c*np.abs(t))**(-1)

    a, b = iflaw < 0, iflaw > 0.5
    aliasing = np.logical_or(a,b)
    if np.any(aliasing):
        print "WARNING: Singal may be undersampled or may have negative frequencies."

    x = np.exp(1j*phi);

    return x, iflaw


def fmlin(N, fnormi=0.0, fnormf=0.5, t0=None):
    """
    Generate a signal with linear frequency modulation.

    Parameters
    ----------
    N : int
        The desired length of the output.

    fnormi : float, optional
        Initial normalized frequency. (default 0.0)

    fnormf : float, optional
        Final normalized frequency. (default 0.5)

    t0 : int
        Time reference for the phase. (default N/2)

    Returns
    -------
    y : 1-D ndarray, shape (N,)
        The modulated signal.

    iflaw : 1-D ndarray, shape (N,)
        The instantaneous amplitude law.

    """

    if t0 is None:
        t0 = np.round(N/2)

    if N <= 0:
        raise TypeError

    elif (np.abs(fnormi) > 0.5) or (np.abs(fnormf) > 0.5):
        raise TypeError

    else:
        y = np.arange(1,N+1)
        y = fnormi * (y - t0) + ((fnormf - fnormi)/(2.0*(N - 1))) * \
            ((y - 1)**2 - (t0 - 1)**2)
        y = np.exp(1j*2.0*np.pi*y)
        y = y/y[t0-1]
        iflaw = np.linspace(fnormi, fnormf, N)

        return y, iflaw


def fmodany(iflaw, t0=1):
    if len(iflaw.shape) > 1:
        if iflaw.shape[1] != 1:
            raise TypeError("iflaw should be a column vector.")

    elif np.amax(np.abs(iflaw)) > 0.5:
        raise TypeError("Elements of iflaw should be within -0.5 and 0.5")

    if (t0 > iflaw.shape[0]) or (t0 == 0):
        raise TypeError("T0 should be between 1 and len(iflaw)")

    y = exp(j*2.0*pi*np.cumsum(iflaw))
    y = y * np.conjugate(y[t0])
    return y

def fmpar(N, p1, p2=None, p3=None):
    a0, a1, a2 = p1
    t = np.arange(N)
    phi = 2*np.pi*(a0*t + (a1/2*t)**2 + (a2/3*t)**3)
    iflaw = a0 + a1*t + (a2*t)**2
    a, b = iflaw < 0, iflaw > 0.5
    aliasing = np.logical_or(a,b)
    if np.any(aliasing):
        print "WARNING: Singal may be undersampled or may have negative frequencies."

    x = np.exp(1j*phi)
    return x

def fmpower(N, k, p1, p2=None):
    f0, c = p1
    t = np.arange(N)
    phi = 2*np.pi*(f0*t + c/(1-k)*np.abs(t)**(1-k))
    iflaw = (f0 + c*np.abs(t)**(-k))

    
    a, b = iflaw < 0, iflaw > 0.5
    aliasing = np.logical_or(a,b)
    if np.any(aliasing):
        print "WARNING: Singal may be undersampled or may have negative frequencies."

    x = np.exp(1j*phi)
    return x

def fmsin(N, fnormin=0.05, fnormax=0.45, period=None, t0=None, fnorm0=None,
          pm1=1):
    fnormid = 0.5*(fnormax + fnormin)
    delta   = 0.5*(fnormax - fnormin)
    phi = -pm1*np.arccos((fnorm0 - fnormid)/delta);
    t = arange(N) - t0
    phase = 2*pi*fnormid*t+delta*period*(np.sin(2*pi*t/period+phi)) - \
            np.sin(phi)
    y = np.exp(1j*phase)
    iflaw = fnormid + delta*cos(2*pi*t/period+phi)
    return y, iflaw






if __name__ == "__main__":
    from matplotlib.pyplot import plot,subplot, show
    y, iflaw =fmlin(128,0.0,0.5)
    subplot(211),plot(np.real(y))
    subplot(212),plot(np.imag(y))
    show()
