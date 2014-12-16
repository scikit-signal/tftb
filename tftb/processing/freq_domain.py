import numpy as np
from scipy import angle


def locfreq(sig):
    """
    Compute the frequency localization characteristics.

    Parameters
    ----------
    sig : array-like
        Input signal.

    Returns
    -------
    fm : float
        Averaged normalized frequency center.

    B : float
        Frequency spreading.

    """
    if sig.ndim > 1:
        if 1 not in sig.shape:
            raise TypeError
        else:
            sig = sig.ravel()
    N = float(len(sig))
    No2r = np.round(N/2)
    No2f = np.floor(N/2)
    Sig = np.fft.fft(sig)
    Sig2 = np.abs(Sig)**2
    Sig2 = Sig2/Sig2.mean()
    freqs = np.r_[np.arange(No2f), np.arange(-No2r, 0)]/N
    fm = np.mean(freqs*Sig2)
    B = 2*np.sqrt(np.pi*np.mean(((freqs-fm)**2)*Sig2))
    return fm, B


def inst_freq(x, t=None, L=1):
    """
    Compute the instantaneous frequency of an analytic signal at specific
    time instants using the trapezoidal integration rule.

    Parameters
    ----------
    x : array-like
        The input analytic signal

    t : array-like, optional
        The time instants at which to calculate the instantaneous frequencies.
        Default: np.arange(2,len(x))

    L : integer, optional
        Non default values are currently not supported.
        If L is 1, the normalized instantaneous frequency is computed. If L > 1,
        the maximum likelihood estimate of the instantaneous frequency of the
        deterministic part of the signal.

    Returns
    -------
    fnorm : 1-D ndarray
        The instantaneous frequencies of the input signal.

    """
    if x.ndim != 1:
        if 1 not in x.shape:
            raise TypeError("Input should be a one dimensional array.")
        else:
            x = x.ravel()
    if t is not None:
        if t.ndim != 1:
            if 1 not in t.shape:
                raise TypeError("Time instants should be a one dimensional "
                                "array.")
            else:
                t = t.ravel()
    else:
        t = np.arange(2,len(x))

    fnorm = 0.5 * (angle(-x[t] * np.conj(x[t-2])) + np.pi) / (2 * np.pi)
    return fnorm