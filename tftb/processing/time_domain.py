import numpy as np


def loctime(sig):
    """
    Compute the time localization characteristics.

    Parameters
    ----------
    sig : array-like
        Input signal.

    Returns
    -------
    tm : float
        Averaged time center.

    T : float
        Time spreading.

    """
    if sig.ndim > 2:
        if 1 not in sig.shape:
            raise TypeError
        else:
            sig = sig.ravel()
    sig2 = np.abs(sig**2)
    sig2 = sig2/sig2.mean()
    t = np.arange(len(sig))
    tm = np.mean(t*sig2)
    T = 2*np.sqrt(np.pi*np.mean(((t-tm)**2)*sig2))
    return tm, T


def group_delay(x, fnorm=None):
    """
    Compute the group delay of a signal at normalized frequencies.

    Paramters
    ---------
    x : array-like
        Time-domain signal

    fnorm : array-like, optional
        Normalized frequency points at which to calculate the group delay.
        By default, this is an array of len(x) points spaced linearly between
        -0.5 and 0.5.

    Returns
    -------
    gd : 1-D ndarray
        The computed group delay of the signal `x`.

    """
    if x.ndim != 1:
        if 1 not in x.shape:
            raise TypeError
        else:
            sig = x.ravel()

    if fnorm is None:
        numerator = np.fft.fft(x * np.arange(1,x.shape[0]+1))
        denominator = np.fft.fft(x)
        window = np.real(numerator / denominator) >= 1
        ratio = np.real(numerator / denominator) * window.astype(int)
        ratio = ratio * (np.real(numerator / denominator) <= (len(x) +
                                                              3)).astype(int)
        gd = np.fft.fftshift(ratio)
    else:
        exponent = np.exp(-1j * 2.0 * np.pi * fnorm.reshape(len(fnorm), 1) * \
                                                              np.arange(len(x)))
        numerator = np.dot(exponent, (x * np.arange(1, x.shape[0]+1)))
        denominator = np.dot(exponent, x)
        window = np.real(numerator / denominator) >= 1
        ratio = np.real(numerator / denominator) * window.astype(int)
        gd = ratio * (np.real(numerator / denominator) <= (len(x) +
                                                              3)).astype(int)
    return gd
