import numpy as np
from scipy import angle


def locfreq(signal):
    """
    Compute the frequency localization characteristics.

    :param sig: input signal
    :type sig: numpy.ndarray
    :return: average normalized frequency center, frequency spreading
    :rtype: tuple
    :Example:
    >>> from tftb.generators import amgauss
    >>> z = amgauss(160, 80, 50)
    >>> fm, B = locfreq(z)
    >>> print("%.4g" % fm)
    -9.183e-14
    >>> print("%.4g" % B)
    0.02
    """
    if signal.ndim > 1:
        if 1 not in signal.shape:
            raise TypeError
        else:
            signal = signal.ravel()
    no2r = np.round(signal.shape[0] / 2.0)
    no2f = np.floor(signal.shape[0] / 2.0)
    sig = np.fft.fft(signal)
    sig = np.abs(sig) ** 2
    sig = sig / sig.mean()
    freqs = np.hstack((np.arange(no2f), np.arange(-no2r, 0))) / signal.shape[0]
    fm = np.mean(freqs * sig)
    bw = 2 * np.sqrt(np.pi * np.mean(((freqs - fm) ** 2) * sig))
    return fm, bw


def inst_freq(x, t=None, L=1):
    """
    Compute the instantaneous frequency of an analytic signal at specific
    time instants using the trapezoidal integration rule.

    :param x: The input analytic signal
    :param t: The time instants at which to calculate the instantaneous frequencies.
    :param L: Non default values are currently not supported.
        If L is 1, the normalized instantaneous frequency is computed. If L > 1,
        the maximum likelihood estimate of the instantaneous frequency of the
        deterministic part of the signal.
    :type x: numpy.ndarray
    :type t: numpy.ndarray
    :type L: int
    :return: instantaneous frequencies of the input signal.
    :rtype: numpy.ndarray
    :Example:
    >>> from tftb.generators import fmsin
    >>> x = fmsin(70, 0.05, 0.35, 25)[0]
    >>> instf, timestamps = inst_freq(x)
    >>> plot(timestamps, instf) #doctest: +SKIP

    .. plot:: docstring_plots/processing/freq_domain/inst_freq.py
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
        t = np.arange(2, len(x))

    fnorm = 0.5 * (angle(-x[t] * np.conj(x[t - 2])) + np.pi) / (2 * np.pi)
    return fnorm, t


def group_delay(x, fnorm=None):
    """
    Compute the group delay of a signal at normalized frequencies.

    :param x: time domain signal
    :param fnorm: normalized frequency at which to calculate the group delay.
    :type x: numpy.ndarray
    :type fnorm: float
    :return: group delay
    :rtype: numpy.ndarray
    :Example:
    >>> import numpy as np
    >>> from tftb.generators import amgauss, fmlin
    >>> x = amgauss(128, 64.0, 30) * fmlin(128, 0.1, 0.4)[0]
    >>> fnorm = np.arange(0.1, 0.38, step=0.04)
    >>> gd = group_delay(x, fnorm)
    >>> plot(gd, fnorm) #doctest: +SKIP

    .. plot:: docstring_plots/processing/freq_domain/group_delay.py
    """
    if x.ndim != 1:
        if 1 not in x.shape:
            raise TypeError
        else:
            x = x.ravel()

    if fnorm is None:
        numerator = np.fft.fft(x * np.arange(1, x.shape[0] + 1))
        denominator = np.fft.fft(x)
        window = np.real(numerator / denominator) >= 1
        ratio = np.real(numerator / denominator) * window.astype(int)
        ratio = ratio * (np.real(numerator / denominator) <= (len(x) +
                                                              3)).astype(int)
        gd = np.fft.fftshift(ratio)
    else:
        exponent = np.exp(-1j * 2.0 * np.pi * fnorm.reshape(len(fnorm), 1) * np.arange(len(x)))
        numerator = np.dot(exponent, (x * np.arange(1, x.shape[0] + 1)))
        denominator = np.dot(exponent, x)
        window = np.real(numerator / denominator) >= 1
        ratio = np.real(numerator / denominator) * window.astype(int)
        gd = ratio * (np.real(numerator / denominator) <= (len(x) + 3)).astype(int)
    return gd
