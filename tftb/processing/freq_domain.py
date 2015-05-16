import numpy as np
from scipy import angle


def locfreq(signal):
    """
    Compute the frequency localization characteristics.

    :param sig: input signal
    :type sig: numpy.ndarray
    :return: average normalized frequency center, frequency spreading
    :rtype: tuple
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
    return fnorm
