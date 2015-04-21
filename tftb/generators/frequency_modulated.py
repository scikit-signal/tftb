import numpy as np
from numpy import pi
import warnings


def fmconst(n_points, fnorm=0.25, t0=None):
    """
    Generate a signal with constant frequency modulation.

    :param n_points: number of points
    :param fnorm: normalized frequency
    :param t0: time center
    :type n_points: int
    :type fnorm: float
    :type t0: float
    :return: frequency modulation signal with frequency fnorm
    :rtype: numpy.ndarray
    """
    if t0 is None:
        t0 = np.round(n_points / 2)

    if n_points <= 0:
        raise TypeError
    elif abs(fnorm) > 0.5:
        raise TypeError
    else:
        tmt0 = np.arange(n_points) - t0
        y = np.exp(1j * 2.0 * np.pi * fnorm * tmt0)
        y = y / y[t0]
        iflaw = fnorm * np.ones((n_points, 1))
        return y, iflaw


def fmhyp(n_points, p1, p2):
    """Signal with hyperbolic frequency modulation.

    :param n_points: number of points.
    :param p1, p2: coefficients of the hyperbolic function.
    :type n_points: int
    :type p1: float
    :type p2: float
    :return: vector containing the modulated signal samples.
    :rtype: numpy.ndarray
    """
    if (len(p1) != 2) or (len(p2) != 2):
        raise TypeError
    elif n_points <= 0:
        raise TypeError

    if (p1[0] > n_points) or (p1[0] < 1):
        raise TypeError
    elif (p2[0] > n_points) or (p2[0] < 1):
        raise TypeError
    elif (p1[1] < 0) or (p2[1] < 0):
        raise TypeError

    c = (p2[1] - p1[1]) / (1.0 / p2[0] - 1 / p1[0])
    f0 = p1[1] - c / p1[0]

    t = np.arange(n_points)
    phi = 2 * np.pi * (f0 * t + c * np.log(np.abs(t)))
    iflaw = (f0 + c * np.abs(t)) ** -1

    a, b = iflaw < 0, iflaw > 0.5
    aliasing = np.logical_or(a, b)
    if np.any(aliasing):
        msg = "Signal may be undersampled or may have negative frequencies."
        warnings.warn(msg, UserWarning)

    x = np.exp(1j * phi)

    return x, iflaw


def fmlin(n_points, fnormi=0.0, fnormf=0.5, t0=None):
    """
    Generate a signal with linear frequency modulation.

    :param n_points: number of points
    :param fnormi: initial normalized frequency
    :param fnormf: final normalized frequency
    :param t0: time center
    :type n_points: int
    :type fnormi: float
    :type fnormf: float
    :type t0: float
    :return: The modulated signal, and the instantaneous amplitude law.
    :rtype: numpy.ndarray
    """
    if t0 is None:
        t0 = np.round(n_points / 2)

    if n_points <= 0:
        raise TypeError

    elif (np.abs(fnormi) > 0.5) or (np.abs(fnormf) > 0.5):
        raise TypeError

    else:
        y = np.arange(1, n_points + 1)
        y = fnormi * (y - t0) + ((fnormf - fnormi) / (2.0 * (n_points - 1))) * \
            ((y - 1) ** 2 - (t0 - 1) ** 2)
        y = np.exp(1j * 2.0 * np.pi * y)
        y = y / y[t0 - 1]
        iflaw = np.linspace(fnormi, fnormf, n_points)
        return y, iflaw


def fmodany(iflaw, t0=1):
    """Arbitrary frequency modulation.

    :param iflaw: Vector of instantaneous frequency law samples.
    :param t0: time center
    :type iflaw: numpy.ndarray
    :type t0: float
    :return: output signal
    :rtype:
    """
    if len(iflaw.shape) > 1:
        if iflaw.shape[1] != 1:
            raise TypeError("iflaw should be a column vector.")

    elif np.amax(np.abs(iflaw)) > 0.5:
        raise TypeError("Elements of iflaw should be within -0.5 and 0.5")

    if (t0 > iflaw.shape[0]) or (t0 == 0):
        raise TypeError("T0 should be between 1 and len(iflaw)")

    y = np.exp(1j * 2.0 * pi * np.cumsum(iflaw))
    y = y * np.conjugate(y[t0])
    return y


def fmpar(n_points, p1, p2=None, p3=None):
    """Parabolic frequency modulated signal.

    :param n_points: number of points
    :param p1, p2, p3: coefficients of the parabolic function.
    :type n_points: int
    :type p1: float
    :type p2: float
    :type p3: float
    :return: Signal with parabolic frequency modulation law.
    :rtype:
    """
    a0, a1, a2 = p1
    t = np.arange(n_points)
    phi = 2 * np.pi * (a0 * t + (a1 / 2 * t) ** 2 + (a2 / 3 * t) ** 3)
    iflaw = a0 + a1 * t + (a2 * t) ** 2
    a, b = iflaw < 0, iflaw > 0.5
    aliasing = np.logical_or(a, b)
    if np.any(aliasing):
        msg = "Signal may be undersampled or may have negative frequencies."
        warnings.warn(msg, UserWarning)
    x = np.exp(1j * phi)
    return x


def fmpower(n_points, k, p1, p2=None):
    """Generate signal with power law frequency modulation.

    :param n_points: number of points.
    :param k: degree of power law.
    :param p1, p2: coefficients of the power law.
    :type n_points: int
    :type k: int
    :type p1: float
    :type p2: float
    :return: vector of modulated signal samples.
    :rtype: numpy.ndarray
    """
    f0, c = p1
    t = np.arange(n_points)
    phi = 2 * np.pi * (f0 * t + c / (1 - k) * np.abs(t) ** (1 - k))
    iflaw = (f0 + c * np.abs(t) ** (-k))
    a, b = iflaw < 0, iflaw > 0.5
    aliasing = np.logical_or(a, b)
    if np.any(aliasing):
        msg = "Signal may be undersampled or may have negative frequencies."
        warnings.warn(msg, UserWarning)

    x = np.exp(1j * phi)
    return x


def fmsin(n_points, fnormin=0.05, fnormax=0.45, period=None, t0=None,
          fnorm0=None, pm1=1):
    """Sinusodial frequency modulation

    :param n_points: number of points
    :param fnormin: smallest normalized frequency
    :param fnormax: highest normalized frequency
    :param period: period of sinusoidal fm
    :param t0: time reference
    :param fnorm0: normalized frequency at time t0
    :param pm1: frequency direction at t0
    :type n_points: int
    :type fnormin: float
    :type fnormax: float
    :type period: int
    :type t0: float
    :type fnorm0: float
    :type pm1: int
    :return: output signal
    :rtype: numpy.ndarray
    """
    fnormid = 0.5 * (fnormax + fnormin)
    delta = 0.5 * (fnormax - fnormin)
    phi = -pm1 * np.arccos((fnorm0 - fnormid) / delta)
    t = np.arange(n_points) - t0
    phase = 2 * pi * fnormid * t + delta * period * \
            (np.sin(2 * pi * t / period + phi)) - np.sin(phi)
    y = np.exp(1j * phase)
    iflaw = fnormid + delta * np.cos(2 * pi * t / period + phi)
    return y, iflaw
