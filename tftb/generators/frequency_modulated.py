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
    :Examples:
    >>> from tftb.generators import amgauss
    >>> z = amgauss(128, 50, 30) * fmconst(128, 0.05, 50)[0]
    >>> plot(real(z)) #doctest: +SKIP

    .. plot:: docstring_plots/generators/frequency_modulated/fmconst.py
    """
    if t0 is None:
        t0 = round(n_points / 2.0)

    if n_points <= 0:
        raise TypeError
    elif abs(fnorm) > 0.5:
        raise TypeError
    else:
        tmt0 = np.arange(n_points) - t0
        y = np.exp(1j * 2.0 * np.pi * fnorm * tmt0)
        y = y / np.real(y[int(t0)])
        iflaw = fnorm * np.ones((n_points,))
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
    :Examples:
    >>> signal, iflaw = fmhyp(128, (1, 0.5), (32, 0.1))
    >>> subplot(211), plot(real(signal)) #doctest: +SKIP
    >>> subplot(212), plot(iflaw)        #doctest: +SKIP

    .. plot:: docstring_plots/generators/frequency_modulated/fmhyp.py
    """
    c = (p2[1] - p1[1]) / (1.0 / p2[0] - 1 / p1[0])
    f0 = p1[1] - c / p1[0]

    t = np.arange(1, n_points + 1)
    phi = 2 * np.pi * (f0 * t + c * np.log(np.abs(t)))
    iflaw = (f0 + c / np.abs(t))

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
    :rtype: tuple(array-like)
    :Examples:
    >>> from tftb.generators import amgauss
    >>> z = amgauss(128, 50, 40) * fmlin(128, 0.05, 0.3, 50)[0]
    >>> plot(real(z)) #doctest: +SKIP

    .. plot:: docstring_plots/generators/frequency_modulated/fmlin.py
    """
    if t0 is None:
        t0 = round(n_points / 2.0)

    if n_points <= 0:
        raise TypeError

    elif (abs(fnormi) > 0.5) or (abs(fnormf) > 0.5):
        raise TypeError

    else:
        y = np.arange(1, n_points + 1)
        y = fnormi * (y - t0) + ((fnormf - fnormi) / (2.0 * (n_points - 1))) * \
            ((y - 1) ** 2 - (t0 - 1) ** 2)
        y = np.exp(1j * 2.0 * np.pi * y)
        y = y / y[int(t0) - 1]
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
    :Examples:
    >>> from tftb.generators import fmlin
    >>> import numpy as np
    >>> y1, ifl1 = fmlin(100)  # A linear instantaneous frequency law.
    >>> y2, ifl2 = fmsin(100)  # A sinusoidal instantaneous frequency law.
    >>> iflaw = np.append(ifl1, ifl2)  # combination of the two
    >>> sig = fmodany(iflaw)
    >>> subplot(211), plot(real(sig)) #doctest: +SKIP
    >>> subplot(212), plot(iflaw)     #doctest: +SKIP

    .. plot:: docstring_plots/generators/frequency_modulated/fmodany.py
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


def fmpar(n_points, coefficients):
    """Parabolic frequency modulated signal.

    :param n_points: number of points
    :param coefficients: coefficients of the parabolic function.
    :type n_points: int
    :type coefficients: tuple
    :return: Signal with parabolic frequency modulation law.
    :rtype: tuple
    :Examples:
    >>> x, iflaw = fmpar(128, (0.4, -0.0112, 8.6806e-05))
    >>> subplot(211), plot(real(x)) #doctest: +SKIP
    >>> subplot(212), plot(iflaw)   #doctest: +SKIP

    .. plot:: docstring_plots/generators/frequency_modulated/fmpar.py
    """
    a0, a1, a2 = coefficients
    t = np.arange(n_points)
    phi = 2 * np.pi * (a0 * t + (a1 / 2 * t) ** 2 + (a2 / 3 * t) ** 3)
    iflaw = a0 + a1 * t + a2 * (t ** 2)
    a, b = iflaw < 0, iflaw > 0.5
    aliasing = np.logical_or(a, b)
    if np.any(aliasing):
        msg = "Signal may be undersampled or may have negative frequencies."
        warnings.warn(msg, UserWarning)
    x = np.exp(1j * phi)
    return x, iflaw


def fmpower(n_points, k, coefficients):
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
    :Examples:
    >>> x, iflaw = fmpower(128, 0.5, (1, 0.5, 100, 0.1))
    >>> subplot(211), plot(real(x)) #doctest: +SKIP
    >>> subplot(212), plot(iflaw)   #doctest: +SKIP

    .. plot:: docstring_plots/generators/frequency_modulated/fmpower.py
    """
    if len(coefficients) == 2:
        f0, c = coefficients
    elif len(coefficients) > 2:
        p1 = coefficients[:2]
        p2 = coefficients[2:]
        c = (p2[1] - p1[1]) / (1.0 / (p2[0] ** k) - 1.0 / (p1[0] ** k))
        f0 = p1[1] - c / (p1[0] ** k)
    t = np.arange(n_points)
    phi = 2 * np.pi * (f0 * t + c / (1 - k) * np.abs(t) ** (1 - k))
    iflaw = (f0 + c * np.abs(t) ** (-k))
    a, b = iflaw < 0, iflaw > 0.5
    aliasing = np.logical_or(a, b)
    if np.any(aliasing):
        msg = "Signal may be undersampled or may have negative frequencies."
        warnings.warn(msg, UserWarning)

    x = np.exp(1j * phi)
    return x, iflaw


def fmsin(n_points, fnormin=0.05, fnormax=0.45, period=None, t0=None,
          fnorm0=None, pm1=1):
    """Sinusodial frequency modulation.

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
    :Examples:
    >>> z = fmsin(140, period=100, t0=20, fnorm0=0.3, pm1=-1)
    >>> plot(real(z)) #doctest: +SKIP

    .. plot:: docstring_plots/generators/frequency_modulated/fmsin.py
    """
    if period is None:
        period = n_points
    if t0 is None:
        t0 = n_points / 2.0
    if fnorm0 is None:
        fnorm0 = 0.5 * (fnormin + fnormax)

    fnormid = 0.5 * (fnormax + fnormin)
    delta = 0.5 * (fnormax - fnormin)
    phi = -pm1 * np.arccos((fnorm0 - fnormid) / delta)
    t = np.arange(n_points) - t0
    phase = 2 * pi * fnormid * t + delta * period * (np.sin(2 * pi * t / period + phi)) -\
        np.sin(phi)
    y = np.exp(1j * phase)
    iflaw = fnormid + delta * np.cos(2 * pi * t / period + phi)
    return y, iflaw


if __name__ == '__main__':
    fmconst(128)
    fmhyp(128, (1, 0.5), (32, 0.1))
    fmlin(128)
    fmodany(np.random.rand(128) - 0.5)
    fmpar(128, (0.4, -0.0112, 8.6806e-05))
    fmpower(128, 0.5, (1, 0.5, 100, 0.1))
    fmsin(128)
