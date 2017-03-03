import numpy as np
from numpy import pi
from scipy.signal import hilbert
from tftb.generators import fmconst


def anaask(n_points, n_comp=None, f0=0.25):
    """Generate an amplitude shift (ASK) keying signal.

    :param n_points: number of points.
    :param n_comp: number of points of each component.
    :param f0: normalized frequency.
    :type n_points: int
    :type n_comp: int
    :type f0: float
    :return: Tuple containing the modulated signal and the amplitude modulation.
    :rtype: tuple(numpy.ndarray)
    :Examples:
    >>> x, am = anaask(512, 64, 0.05)
    >>> subplot(211), plot(real(x)) #doctest: +SKIP
    >>> subplot(212), plot(am)      #doctest: +SKIP

    .. plot:: docstring_plots/generators/analytic_signals/anaask.py
    """
    if n_comp is None:
        n_comp = round(n_points / 2.0)
    if (f0 < 0) or (f0 > 0.5):
        raise TypeError("f0 must be between 0 and 0.5")
    m = int(np.ceil(n_points / n_comp))
    jumps = np.random.rand(m)
    am = np.repeat(jumps, n_comp)[:n_points]
    fm, _ = fmconst(n_points, f0, 1)
    y = am * fm
    return y, am


def anabpsk(n_points, n_comp=None, f0=0.25):
    """Binary phase shift keying (BPSK) signal.

    :param n_points: number of points.
    :param n_comp: number of points in each component.
    :param f0: normalized frequency.
    :type n_points: int
    :type n_comp: int
    :type f0: float
    :return: BPSK signal
    :rtype: numpy.ndarray
    :Examples:
    >>> x, am = anabpsk(300, 30, 0.1)
    >>> subplot(211), plot(real(x)) #doctest: +SKIP
    >>> subplot(212), plot(am)      #doctest: +SKIP

    .. plot:: docstring_plots/generators/analytic_signals/anabpsk.py
    """
    if n_comp is None:
        n_comp = round(n_points / 5.0)
    if (f0 < 0) or (f0 > 0.5):
        raise TypeError("f0 must be between 0 and 0.5")
    m = int(np.ceil(n_points / n_comp))
    jumps = 2.0 * np.round(np.random.rand(m)) - 1
    am = np.repeat(jumps, n_comp)[:n_points]
    y = am * fmconst(n_points, f0, 1)[0]
    return y, am


def anafsk(n_points, n_comp=None, Nbf=4):
    """Frequency shift keying (FSK) signal.

    :param n_points: number of points.
    :param n_comp: number of points in each components.
    :param Nbf: number of distinct frequencies.
    :type n_points: int
    :type n_comp: int
    :type Nbf: int
    :return: FSK signal.
    :rtype: numpy.ndarray
    :Examples:
    >>> x, am = anafsk(512, 54.0, 5.0)
    >>> subplot(211), plot(real(x)) #doctest: +SKIP
    >>> subplot(212), plot(am)      #doctest: +SKIP

    .. plot:: docstring_plots/generators/analytic_signals/anafsk.py
    """
    if n_comp is None:
        n_comp = int(round(n_points / 5.0))
    m = np.ceil(n_points / n_comp)
    m = int(np.ceil(n_points / n_comp))
    freqs = 0.25 + 0.25 * (np.floor(Nbf * np.random.rand(m, 1)) / Nbf - (Nbf - 1) / (2 * Nbf))
    iflaw = np.repeat(freqs, n_comp).ravel()
    y = np.exp(1j * 2 * pi * np.cumsum(iflaw))
    return y, iflaw


def anapulse(n_points, ti=None):
    """Analytic projection of unit amplitude impulse signal.

    :param n_points: Number of points.
    :param ti: time position of the impulse.
    :type n_points: int
    :type ti: float
    :return: analytic impulse signal.
    :rtype: numpy.ndarray
    :Examples:
    >>> x = 2.5 * anapulse(512, 301)
    >>> plot(real(x)) #doctest: +SKIP

    .. plot:: docstring_plots/generators/analytic_signals/anapulse.py
    """
    if ti is None:
        ti = np.round(n_points / 2)
    t = np.arange(n_points)
    x = t == ti
    y = hilbert(x.astype(float))
    return y


def anaqpsk(n_points, n_comp=None, f0=0.25):
    """Quaternary Phase Shift Keying (QPSK) signal.

    :param n_points: number of points.
    :param n_comp: number of points in each component.
    :param f0: normalized frequency
    :type n_points: int
    :type n_comp: int
    :type f0: float
    :return: complex phase modulated signal of normalized frequency f0 and
        initial phase.
    :rtype: tuple
    :Examples:
    >>> x, phase = anaqpsk(512, 64.0, 0.05)
    >>> subplot(211), plot(real(x)) #doctest: +SKIP
    >>> subplot(212), plot(phase)   #doctest: +SKIP

    .. plot:: docstring_plots/generators/analytic_signals/anaqpsk.py
    """
    if n_comp is None:
        n_comp = round(n_points / 5.0)
    if (f0 < 0) or (f0 > 0.5):
        raise TypeError("f0 must be between 0 and 0.5")
    m = int(np.ceil(n_points / n_comp))
    jumps = np.floor(4 * np.random.rand(m))
    jumps[jumps == 4] = 3
    pm0 = (np.pi * np.repeat(jumps, n_comp) / 2).ravel()
    tm = np.arange(n_points) - 1
    pm = 2 * np.pi * f0 * tm + pm0
    y = np.exp(1j * pm)
    return y, pm0


def anasing(n_points, t0=None, h=0.0):
    """Lipschitz singularity.
    Refer to the wiki page on `Lipschitz condition`, good test case.

    :param n_points: number of points in time.
    :param t0: time localization of singularity
    :param h: strength of the singularity
    :type n_points: int
    :type t0: float
    :type h: float
    :return: N-point Lipschitz singularity centered around t0
    :rtype: numpy.ndarray
    :Examples:
    >>> x = anasing(128)
    >>> plot(real(x)) #doctest: +SKIP

    .. plot:: docstring_plots/generators/analytic_signals/anasing.py
    """
    if t0 is None:
        t0 = n_points / 2.0
    if h <= 0:
        start, end = 1.0 / n_points, 0.5 - 1.0 / n_points
        N = end / start
        f = np.linspace(start, end, N)
        y = np.zeros((int(n_points / 2.0),), dtype=complex)
        y[1:int(n_points / 2)] = (f ** (-1 - h)) * np.exp(-1j * 2 * pi * f * (t0 - 1))
        x = np.real(np.fft.ifft(y, n_points))
        x = x / x.max()
        x = x - np.sign(x.min()) * np.abs(x.min())
    else:
        t = np.arange(n_points)
        x = np.abs(t - t0) ** h
        x = x.max() - x
    x = hilbert(x)
    return x


def anastep(n_points, ti=None):
    """Analytic projection of unit step signal.

    :param n_points: Number of points.
    :param ti: starting position of unit step.
    :type n_points: int
    :type ti: float
    :return: output signal
    :rtype: numpy.ndarray
    :Examples:
    >>> x = anastep(256, 128)
    >>> plot(real(x)) #doctest: +SKIP

    .. plot:: docstring_plots/generators/analytic_signals/anastep.py
    """
    if ti is None:
        ti = np.round(n_points / 2)
    t = np.arange(n_points)
    x = t > ti
    y = hilbert(x.astype(float))
    return y

if __name__ == '__main__':
    anaask(128)
    anabpsk(128, 128)
    anafsk(128)
    anapulse(128)
    anasing(128)
    anastep(128)
