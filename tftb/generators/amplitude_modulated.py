import numpy as np
from math import sqrt


def amgauss(n_points, t0=None, spread=None):
    """Generate a Gaussian amplitude modulated signal.

    :param n_points: Number of points in the output.
    :param t0: Center of the Gaussian function. (default: t0 / 2)
    :param spread: Standard deviation of the Gaussian. (default 2 *
        sqrt(n_points)
    :type n_points: int
    :type t0: float
    :type spread: float
    :return: Gaussian function centered at time ``t0``.
    :rtype: numpy.ndarray
    :Example:
    >>> x = amgauss(160)
    >>> plot(x) #doctest: +SKIP

    .. plot:: docstring_plots/generators/amplitude_modulated/amgauss1.py
    >>> x = amgauss(160, 90)
    >>> plot(x) #doctest: +SKIP

    .. plot:: docstring_plots/generators/amplitude_modulated/amgauss2.py
    """
    if t0 is None:
        t0 = round(n_points / 2.0)

    if spread is None:
        spread = 2 * sqrt(n_points)

    if n_points <= 0:
        raise TypeError("n_points should be >= 0")
    else:
        tmt0 = np.arange(1, n_points + 1, dtype=float) - t0
        y = np.exp(-((tmt0 / spread) ** 2) * np.pi)
        return y


def amexpos(n_points, t0=None, spread=None, kind="bilateral"):
    """Exponential amplitude modulation.

    `amexpos` generates an exponential amplitude modulation starting at time
    `t0` and spread proportioanl to `spread`.

    :param n_points: Number of points.
    :param kind: "bilateral" (default) or "unilateral"
    :param t0: Time center.
    :param spread: Standard deviation.
    :type n_points: int
    :type kind: str
    :type t0: float
    :type spread: float
    :return: exponential function
    :rtype: numpy.ndarray
    :Examples:
    >>> x = amexpos(160)
    >>> plot(x) #doctest: +SKIP

    .. plot:: docstring_plots/generators/amplitude_modulated/amexpos_bilateral.py
    >>> x = amexpos(160, kind='unilateral')
    >>> plot(x) #doctest: +SKIP

    .. plot:: docstring_plots/generators/amplitude_modulated/amexpos_unilateral.py
    """
    if t0 is None:
        t0 = round(n_points / 2.0)
    if spread is None:
        spread = 2 * sqrt(n_points)

    if n_points <= 0:
        raise TypeError
    else:
        tmt0 = np.arange(n_points) - t0
        if kind == "bilateral":
            y = np.exp(-sqrt(2 * np.pi) * np.abs(tmt0) / spread)
        else:
            y = np.exp(-sqrt(np.pi) * tmt0 / spread) * (tmt0 >= 0.0)
        return y


def amrect(n_points, t0=None, spread=None):
    """Generate a rectangular amplitude modulation.

    :param n_points: Number of points in the function.
    :param t0: Time center
    :param spread: standard deviation of the function.
    :type n_points: int
    :type t0: float
    :type spread: float
    :return: A rectangular amplitude modulator.
    :rtype: numpy.ndarray.
    :Examples:
    >>> x = amrect(160, 90, 40.0)
    >>> plot(x) #doctest: +SKIP

    .. plot:: docstring_plots/generators/amplitude_modulated/amrect1.py
    """
    if t0 is None:
        t0 = round(n_points / 2.0)
    if spread is None:
        spread = 2 * sqrt(n_points)

    if n_points <= 0:
        raise TypeError
    else:
        tmt0 = np.arange(n_points) - t0
        y = np.abs(tmt0) <= 0.5 * spread * sqrt(3.0 / np.pi)
        return y


def amtriang(n_points, t0=None, spread=None):
    """Generate a triangular amplitude modulation.

    :param n_points: Number of points in the function.
    :param t0: Time center
    :param spread: standard deviation of the function.
    :type n_points: int
    :type t0: float
    :type spread: float
    :return: A triangular amplitude modulator.
    :rtype: numpy.ndarray.
    :Examples:
    >>> x = amtriang(160)
    >>> plot(x) #doctest: +SKIP

    .. plot:: docstring_plots/generators/amplitude_modulated/amtriang1.py
    """
    if t0 is None:
        t0 = round(n_points / 2.0)
    if spread is None:
        spread = 2 * sqrt(n_points)

    if n_points <= 0:
        raise TypeError
    else:
        tmt0 = np.arange(n_points) - t0
        L = sqrt(10.0 / np.pi) * spread / 2.0
        t = np.amin(np.c_[L + tmt0, L - tmt0], axis=1)
        t = np.c_[t, np.zeros(t.shape)]
        y = np.amax(t, axis=1) / L

        return y


if __name__ == '__main__':
    amgauss(128)
    amexpos(128)
    amrect(128)
    amtriang(128)
