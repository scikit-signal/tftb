import numpy as np


def amgauss(n_points, t0=None, spread=None):
    """Generate a Gaussian amplitude modulator.

    :param n_points: Number of points in the output.
    :param t0: Center of the Gaussian function.
    :param spread: Standard deviation of the Gaussian.
    :type n_points: int
    :type t0: float
    :type spread: float
    :return: Gaussian function centered at time `t0`.
    :rtype: numpy.ndarray
    """
    if t0 is None:
        t0 = np.round(n_points / 2)

    if spread is None:
        spread = 2 * np.sqrt(n_points)

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
    """
    if t0 is None:
        t0 = np.round(n_points / 2)
    if spread is None:
        spread = 2 * np.sqrt(n_points)

    if n_points <= 0:
        raise TypeError
    else:
        tmt0 = np.arange(n_points) - t0
        if kind == "bilateral":
            y = np.exp(-np.sqrt(2 * np.pi) * np.abs(tmt0) / spread)
        else:
            y = np.exp(-np.sqrt(np.pi) * tmt0 / spread) * (tmt0 >= 0.0)
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
    """
    if t0 is None:
        t0 = np.round(n_points / 2)
    if spread is None:
        spread = 2 * np.sqrt(n_points)

    if n_points <= 0:
        raise TypeError
    else:
        tmt0 = np.arange(n_points) - t0
        y = np.abs(tmt0) <= 0.5 * spread * np.sqrt(3.0 / np.pi)
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
    """
    if t0 is None:
        t0 = np.round(n_points / 2)
    if spread is None:
        spread = 2 * np.sqrt(n_points)

    if n_points <= 0:
        raise TypeError
    else:
        tmt0 = np.arange(n_points) - t0
        L = np.sqrt(10.0 / np.pi) * spread / 2.0
        t = np.amin(np.vstack((L + tmt0, L - tmt0)).T, axis=1)
        t = np.hstack((t.reshape((len(t), 1)), np.zeros((len(t), 1))))
        y = np.amax(t, axis=1) / L

        return y
