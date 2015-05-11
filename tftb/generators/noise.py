import numpy as np
from utils import nextpow2
from scipy.signal import hilbert


def noisecu(n_points):
    """Compute analytic complex uniform white noise.

    :param n_points: Length of the noise signal.
    :type n_points: int
    :return: analytic complex uniform white noise signal of length N
    :rtype: numpy.ndarray
    """
    if n_points <= 2:
        noise = (np.random.rand(n_points, 1) - 0.5 + 1j * (np.random.rand(n_points, 1) - 0.5)) * np.sqrt(6)
    else:
        noise = np.random.rand(2 ** nextpow2(n_points),) - 0.5
        noise = hilbert(noise) / noise.std() / np.sqrt(2)
        inds = noise.shape[0] - np.arange(n_points - 1, 1, step=-1)
        noise = noise[inds]
    return noise


def noisecg(n_points, a1=None, a2=None):
    """
    Generate analytic complex gaussian noise with mean 0.0 and variance 1.0.

    :param n_points: Length of the desired output signal.
    :param a1:
        Coefficients of the filter through which the noise is passed.
    :param a2:
        Coefficients of the filter through which the noise is passed.
    :type n_points: int
    :type a1: float
    :type a2: float
    :return: Analytic complex Gaussian noise of length n_points.
    :rtype: numpy.ndarray
    """
    assert n_points > 0
    if n_points <= 2:
        noise = (np.random.randn(n_points, 1) + 1j * np.random.randn(n_points, 1)) / np.sqrt(2)
    else:
        noise = np.random.normal(size=(2 ** nextpow2(n_points),))
        noise = hilbert(noise) / noise.std() / np.sqrt(2)
        noise = noise[len(noise) - np.arange(n_points - 1, -1, -1) - 1]
    return noise


if __name__ == "__main__":
    n = noisecg(128)
