import numpy as np


def sigmerge(x1, x2, ratio=0.0):
    """
    Add two signals with a specific energy ratio in decibels.

    :param x1: 1D numpy.ndarray
    :param x2: 1D numpy.ndarray
    :param ratio: Energy ratio in decibels.
    :type x1: numpy.ndarray
    :type x2: numpy.ndarray
    :type ratio: float
    :return: The merged signal
    :rtype: numpy.ndarray
    """
    assert x1.ndim == 1
    assert x2.ndim == 1
    assert type(ratio) in (float, int)
    ex1 = np.mean(np.abs(x1) ** 2)
    ex2 = np.mean(np.abs(x2) ** 2)
    h = np.sqrt(ex1 / (ex2 * 10 ** (ratio / 10.0)))
    sig = x1 + h * x2
    return sig
