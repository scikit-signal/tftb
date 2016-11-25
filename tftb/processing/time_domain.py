import numpy as np


def loctime(sig):
    """
    Compute the time localization characteristics.

    :param sig: input signal
    :type sig: numpy.ndarray
    :return: Average time center and time spreading
    :rtype: tuple
    :Example:
    >>> from tftb.generators import amgauss
    >>> x = amgauss(160, 80.0, 50.0)
    >>> tm, T = loctime(x)
    >>> print("%.2f" % tm)
    79.00
    >>> print("%.2f" % T)
    50.00
    """
    if sig.ndim > 2:
        if 1 not in sig.shape:
            raise TypeError
        else:
            sig = sig.ravel()
    sig2 = np.abs(sig**2)
    sig2 = sig2 / sig2.mean()
    t = np.arange(len(sig))
    tm = np.mean(t * sig2)
    T = 2 * np.sqrt(np.pi * np.mean(((t - tm)**2) * sig2))
    return tm, T
