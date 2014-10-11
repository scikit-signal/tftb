import numpy as np


def nextpow2(n):
    m_f = np.log2(n)
    m_i = np.ceil(m_f)
    return m_i


def scale():
    pass


def sigmerge(x1, x2, ratio=0.0):
    """
    Add two signals with a specific energy ratio in decibels.

    Parameters
    ----------
    x1, x2 : 1-D ndarrays
        Signals to be added.

    ratio : float, optional
        Energy ratio in decibels.

    Returns
    -------
    sig : 1-D ndarray
        The merged signal

    """
    assert x1.ndim == 1
    assert x2.ndim == 1
    assert type(ratio) in (float, int)
    ex1 = np.mean(np.abs(x1)**2)
    ex2 = np.mean(np.abs(x2)**2)
    h = np.sqrt(ex1/(ex2*10**(ratio/10.0)))
    sig = x1 + h * x2
    return sig


if __name__ == "__main__":
    print nextpow2(128)
