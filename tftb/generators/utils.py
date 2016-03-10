import numpy as np
from scipy.signal import hilbert
from scipy.integrate import trapz
# from tftb.utils import nextpow2


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


def scale(X, a, fmin, fmax, N):
    """Scale a signal with the Mellin transform.

    :param X: signal to be scaled.
    :param a: scale factor
    :param fmin: lower frequency bound
    :param fmax: higher frequency bound
    :param N: number of analyzed voices
    :type X: array-like
    :type a: float
    :type fmin: float
    :type fmax: float
    :type N: int
    :return: A-scaled version of X.
    :rtype: array-like
    """
    Z = hilbert(np.real(X))
    T = X.shape[0]
    M = (T + np.remainder(T, 2)) / 2

    # B = fmax - fmin
    # R = B / ((fmin + fmax) / 2)
    # Nq = np.ceil((B * T * (1 + 2.0 / R)) * np.log())
    # Nmin = Nq - np.remainder(Nq, 2)
    # Ndflt = 2 ** nextpow2(Nmin)

    # Geometric sampling of the analyzed spectrum
    k = np.arange(N)
    q = (fmax / float(fmin)) ** (1.0 / (N - 1.0))
    geo_f = fmin * np.exp((k - 1) * np.log(q))
    t = np.arange(X.shape[0]) - M - 1
    tfmatx = np.exp(-2 * 1j * np.dot(t.reshape(128, 1), geo_f.reshape(1, 128)) * np.pi)
    ZS = np.dot(Z, tfmatx)
    ZS = np.hstack((ZS, np.zeros((N,))))

    # Mellin transform computation of the analyzed signal
    p = np.arange(2 * N)
    MS = np.fft.fftshift(np.fft.ifft(ZS, axis=0))
    beta = (p / float(N) - 1.0) / (2 * np.log(q))

    # Inverse mellin and fourier transform.
    Mmax = np.amax(np.ceil(X.shape[0] / 2.0 * a))
    if isinstance(a, np.ndarray):
        S = np.zeros((2 * Mmax, a.shape[0]), dtype=complex)
    else:
        S = np.zeros((2 * Mmax,), dtype=complex)
        ptr = 0
        DMS = np.exp(-2 * np.pi * 1j * beta * np.log(a)) * MS
        DS = np.fft.fft(np.fft.fftshift(DMS), axis=0)
        Mcurrent = np.ceil(a * X.shape[0] / 2)
        t = np.arange(-Mcurrent, Mcurrent) - 1
        itfmatx = np.exp(2 * 1j * np.pi * np.dot(t.reshape((256, 1)),
                                                 geo_f.reshape((1, 128))))
        dilate_sig = np.zeros((2 * Mcurrent,), dtype=complex)
        for kk in range(2 * int(Mcurrent)):
            dilate_sig[kk] = trapz(itfmatx[kk, :] * DS[:N], geo_f)
        S[(Mmax - Mcurrent):(Mmax + Mcurrent)] = dilate_sig
        ptr += 1

    S = S * np.linalg.norm(X) / np.linalg.norm(S)
    return S

if __name__ == '__main__':
    from tftb.generators import altes
    sig = altes(256, 0.1, 0.45, 10000)
