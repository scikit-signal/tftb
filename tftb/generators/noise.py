import numpy as np
from tftb.utils import nextpow2
from scipy.signal import hilbert


def noisecu(n_points):
    """Compute analytic complex uniform white noise.

    :param n_points: Length of the noise signal.
    :type n_points: int
    :return: analytic complex uniform white noise signal of length N
    :rtype: numpy.ndarray
    :Examples:
    >>> import matplotlib.pyplot as plt
    >>> import numpy as np
    >>> noise = noisecu(512)
    >>> print("%.2f" % abs((noise ** 2).mean()))
    0.00
    >>> print("%.1f" % np.std(noise) ** 2)
    1.0
    >>> plt.subplot(211), plt.plot(real(noise))                                #doctest: +SKIP
    >>> plt.subplot(212),  #doctest: +SKIP
    >>> plt.plot(linspace(-0.5, 0.5, 512), abs(fftshift(fft(noise))) ** 2) #doctest: +SKIP

    .. plot:: docstring_plots/generators/noise/noisecu.py
    """
    if n_points <= 2:
        noise = (np.random.rand(n_points, 1) - 0.5 + 1j * (np.random.rand(n_points, 1) - 0.5)) * \
            np.sqrt(6)
    else:
        noise = np.random.rand(2 ** int(nextpow2(n_points)),) - 0.5
        noise = hilbert(noise) / noise.std() / np.sqrt(2)
        inds = noise.shape[0] - np.arange(n_points - 1, -1, step=-1) - 1
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
    :Examples:
    >>> import matplotlib.pyplot as plt
    >>> import numpy as np
    >>> noise = noisecg(512)
    >>> print("%.1f" % abs((noise ** 2).mean()))
    0.0
    >>> print("%.1f" % np.std(noise) ** 2)
    1.0
    >>> plt.subplot(211), plt.plot(real(noise))                                #doctest: +SKIP
    >>> plt.subplot(212), #doctest: +SKIP
    >>> plt.plot(linspace(-0.5, 0.5, 512), abs(fftshift(fft(noise))) ** 2) #doctest: +SKIP

    .. plot:: docstring_plots/generators/noise/noisecg.py
    """
    assert n_points > 0
    if n_points <= 2:
        noise = (np.random.randn(n_points, 1.) + 1j * np.random.randn(n_points, 1.)) / np.sqrt(2.)
    else:
        noise = np.random.normal(size=int(2. ** nextpow2(float(n_points)),))
        noise = hilbert(noise) / noise.std() / np.sqrt(2.)
        noise = noise[len(noise) - np.arange(n_points - 1, -1, -1) - 1]
    return noise


def dopnoise(n_points, s_freq, f_target, distance, v_target,
             time_center=None, c=340):
    """Generate complex noisy doppler signal, normalized to have unit energy.

    :param n_points: Number of points.
    :param s_freq: Sampling frequency.
    :param f_target: Frequency of target.
    :param distance: Distnace from line to observer.
    :param v_target: velocity of target relative to observer.
    :param time_center: Time center. (Default n_points / 2)
    :param c: Wave velocity (Default 340 m/s)
    :type n_points: int
    :type s_freq: float
    :type f_target: float
    :type distance: float
    :type v_target: float
    :type time_center: float
    :type c: float
    :return: tuple (output signal, instantaneous frequency law.)
    :rtype: tuple(array-like)
    :Example:
    >>> import numpy as np
    >>> from tftb.processing import inst_freq
    >>> z, iflaw = dopnoise(500, 200.0, 60.0, 10.0, 70.0, 128.0)
    >>> subplot(211), plot(real(z))                #doctest: +SKIP
    >>> ifl = inst_freq(z, np.arange(11, 479), 10)
    >>> subplot(212), plot(iflaw, 'r', ifl, 'g')   #doctest: +SKIP

    .. plot:: docstring_plots/generators/noise/dopnoise.py
    """
    if time_center is None:
        time_center = np.floor(n_points / 2.0)

    r = 0.9
    rr = r ** 2
    r2 = 2 * r
    vv = v_target ** 2
    x = np.random.randn(2 * n_points,)
    tmt0 = (np.arange(1, 2 * n_points + 1, dtype=float) - time_center - n_points) / s_freq
    dist = np.sqrt(distance ** 2 + (v_target * tmt0) ** 2)
    iflaw = (1 - vv * tmt0 / dist / c) * f_target / s_freq
    y = np.zeros((2 * n_points,))
    for t in range(2, 2 * n_points):
        y[t] = x[t] - rr * (x[t - 2] + y[t - 2]) + r2 * np.cos(2.0 * np.pi * iflaw[t]) * y[t - 1]
    y = hilbert(y[(n_points + 1): (2 * n_points + 1)]) /\
        np.sqrt(dist[(n_points + 1): (2 * n_points + 1)])
    y = y / np.sqrt(np.sum(np.abs(y) ** 2))
    iflaw = iflaw[(n_points + 1):(2 * n_points + 1)]
    return y, iflaw


if __name__ == "__main__":
    n = noisecg(128)
