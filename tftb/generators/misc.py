import numpy as np
from numpy import pi
from math import sqrt
from tftb.generators.amplitude_modulated import amgauss
from tftb.generators.frequency_modulated import fmconst


def altes(n_points, fmin=0.05, fmax=0.5, alpha=300):
    """Generate the Altes signal in time domain.

    :param n_points: Number of points in time.
    :param fmin: Lower frequency bound.
    :param fmax: Higher frequency bound.
    :param alpha: Attenuation factor of the envelope.
    :type n_points: int
    :type fmin: float
    :type fmax: float
    :type alpha: float
    :return: Time vector containing the Altes signal samples.
    :rtype: numpy.ndarray
    :Examples:
    >>> x = altes(128, 0.1, 0.45)
    >>> plot(x) #doctest: +SKIP

    .. plot:: docstring_plots/generators/misc/altes.py
    """
    g = np.exp((np.log(fmax / fmin)) ** 2 / (8 * np.log(alpha)))
    nu0 = np.sqrt(fmin * fmax)
    beta = np.sqrt(2 * np.log(g) * np.log(alpha))
    t0 = n_points / (np.exp(beta) - np.exp(-beta))
    t1 = t0 * np.exp(-beta)
    t2 = t0 * np.exp(beta)
    b = -t0 * nu0 * g * np.log(g)
    t = np.linspace(t1, t2, n_points + 1)[:n_points]
    x = np.exp(-(np.log(t / t0) ** 2) / (2 * np.log(g))) * \
        np.cos(2 * pi * b * np.log(t / t0) / np.log(g))
    x = x / np.linalg.norm(x)
    return x


def atoms(n_points, coordinates):
    """Compute linear combination of elementary Gaussian atoms.

    :param n_points: Number of points in a signal
    :param coordinates: matrix of time-frequency centers
    :type n_points: int
    :type coordinates: array-like
    :return: signal
    :rtype: array-like
    :Examples:
    >>> import numpy as np
    >>> coordinates = np.array([[32.0, 0.3, 32.0, 1.0], [0.15, 0.15, 48.0, 1.22]])
    >>> sig = atoms(128, coordinates)
    >>> plot(real(sig)) #doctest: +SKIP

    .. plot:: docstring_plots/generators/misc/atoms.py
    """
    if coordinates.ndim == 1:
        coordinates = coordinates.reshape((coordinates.shape[0], 1))
    signal = np.zeros((n_points,), dtype=complex)
    n_atoms = coordinates.shape[0]
    for k in range(n_atoms):
        t0 = round(np.max((np.min((coordinates[k, 0], n_points)), 1)))
        f0 = np.max((np.min((coordinates[k, 1], 0.5)), 0.0))
        T = coordinates[k, 2]
        A = coordinates[k, 3]
        signal += A * amgauss(n_points, t0, T) * fmconst(n_points, f0, t0)[0]
    return signal


def doppler(n_points, s_freq, f0, distance, v_target, t0=None, v_wave=340.0):
    """Generate complex Doppler signal

    :param n_points: Number of points
    :param s_freq: Sampling frequency
    :param f0: Target frequency
    :param distance: distance from line to observer
    :param v_target: Target velocity
    :param t0: Time center
    :param v_wave: wave velocity.
    :type n_points: int
    :type s_freq: float
    :type f0: float
    :type distance: float
    :type v_target: float
    :type t0: float
    :type v_wave: float
    :return: Tuple containing output frequency modulator, output amplitude \
        modulator, output instantaneous frequency law.
    :rtype: tuple
    :Example:
    >>> fm, am, iflaw = doppler(512, 200.0, 65.0, 10.0, 50.0)
    >>> subplot(211), plot(real(am * fm)) #doctest: +SKIP
    >>> subplot(211), plot(iflaw)         #doctest: +SKIP

    .. plot:: docstring_plots/generators/misc/doppler.py
    """
    if t0 is None:
        t0 = n_points / 2

    if distance <= 0.0:
        raise TypeError("distance must be strictly positive.")
    elif s_freq < 0.0:
        raise TypeError("Sampling frequency must be positive.")
    elif (t0 < 1) or (t0 > n_points):
        raise TypeError("T0 must be between 1 and n_points")
    elif (f0 < 0) or (f0 > s_freq / 2):
        raise TypeError("F0 must be between 0 and s_freq/2")
    elif v_target < 0:
        raise TypeError("v_target must be positive")

    tmt0 = (np.arange(1, n_points + 1) - t0) / s_freq
    dist = np.sqrt(distance ** 2 + (v_target * tmt0) ** 2)
    fm = np.exp(1j * 2 * pi * f0 * (tmt0 - dist / v_wave))
    if np.abs(f0) < np.spacing(1):
        am = 0
    else:
        am = 1. / np.sqrt(dist)
    iflaw = (1 - v_target ** 2 * tmt0 / dist / v_wave) * f0 / s_freq
    return fm, am, iflaw


def klauder(n_points, attenuation=10.0, f0=0.2):
    """Klauder wavelet in time domain.

    :param n_points: Number of points in time.
    :param attenuation: attenuation factor of the envelope.
    :param f0: central frequency of the wavelet.
    :type n_points: int
    :type attenuation: float
    :type f0: float
    :return: Time row vector containing klauder samples.
    :rtype: numpy.ndarray
    :Example:
    >>> x = klauder(128)
    >>> plot(x) #doctest: +SKIP

    .. plot:: docstring_plots/generators/misc/klauder.py
    """

    assert n_points > 0
    assert ((f0 < 0.5) and (f0 > 0))

    f = np.linspace(0., 0.5, int(n_points / 2 + 1))
    mod = np.exp(-2. * np.pi * attenuation * f) * f ** (2. * np.pi * attenuation * f0 - 0.5)
    wave = mod
    wave[0] = 0
    a, b = wave[:int(n_points / 2)], wave[1:int(n_points / 2 + 1)][::-1]
    wave = np.hstack((a, b))
    wavet = np.fft.ifft(wave)
    wavet = np.fft.fftshift(wavet)
    x = np.real(wavet) / np.linalg.norm(wavet)
    return x


def mexhat(nu=0.05):
    """Mexican hat wavelet in time domain.

    :param nu: Central normalized frequency of the wavelet. Must be a real \
        number between 0 and 0.5
    :type nu: float
    :return: time vector containing mexhat samples.
    :rtype: numpy.ndarray
    :Example:
    >>> plot(mexhat()) #doctest: +SKIP

    .. plot:: docstring_plots/generators/misc/mexhat.py
    """
    assert (nu <= 0.5) and (nu >= 0)
    n_points = 1.5
    alpha = np.pi ** 2 * nu ** 2
    n = np.ceil(n_points / nu)
    t = np.arange(-n, n + 1, step=1)
    h = nu * sqrt(pi) / 2 * np.exp(-alpha * t ** 2) * (1 - 2 * alpha * t ** 2)
    return h


def gdpower(n_points, degree=0.0, rate=1.0):
    """Generate a signal with a power law group delay.

    :param n_points: Number of points in time.
    :param degree: degree of the power law.
    :param rate: rate-coefficient of the power law GD.
    :type n_points: int
    :type degree: int		# yoder: i'm pretty sure this is supposed to be a float. at least in test_misc.py, the test for this funct. passes degree=0.5
    :type rate: float
    :return: Tuple of time row containing modulated samples, group delay, \
            frequency bins.
    :rtype: tuple
    """
    n_points = int(n_points)
    degree = float(degree)
    rate = float(rate)
    t0 = 0
    lnu = int(np.ceil(n_points / 2.0))
    nu = np.linspace(0, 0.5, lnu + 1)
    nu = nu[1:]
    am = nu ** ((degree - 2.0) / 6.0)

    if rate == 0.:
        raise TypeError("rate must be non-zero")

    tfx = np.zeros((n_points,), dtype=complex)

    if (degree < 1.) and (degree != 0):
        d = float(n_points) ** (degree * rate)
        t0 = n_points / 10.0
        tfx[:lnu] = np.exp(-1.0j * 2. * pi * (t0 * nu + d * nu ** degree / degree)) * am
        x = np.fft.ifft(tfx)
    elif degree == 0:
        d = rate
        t0 = n_points / 10.0
        tfx[:lnu] = np.exp(-1j * 2 * np.pi * (t0 * nu + d * np.log(nu))) * am
        x = np.fft.ifft(tfx)
    elif degree == 1.:
        from tftb.generators.analytic_signals import anapulse
        t0 = n_points
        x = anapulse(n_points, t0)
    elif degree > 1.:
        d = n_points * 2. ** (degree - 1.) * rate
        tfx[:lnu] = np.exp(-1.0j * 2.0 * pi * (t0 * nu + d * nu ** degree / degree)) * am
        x = np.fft.ifft(tfx)
    else:
        t0 = n_points / 10
        d = n_points * 2.0 ** (degree - 1.0) * rate
        tfx[:lnu] = np.exp(-1j * 2 * pi * (t0 * nu + d * np.log(nu))) * am
        x = np.fft.ifft(tfx)

    if degree != 1.0:
        gpd = t0 + np.abs(np.sign(rate) - 1.0) / 2.0 * (n_points + 1.0) + d * nu ** (degree - 1.0)
    else:
        gpd = t0 * np.ones((n_points / 2.0,))

    x = x - x.mean()
    x = x / np.linalg.norm(x)

    return x, gpd, nu

if __name__ == "__main__":
    altes(128)
    coordinates = np.array([[32.0, 0.3, 32.0, 1.0],
                            [0.15, 0.15, 48.0, 1.22]])
    atoms(128, coordinates)
    doppler(128, 200.0, 65.0, 10.0, 50.0)
    klauder(128)
    mexhat()
    gdpower(128)
