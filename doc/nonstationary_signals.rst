======================
Non Stationary Signals
======================

A Note on Stationarity
----------------------

As per Wikipedia, a "stationary" process is one whose joint probability
distribution does not change with time (or space). Let's try and see what a
stationary process looks like. Consider a signal generated as follows:

    >>> import numpy as np
    >>> import matplotlib.pyplot as plt
    >>> fs = 32768
    >>> ts = np.linspace(0, 1, fs)
    >>> y1 = np.sin(2 * np.pi * 697 * ts)
    >>> y2 = np.sin(2 * np.pi * 1336 * ts)
    >>> y = (y1 + y2) / 2
    >>> plt.plot(ts, y)
    >>> plt.xlim(0, 0.1)
    >>> plt.show()

    .. plot:: misc_plots/touchtone.py

The plot shows a slice of a touchtone signal with the duration of a tenth of a
second. This is the signal you hear when you press the "2" key on a telephone
dialing pad. (You can save the generated signal as a WAV file as follows:

    >>> from scipy.io import wavfile
    >>> wavfile.write("tone.wav", fs, y)

and listen to the file ``tone.wav`` with your favourite music player.)

Since the signal is composed of two sinusoids, ``y1`` and ``y2``, we would
expect it to be stationary. Let's try and assert this qualitatively. Let's try
to plot the signal in its phase space. In order to do so, we will first need to
construct an analytic representation of the signal. This can be achieved by
taking the Hilbert transform of the signal. For simplicity, we shall only
consider a part of the original signal


    >>> y = y[:(fs / 16)]
    >>> y_analytic = hilbert(y)
    >>> plt.plot(np.real(y_analytic), np.imag(y_analytic))
    >>> plt.xlabel("Real part")
    >>> plt.ylabel("Imaginary part")
    >>> plt.show()
    
    .. plot:: misc_plots/stationary_phase_plot.py

This visualization can be interpreted as follows. Imagine that there is a
vector centered at the origin of this plot which traces out the signal as it
rotates about the origin. Then, at any time :math:`t`, the angle which the
vector makes with the real axis is the instantaneous phase of the signal,
:math:`\theta(t)`. The angular speed with which the phasor rotates is the
instantaneous frequency of the signal:

    .. math::

      \omega(t) = \frac{d}{dt} \theta(t)

Now, let's compare this phase plot with that of a known nonstationary signal.


    >>> from tftb.generators import fmlin, amgauss
    >>> y_ns, _ = fmlin(2048)  # Already analytic, no need of Hilbert transorm
    >>> y_nonstat = y_ns * amgauss(2048)
    >>> plt.plot(np.real(y_nonstat), np.imag(y_nonstat))
    >>> plt.xlabel("Real part")
    >>> plt.ylabel("Imaginary part")
    >>> plt.show()

    .. plot:: misc_plots/nonstationary_phase_plot.py

Notice that the second plot has a lot of rough edges and sharp angles. This
means that when the signal vector rotates through the phase space, it will have
to make sharp jumps in order to trace the signal. Moreover, by the definition
of instantaneous frequency, when the instantaneous phase is not differentiable,
the instantaneous frequency will be indeterminate. By contrast, the phase plot
for the stationary signal is a lot smoother, and we can expect that the
instantaneous frequency will be finite at all times. Physically, this means
that in the nonstationary signal, the variation in frequency components has no
structure, and these components can change arbitrarily.

This phenomenon of arbitrary, unstructured changes in frequency over time is a
symptom of nonstationarity, and will become increasingly relevant as we
proceed.
