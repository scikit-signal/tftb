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
expect it to be stationary. Let's try and assert this experimentally. If we can
verify that the mean and variance of the signal are constant even at different
instants in time, we can conclude that the signal is stationary. The localized
mean of the signal can simply be calculated by convolving the signal with a
constant window function of adequate length.
