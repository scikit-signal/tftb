import numpy as np
import matplotlib.pyplot as plt
from tftb.processing.linear import ShortTimeFourierTransform

# from tftb.processing import ShortTimeFourierTransform
from matplotlib.pyplot import cm

f1, f2 = 500, 1000
t1, t2 = 0.192, 0.196
f_sample = 8000
n_points = 2048
ts = np.arange(n_points, dtype=float) / f_sample
signal = np.sin(2 * np.pi * f1 * ts) + np.sin(2 * np.pi * f2 * ts)
signal[int(t1 * f_sample) - 1] += 3
signal[int(t2 * f_sample) - 1] += 3

wlengths = [2, 4, 8, 16]
nf = [(w * 0.001 * f_sample) + 1 for w in wlengths]
fig = plt.figure(figsize=(10, 8))
extent = [0, ts.max(), 0, 2000]
for i, wlen in enumerate(wlengths):
    window = np.ones((int(nf[i]),), dtype=float)

    n_fbins = n_points
    nperseg = int(nf[i])
    noverlap = nperseg - 1
    nfft = n_points
    stft = ShortTimeFourierTransform(signal, timestamps=None, n_fbins=n_fbins)
    tfr, ts, freqs = stft.run(
        nfft=nfft,
        nperseg=nperseg,
        noverlap=noverlap,
        return_onesided=False,
        window=window,
        scaling="psd")

    ax = fig.add_subplot(4, 1, i + 1)
    # stft.plot(ax=ax, show_tf=True, cmap=cm.gray)

    # stft = ShortTimeFourierTransform(signal, fwindow=window)
    # stft.run()
    stft.plot(ax=ax, default_annotation=False, show=False,
              extent=extent)
    ax.set_yticklabels([])
    if i != 3:
        ax.set_xticklabels([])
    ax.set_title("window duration = {} ms".format(wlen))
plt.subplots_adjust(hspace=0.5)
plt.show()
