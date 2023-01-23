from tftb.processing.linear import ShortTimeFourierTransform
from tftb.generators import fmconst
import numpy as np
from scipy.signal import hamming
from matplotlib.pyplot import cm

sig = np.r_[fmconst(128, 0.2)[0], fmconst(128, 0.4)[0]]

nperseg = 65
window = hamming(nperseg)
n_fbins = 128
noverlap = nperseg - 1
nfft = 128
stft = ShortTimeFourierTransform(sig, timestamps=None, n_fbins=n_fbins)
tfr, ts, freqs = stft.run(
    nfft=nfft,
    nperseg=nperseg,
    noverlap=noverlap,
    return_onesided=False,
    window=window,
    scaling="psd")
stft.plot(show_tf=True, cmap=cm.gray)
