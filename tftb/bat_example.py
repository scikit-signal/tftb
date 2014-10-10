from scipy.io import loadmat
import matplotlib.pyplot as plt
import os
import numpy as np

filename = "/Users/jaidevd/GitHub/scikit-signal/tftb/tftb-0.2/data"
bat = loadmat(os.path.join(filename,'bat.mat'))
bat = bat.get('bat').ravel()

t0 = np.linspace(0,2500./2304,2500)
t0 = t0[:len(bat)]

# Spectrum of the bat
dsp = np.fft.fftshift(np.abs(np.fft.fft(bat))**2)
f0 = np.arange(-1250,1250)*230.4/2500
f0 = f0[:len(dsp)]
plt.plot(f0,dsp.ravel())
plt.show()