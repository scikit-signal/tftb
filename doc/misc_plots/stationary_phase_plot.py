import numpy as np
import matplotlib.pyplot as plt
from scipy.signal import hilbert

fs = 32768
ts = np.linspace(0, 1, fs)
y1 = np.sin(2 * np.pi * 697 * ts)
y2 = np.sin(2 * np.pi * 1336 * ts)
y = (y1 + y2) / 2


y = y[:int(fs / 16)]
y_analytic = hilbert(y)
plt.plot(np.real(y_analytic), np.imag(y_analytic))
plt.xlabel("Real part")
plt.ylabel("Imaginary part")
plt.show()
