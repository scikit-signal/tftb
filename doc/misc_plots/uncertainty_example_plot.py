import numpy as np
import matplotlib.pyplot as plt
f1, f2 = 500, 1000
t1, t2 = 0.192, 0.196
f_sample = 8000
n_points = 2048
ts = np.arange(n_points, dtype=float) / f_sample
signal = np.sin(2 * np.pi * f1 * ts) + np.sin(2 * np.pi * f2 * ts)
signal[int(t1 * f_sample) - 1] += 3
signal[int(t2 * f_sample) - 1] += 3
plt.plot(ts, signal)
plt.show()
