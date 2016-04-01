from tftb.generators import fmlin, amgauss
import numpy as np
import matplotlib.pyplot as plt

y_nonstat, _ = fmlin(2048)  # Already analytic, no need of Hilbert transorm
y_nonstat *= amgauss(2048)
plt.plot(np.real(y_nonstat), np.imag(y_nonstat))
plt.xlabel("Real part")
plt.ylabel("Imaginary part")
plt.show()
