import numpy as np
from generators.frequency_modulated import fmlin
from matplotlib.pyplot import plot, grid, xlabel, ylabel, title, show, figure,axis

def energy_spectrum():
    sig1, iflaw = fmlin(128,0,0.5)
    plot(np.real(sig1))
    xlabel("Time")
    ylabel("Real part")
    title("Linear frequency modulation")
    dsp1 = np.fft.fftshift(np.abs(np.fft.fft(sig1))**2)
    figure()
    plot(np.arange(-64,64)/128.0,dsp1)
    axis([-0.5,0.5,0,400])
    grid()
    xlabel("Normalized frequency")
    ylabel("Squared modulus")
    title("Spectrum")
    show()


if __name__ == "__main__":
    energy_spectrum()

