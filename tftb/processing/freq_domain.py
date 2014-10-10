import numpy as np

def locfreq(sig):
    if sig.ndim > 1:
        if 1 not in sig.shape:
            raise TypeError
        else:
            sig = sig.ravel()
    N = float(len(sig))
    No2r = np.round(N/2)
    No2f = np.floor(N/2)
    Sig = np.fft.fft(sig)
    Sig2 = np.abs(Sig)**2
    Sig2 = Sig2/Sig2.mean()
    freqs = np.linspace(0,-1,-1./(No2f-1-No)+1,endpoint=True)/N
    fm = np.mean(freqs*Sig2)
    B = 2*np.sqrt(np.pi*np.mean(((freqs-fm)**2)*Sig2))
    return fm, B
