import numpy as np


def amgauss(N, t0=None, T=None):
    if t0 is None:
        t0 = np.round(N/2)

    if T is None:
        T = 2*np.sqrt(N);


    if N <= 0:
        raise TypeError
    else:
        tmt0 = np.arange(N) - t0;
        y = np.exp(-((tmt0/T)**2)*np.pi)
        return y


def amexpos(N, kind="unilateral", t0=None, T=None):
    if t0 is None:
        t0 = np.round(N/2)
    if T is None:
        T = 2*np.sqrt(N)

    if N <= 0:
        raise TypeError
    else:
        tmt0 = np.arange(N) - t0
        if kind == "bilateral":
            y = np.exp(-np.sqrt(2*np.pi)*np.abs(tmt0)/T)
        else:
            y = np.exp(-np.sqrt(np.pi)*tmt0/T)*(tmt0>=0.0)
        return y

def amrect(N, t0=None, T=None):
    if t0 is None:
        t0 = np.round(N/2)
    if T is None:
        T = 2*np.sqrt(N)

    if N <= 0:
        raise TypeError
    else:
        tmt0 = np.arange(N) - t0
        y = np.abs(tmt0) <= 0.5*T*np.sqrt(3.0/np.pi)
        return y

def amtriang(N, t0=None, T=None):
    if t0 is None:
        t0 = np.round(N/2)
    if T is None:
        T = 2*np.sqrt(N)

    if N <= 0:
        raise TypeError
    else:
        tmt0 = np.arange(N) - t0
        L = np.sqrt(10.0/np.pi)*T/2.0
        t = np.amin(np.vstack((L+tmt0, L-tmt0)).T, axis=1)
        t = np.hstack((t.reshape((len(t),1)),np.zeros((len(t),1))))
        y = np.amax( t, axis=1)/L

        
        return y

if __name__ == "__main__":
    from matplotlib.pyplot import plot, show
    y = amtriang(64)
    plot(y)
    show()
