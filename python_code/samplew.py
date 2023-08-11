import numpy as np
from scipy.special import i0, i1
from scipy.interpolate import interp1d

def samplew(size, kappa):
    u = np.random.rand(size)
    w = np.arange(-1, 1, 0.0001)
    first_two_terms = 1 + (w * np.sqrt(1 - w**2) - np.arccos(w)) / np.pi
    mplus2_term = lambda m: (m + 1) / np.pi * (np.sin((m + 2) * np.arccos(w)) / (m + 2) - np.sin(m * np.arccos(w)) / m) * i1(m + 1, kappa) / i1(1, kappa)
    Fw = first_two_terms
    mvalue = 1
    if kappa:
        while True:
            newterm = mplus2_term(mvalue)
            newterm_magnitude = np.sqrt(np.sum(newterm**2))
            if np.isnan(newterm_magnitude):
                raise ValueError("Nan values detected, kappa is too high!")
            elif abs(newterm_magnitude) < 1e-10:
                break
            Fw += newterm
            mvalue += 1
    Fw[np.abs(Fw) < 1e-10] = 0
    if not np.all(np.diff(Fw) >= 0):
        raise ValueError("Cumulative function is not monotonic!")
    smallseries = np.linspace(0, 1e-15, len(w))
    data = np.column_stack((w, (Fw + smallseries) / (Fw[-1] + smallseries[-1])))
    sampledvalues = interp1d(data[:, 1], data[:, 0])(u)
    return sampledvalues
