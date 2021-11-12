import sys
import math
import numpy as np
import matplotlib.pyplot as plt
from scipy import signal
from scipy import fft



def Fourier_filtering():
    return

def FFT (sig):
    N = len(sig)
    sig_ft = fft(sig)[1:N/2]
    periods = [float(N)/k for k in range(1, N/2)]
    amplts = np.abs(sig_ft)/float(N)
    phases = np.arctan2(sig_ft.imag, sig_ft.real) / np.pi # unit of pi
    shifts = np.asarray([(phases[k-1] * N) / (2*np.pi*k) for k in range(1, N/2)]) /np.pi # unit of pi
    return periods, amplts, phases
