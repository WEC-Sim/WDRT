import numpy as np
import matplotlib.pyplot as plt
import WDRT.shortTermExtreme as ecm
import os

# load global peaks
t_peaks_file = os.path.join('data', 't.dat')
peaks_file   = os.path.join('data', 'peaks.dat')
t_peaks = np.loadtxt(t_peaks_file)
peaks = np.loadtxt(peaks_file)/1000.

t_st = 1. * 60. * 60

f1, f2, ev = ecm.compare_methods(peaks, t_peaks, t_st, methods=[1, 2, 3, 4, 5],
                                 colors=['g', 'b', 'r', 'k', 'k'],
                                 lines=['-', '-', '-', '-', '--'])
plt.show()
