

import numpy as np
import matplotlib.pyplot as plt
import WDRT.shortTermExtreme as ecm

t_peaks = np.loadtxt('data/t.dat')
peaks = np.loadtxt('data/peaks.dat')/1000.
t_st = 1.*60.*60

f1,f2,ev = ecm.compare_methods(peaks, t_peaks, t_st, methods=[1,2,3,4,5],colors=['g','b','r','k','k'],lines=['-','-','-','-','--'])
plt.show()
