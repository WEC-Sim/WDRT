import numpy as np
import scipy.io
import WDRT.ESSC as ESSC
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import scipy.interpolate as interp
import WDRT.longTermExtreme as lte
import WDRT.shortTermExtreme as ste
import h5py
import os


# Load data from example_envSampling.py
#envFile = h5py.File(os.path.join(r'data', 'NDBC46022.h5'), 'r')
envFile = h5py.File(os.path.join(r'testNDBC46022.h5'), 'r')
Hs_Return = np.array(envFile['ReturnContours/Hs_Return'])
T_Return = np.array(envFile['ReturnContours/T_Return'])
Hs_sample = np.array(envFile['Samples_ContourApproach/Hs_SampleCA'])
T_sample = np.array(envFile['Samples_ContourApproach/T_SampleCA'])

# Load data from modeling
modResFile = h5py.File(os.path.join(
    r'data', 'longTerm_contourApproach.h5'), 'r')
t = np.array(modResFile['time'])
tSim = t[-1]
n = len(modResFile['x'])
peaks = []
mmax = []
x = []
for ii in range(n):
    ss = 'ss_%03d' % ii
    x.append(np.array(modResFile['/x/' + ss]))
    peaks.append(ste.globalPeaks(t, x[ii])[1])
    mmax.append(np.max(peaks[ii]))

# Short-term extreme response at each sea state
tPer = 1 * 60 * 60  # storm period
x_t = np.linspace(0, 1.5 * np.max(mmax), 100)
edist = []
ev = []
r975 = []
for ii in range(n):
    edist.append(ste.extremeDistribution_WeibullTailFit(
        x=peaks[ii], x_e=x_t, t_x=tSim, t_st=tPer)[0])
    ev.append(edist[ii].getExpVal())

# Find largest response
mi = np.argmax(ev)
x0 = edist[mi].getRthVal(0.00001)
x1 = edist[mi].getRthVal(1 - 0.00001)
x = np.linspace(x0, x1, 500)
r95 = edist[mi].getRthVal(0.95)
print('design state (Hs, Te): (%.1f, %.1f)' % (Hs_sample[mi], T_sample[mi]))
print('extreme value: %e' % (r95))

# Plot data
plt.figure()
plt.plot(T_Return, Hs_Return, 'k-', label='100 year contour')
plt.plot(T_sample, Hs_sample, 'y^', label='full sea state samples')
plt.legend(loc='lower right', fontsize='small')
plt.grid(True)
plt.xlabel('Energy period, $T_e$ [s]')
plt.ylabel('Sig. wave height, $H_s$ [m]')

fig = plt.figure()
ax = fig.add_subplot(111, projection='3d', autoscale_on=True)
ax.plot(T_Return, Hs_Return, np.zeros_like(T_Return), 'k')
ax.plot(T_sample, Hs_sample, np.zeros_like(Hs_sample), 'yo')
for i in range(len(T_sample)):
    ax.plot([T_sample[i], T_sample[i]], [Hs_sample[i], Hs_sample[i]],
            [0, ev[i]], '-', linewidth=2, color='b', alpha=0.5)
ax.plot(T_sample, Hs_sample, ev, '-o')
ax.set_xlabel('Energy period, $T_e$ (s)')
ax.set_ylabel('Sig. wave height, $H_s$ (m)')
ax.set_zlabel('$x$')

plt.figure()
plt.plot(x, edist[mi].cdf(x))
plt.plot(plt.xlim(), 0.95 * np.ones(2), 'r--')
plt.plot(r95 * np.ones(2), plt.ylim(), 'r--')
plt.grid(True)
plt.xlabel('$x$')
plt.ylabel('$CDF(x)$')

plt.show()
