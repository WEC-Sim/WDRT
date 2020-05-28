import numpy as np
import WDRT.shortTermExtreme as ste
import WDRT.longTermExtreme as lte
import matplotlib.pyplot as plt
import h5py
import os


# Load data from example_envSampling.py
envFile = h5py.File(os.path.join(r'data', 'NDBC46022.h5'), 'r')
Hs_sample = np.array(envFile['Samples_FullSeaState/Hs_SampleFSS'])
T_sample = np.array(envFile['Samples_FullSeaState/T_SampleFSS'])
Weight_sample = np.array(envFile['Samples_FullSeaState/Weight_SampleFSS'])
Hs = np.array(envFile['buoy_Data/Hs'])
T = np.array(envFile['buoy_Data/Te'])
Hs_Return = np.array(envFile['ReturnContours/Hs_Return'])
T_Return = np.array(envFile['ReturnContours/T_Return'])

# Load data from modeling (if the cwd is not examples adjust the data path)
modResFile = h5py.File(os.path.join('data', 'longTerm_FullSeaState.h5'), 'r')
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
ccdf = []
for ii in range(n):
    edist.append(ste.extremeDistribution_WeibullTailFit(
        x=peaks[ii], x_e=x_t, t_x=tSim, t_st=tPer)[0])
    ccdf.append(edist[ii].ccdf)

# Long-term extreme response
LTS = lte.fullLongTermSurvival(Fr=ccdf, fs=Weight_sample)

# Plotting
plt.figure().canvas.set_window_title('Sea state sampling')
plt.plot(T, Hs, 'bo', alpha=0.25, label='data')
plt.plot(T_sample, Hs_sample, 'ro', label='samples')
plt.plot(T_Return, Hs_Return, label='100 year contour')
plt.legend(loc='best')
plt.grid(True)
plt.xlabel('Energy period, $T_e$ [s]')
plt.ylabel('Sig. wave height, $H_s$ [m]')

plt.figure().canvas.set_window_title('Long-term response')
plt.plot(x_t, LTS(x_t), 'k-', label='Full sea state survival')
for ii in range(n):
    if ii == 0:
        plt.plot(x_t, ccdf[ii](x_t), 'r-', alpha=0.1, label='Individual sea state survivals')
    else:
        plt.plot(x_t, ccdf[ii](x_t), 'r-', alpha=0.1)
plt.gca().set_yscale('log')
plt.grid(True)
yLs = [1, 5, 10, 25, 50, 100]
for yL in yLs:
    plt.plot(plt.xlim(), (1 / (yL * 365.25 * 24 * 1)) * np.ones(2), 'r--')
    plt.text(plt.xlim()[0] + 0.05 * plt.xlim()[1],
             (1 / (yL * 365.25 * 24 * 1)),
             '%i-year return' % (yL), fontsize=10)
plt.xlabel('$x$')
plt.ylabel('$CCDF(x)$')
plt.xlim((0, 30))
plt.ylim((1e-1 / (yLs[-1] * 365.25 * 24 * 1), 1))
plt.legend(loc='upper right')

plt.show()
