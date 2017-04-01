
import numpy as np
import matplotlib.pyplot as plt
import WDRT.shortTermExtreme as ecm
import WDRT.fatigue as fatigue

# load response time series
data = ecm.loadtxt('data/data.csv', delimiter=',')
t = data['t']
response = data['data']

# find global peaks
t_peaks, peaks = ecm.globalPeaks(t, response)

# plot
plt.figure()
plt.hold(True)
plt.plot(t, response, 'k-')
plt.plot(t_peaks, peaks, 'go')
plt.plot([0, t[-1]], [0, 0], 'k--')
plt.xlabel('Time, $t$ [s]')
plt.ylabel('Response, $x$')
plt.grid(True)
plt.ticklabel_format(style='sci', axis='y', scilimits=(0, 0))

# get the 1-hour extreme distribution using the Weibull tail fit method
x_e = np.linspace(0, 2 * np.max(peaks), 10000)
t_x = (t[-1] - t[0])
t_st = 1. * 60. * 60.
stextreme_dist, peaks_dist, _, _, _ = ecm.extremeDistribution_WeibullTailFit(x=peaks, x_e=x_e, t_x=t_x, t_st=t_st)

# plot
plt.figure()
ax = plt.subplot(2, 1, 1)
plt.hold(True)
plt.plot(x_e, peaks_dist.pdf(x_e), 'g-', label='Peak distribution')
plt.plot(x_e, stextreme_dist.pdf(x_e), 'r-', label='Extreme distribution')
xlim = ax.get_xlim()
ylim = ax.get_ylim()
plt.ylim([0, ylim[1]])
plt.xlim([0, xlim[1]])
plt.ylabel('$PDF(x)$')
plt.ylabel('Response, $x$')
plt.grid(True)
plt.ticklabel_format(style='sci', axis='x', scilimits=(0, 0))
plt.ticklabel_format(style='sci', axis='y', scilimits=(0, 0))
plt.legend()

ax = plt.subplot(2, 1, 2)
plt.hold(True)
plt.plot(x_e, peaks_dist.cdf(x_e), 'g-')
plt.plot(x_e, stextreme_dist.cdf(x_e), 'r-')
xlim = ax.get_xlim()
ylim = ax.get_ylim()
plt.ylim([0, ylim[1]])
plt.xlim([0, xlim[1]])
plt.xlabel('Response, $x$')
plt.ylabel('$CDF(x)$')
plt.grid(True)
plt.ticklabel_format(style='sci', axis='x', scilimits=(0, 0))

# goodness of fit plots
gof_plots = ecm.goodnessOfFitPlots(data=peaks, prob_func=peaks_dist, np_return=1000001, x_pdf=x_e, bins_pdf=20)

plt.show()
