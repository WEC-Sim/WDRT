import numpy as np
import matplotlib.pyplot as plt
import WDRT.shortTermExtreme as ecm
import WDRT.fatigue as fatigue
import os

method = 1
	# 1 - All peaks Weibull
	# 2 - Weibull tail fit
	# 3 - Peaks over threshold
	# 4 - Block maxima GEV
	# 5 - Block maxima Gumbel

# load global peaks
t_peaks_file = os.path.join('data', 't.dat')
peaks_file   = os.path.join('data', 'peaks.dat')
t_peaks = np.loadtxt(t_peaks_file)
peaks = np.loadtxt(peaks_file)/1000.

# get the 1-hour extreme distribution using the method selected above
x_e = np.linspace(0, 2 * np.max(peaks), 10000)
t_x = (t_peaks[-1]-t_peaks[0]) + ((t_peaks[-1]-t_peaks[0])/(1.*len(peaks)))
t_st = 1. * 60. * 60.
if method==1:
	stextreme_dist, peaks_dist, _ = ecm.extremeDistribution_Weibull(x=peaks, x_e=x_e, t_x=t_x, t_st=t_st)
elif method==2:
	stextreme_dist, peaks_dist, _, _, _ = ecm.extremeDistribution_WeibullTailFit(x=peaks, x_e=x_e, t_x=t_x, t_st=t_st)
elif method==3:
	thresh = np.mean(peaks) + 1.4*np.std(peaks)
	thresh_x = np.min(x_e[x_e>thresh])
	stextreme_dist, peaks_dist, pot_dist, _ = ecm.extremeDistribution_peaksOverThreshold(x=peaks, x_e=x_e, t_x=t_x, t_st=t_st, u=thresh)
elif method==4:
	stextreme_dist,_,bm = ecm.extremeDistribution_blockMaximaGEV(x=peaks, t=t_peaks, t_st=t_st)
elif method == 5:
	stextreme_dist,_,bm = ecm.extremeDistribution_blockMaximaGumb(x=peaks, t=t_peaks, t_st=t_st)

# goodness of fit plots
if method==1 or method==2:
	bm = ecm.blockMaxima(x=peaks, t=t_peaks, t_st=t_st)
	_ = ecm.goodnessOfFitPlots(data=peaks, prob_func=peaks_dist, np_return=1000001, x_pdf=x_e, bins_pdf=20, response_name='PTO Force', response_name_2='Peaks',response_units='kN')
if not method==3:
	fig_gof = ecm.goodnessOfFitPlots(data=bm, prob_func=stextreme_dist, np_return=10001, x_pdf=x_e, bins_pdf=20, response_name='PTO Force', response_name_2='1-hr Extreme',response_units='kN')
if method==3:
	bm = ecm.blockMaxima(x=peaks, t=t_peaks, t_st=t_st)
	_ = ecm.goodnessOfFitPlots(data=peaks[peaks>thresh_x], prob_func=peaks_dist, np_return=100001, x_pdf=x_e[x_e>thresh_x], bins_pdf=20,m_prob=1.*len(peaks[peaks<thresh_x]), response_name='PTO Force', response_name_2='Peaks',response_units='kN')
	_ = ecm.goodnessOfFitPlots(data=peaks[peaks>thresh]-thresh, prob_func=pot_dist, np_return=100001, x_pdf=x_e[x_e>thresh]-thresh, bins_pdf=20, response_name='PTO Force', response_name_2='Peaks Over Threshold',response_units='kN')
	fig_gof = ecm.goodnessOfFitPlots(data=bm, prob_func=stextreme_dist, np_return=10001, x_pdf=x_e[x_e>thresh_x], bins_pdf=20, response_name='PTO Force', response_name_2='1-hr Extreme',response_units='kN')

# plot
plt.figure()
if method==3:
	plt.plot(t_peaks[peaks<thresh], peaks[peaks<thresh], 'ko', alpha=0.2)
	plt.plot(t_peaks[peaks>thresh], peaks[peaks>thresh], 'go')
	plt.plot([0, t_peaks[-1]], [thresh, thresh], 'r--')
else:
	plt.plot(t_peaks, peaks, 'go')
plt.plot([0, t_peaks[-1]], [0, 0], 'k--')
plt.xlabel('Time, $t$ [s]')
plt.ylabel('Response, $x$')
plt.xlim([0,3600*2])
plt.grid(True)
plt.ticklabel_format(style='sci', axis='y', scilimits=(0, 0))

plt.figure()
ax = plt.subplot(2, 1, 1)
if method==1 or method==2:
	plt.plot(x_e, peaks_dist.pdf(x_e), 'g-', label='Peak distribution')
if not method==3:
	plt.plot(x_e, stextreme_dist.pdf(x_e), 'r-', label='Extreme distribution')
if method==3:
	plt.plot(x_e[x_e>thresh_x], peaks_dist.pdf(x_e[x_e>thresh_x]), 'g-', label='Peak distribution')
	plt.plot(x_e[x_e>thresh_x], stextreme_dist.pdf(x_e[x_e>thresh_x]), 'r-', label='Extreme distribution')
xlim = ax.get_xlim()
ylim = ax.get_ylim()
if method==3:
	plt.plot([thresh, thresh], [0, ylim[1]], 'k--')
plt.ylim([0,ylim[1]])
plt.xlim([0,xlim[1]])
plt.xlabel('Response, $x$')
plt.ylabel('$PDF(x)$')
plt.grid(True)
plt.ticklabel_format(style='sci', axis='x', scilimits=(0, 0))
plt.ticklabel_format(style='sci', axis='y', scilimits=(0, 0))
plt.legend()

ax = plt.subplot(2, 1, 2)
if method==1 or method==2:
	plt.plot(x_e, peaks_dist.cdf(x_e), 'g-')
if not method==3:
	plt.plot(x_e, stextreme_dist.cdf(x_e), 'r-')
if method==3:
	plt.plot(x_e[x_e>thresh_x], peaks_dist.cdf(x_e[x_e>thresh_x]), 'g-')
	plt.plot(x_e[x_e>thresh_x], stextreme_dist.cdf(x_e[x_e>thresh_x]), 'r-')
xlim = ax.get_xlim()
ylim = ax.get_ylim()
if method==3:
	plt.plot([thresh, thresh], [0, ylim[1]], 'k--')
plt.ylim([0,ylim[1]])
plt.xlim([0,xlim[1]])
plt.xlabel('Response, $x$')
plt.ylabel('$CDF(x)$')
plt.grid(True)
plt.ticklabel_format(style='sci', axis='x', scilimits=(0, 0))

plt.show()
