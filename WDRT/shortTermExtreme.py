# Copyright 2016 Sandia Corporation and the National Renewable Energy
# Laboratory
#
# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
#
#     http://www.apache.org/licenses/LICENSE-2.0
#
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.

'''

This module contains tools for analysis of short-term (sea-state specific)
extremes, as well as some usefull statistical methods.

'''

import numpy as np
import scipy.interpolate as interp
import scipy.stats as stats
import scipy.optimize as optim
import matplotlib.pyplot as plt


class ecmDist():
    '''Class for extreme response distributions.

    Class provides basic functionality for statistical distributions as needed
    for extreme value analysis. Can be initialized via either PDF or CDF.

    Note: For pdf, cdf, and ppf toggle ecmDist.{function}.bounds_error to give
    an error or a NaN if input is outside of bounds(x).

    Attributes
    ----------
        pdf : scipy.interpolate.interpolate.interp1d
            Probability density function
        cdf : scipy.interpolate.interpolate.interp1d
            Cumulative distribution function
        ppf : scipy.interpolate.interpolate.interp1d
            Percent point function (inverse cdf)
        getExpVal: method
            Expected value of distribution
        gerRthValue: method
            Rth percentile value of distribution
    '''

    def __init__(self, x, pdf=None, cdf=None):
        '''Class initialization function

        Parameters
        ----------
            x : np.array
                Response variable
            pdf : np.array
                Probability density function for x
            cdf : np.array
                Cumulative distribution function (a.k.a. distribution function)
                for x

        Example
        -------
        To create an ecmDist object, you can use either the PDF::

        >>> myDist1 = ste.ecmDist(x=my_x, pdf=my_PDF)
        or the CDF::

        >>> myDist2 = ste.ecmDist(x=my_x, cdf=my_CDF)
        '''
        # if only CDF is entered
        if pdf is None and not(cdf is None):
            self.cdf = interp.interp1d(x, cdf, bounds_error=True)
            pdf = cdf2pdf(x, cdf)
            self.pdf = interp.interp1d(x, pdf, bounds_error=True)
        # if only PDF is entered
        elif cdf is None and not(pdf is None):
            self.pdf = interp.interp1d(x, pdf, bounds_error=True)
            cdf = pdf2cdf(x, pdf)
            self.cdf = interp.interp1d(x, cdf, bounds_error=True)
        # if both CDF and PDF are entered
        elif any(not(x is None) for x in [pdf, cdf]):
            self.pdf = interp.interp1d(x, pdf, bounds_error=True)
            self.cdf = interp.interp1d(x, cdf, bounds_error=True)
        else:
            ValueError('pdf and/or cdf arguments must not be None')

        self.ppf = interp.interp1d(cdf, x, bounds_error=True)
        self.ccdf = interp.interp1d(x, 1 - cdf, bounds_error=True)

    def getExpVal(self):
        '''Expected value from ecmDist.

        Finds expected value (mean) from distribution.

        Returns
        -------
            ev : float
                Expected value (mean)

        Example
        -------
        In order to find the expected value from the distribution

        >>> vExp = myDist.getExpVal()
        '''
        m0 = np.trapz(self.pdf.y, self.pdf.x)
        m1 = np.trapz(self.pdf.y * self.pdf.x, self.pdf.x)
        ev = m1 / m0
        return ev

    def getRthVal(self, r):
        '''Rth percentile value from ecmDist.

        Finds value of cumulative distribution function for a given percentile.

        Parameters
        ----------
            r : float
                Target percentile (0 < r < 1)

        Returns
        -------
            xr : float
                Response for rth percentile

        Example
        -------
        In order to find the 95th percentile value from the distribution::

        >>> v95 = myDist.getRthVal(r=0.95)
        '''
        m = np.argmin(abs(self.cdf.y - r))
        xr = self.cdf.x[m]
        return xr


def empiricalPdf(x, bins, m=0):
    '''Returns an empirical probability density function.

    Parameters
    ----------
        x : np.array
            Response variable
        bins : int
            Number of bins
        m : float
            Number of points below min(x).
            Use if x consists of only values above a certain minimum and PDF
            for entire population wanted.

    Returns
    -------
        x_pdf : np.array
            Bin (x) locations
        pdf : np.array
            PDF values at x_pdf
    '''
    pdf, y = np.histogram(x, bins=bins, density=True)
    N = len(x)+m
    area = 1. - 1.*m/(1.*N)
    pdf = pdf * area
    x_pdf = []
    for ii in range(len(y) - 1):
        x_pdf.append(y[ii] + 0.5 * (y[ii + 1] - y[ii]))
    x_pdf = np.array(x_pdf)
    return x_pdf, pdf


def empiricalCdf(x,m=0):
    '''Returns an empirical cumulative distribution function.

    Parameters
    ----------
        x : np.array
            Response variable
        m : float
            Number of points below min(x).
            Use if x consists of only values above a certain minimum and CDF
            for entire population wanted.

    Returns
    -------
        x_cdf : np.array
            x locations
        cdf : np.array
            CDF values at x_cdf
    '''
    x_cdf = np.sort(x)
    N = len(x_cdf) + m
    i = np.arange(N) + 1.
    cdf = (i) / (N + 1.)
    return x_cdf, cdf[m:]


def cdf2pdf(x, cdf):
    '''Calculates a probability density function from the cumulative
    distribution function via numerical differentiation.

    Parameters
    ----------
        x : np.array
            Response variable
        cdf : np.array
            Cumulative distribution function values at x

    Returns
    -------
        pdf : np.array
            Probability density function values at x
    '''
    spl = interp.UnivariateSpline(x, cdf)
    spl.set_smoothing_factor(0)
    spl_1 = spl.derivative(n=1)
    pdf = spl_1(x)
    return pdf


def pdf2cdf(x, pdf):
    '''Calculates a cumulative distribution function from the probability
    density function via numerical integration.

    Parameters
    ----------
        x : np.array
            Response variable
        pdf : np.array
            Probability density function values at x

    Returns
    -------
        cdf : np.array
            Cumulative distribution function values at x
    '''
    spl = interp.UnivariateSpline(x, pdf)
    spl.set_smoothing_factor(0)
    spl_1 = spl.derivative(n=-1)
    cdf = spl_1(x)
    return cdf


def loadtxt(filename, delimiter='ws'):
    '''Loads data from a ASCII file.
    First line is header with same delimiter as data, followed by data.

    Parameters
    ----------
        filename : str
            Filename
        delimiter: str
            Data delimiter. Default 'ws' for whitespace, use ',' for CSV files.

    Returns
    -------
        data : dictionary
            Data as a dictionary:
                keys: str
                    Header name for each column.
                values: np.array
                    Column data
    '''
    f = open(filename,'r')
    if delimiter=='ws':
        header = f.readline().rstrip().split()
        data_raw = np.loadtxt(filename, skiprows=1)
    else:
        header = f.readline().rstrip().split(delimiter)
        data_raw = np.loadtxt(filename, delimiter=delimiter, skiprows=1)
    f.close()
    data = {}
    for icol in range(len(header)):
        data[header[icol].strip()] = data_raw[:,icol]
    return data


def globalPeaks(t, data):
    '''Finds the global peaks (maxima between consecutive zero up-crossings)
    of a zero-centered response time-series

    Parameters
    ----------
        t : np.array
            Time vector
        data : np.array
            Response time-series

    Returns
    -------
        t_peaks : np.array
            Time vector for peaks
        peaks : np.array
            Peaks of the response time-series
    '''
    # eliminate zeros
    zeroMask = (data == 0)
    data[zeroMask] = 0.5 * np.min(np.abs(data))
    # zero up-crossings
    diff = np.diff(np.sign(data))
    zeroUpCrossings_mask = (diff == 2) | (diff == 1)
    zeroUpCrossings_index = np.where(zeroUpCrossings_mask)[0]
    zeroUpCrossings_index = np.append(zeroUpCrossings_index, len(data) - 1)
    # global peaks
    N = len(zeroUpCrossings_index)
    peaks = np.array([])
    t_peaks = np.array([])
    for i in range(N - 1):
        peak_index = np.argmax(
            data[zeroUpCrossings_index[i]:zeroUpCrossings_index[i + 1]])
        t_peaks = np.append(t_peaks, t[zeroUpCrossings_index[i] + peak_index])
        peaks = np.append(peaks, data[zeroUpCrossings_index[i] + peak_index])
    # return
    return t_peaks, peaks


def blockMaxima(x, t, t_st):
    '''Finds the block maxima of a time-series.

    Parameters
    ----------
        x : np.array
            Independent random variable (global peaks)
        t : np.array
            Time vector corresponding to x
        t_st : float
            Short-term period

    Returns
    -------
        block_maxima: np.array
            Block maxima (i.e. largest peak in each block)
    '''

    nblock  = int(t[-1] / t_st)
    block_maxima = np.zeros(int(nblock))
    for iblock in range(nblock):
        ix = x[(t >= iblock * t_st) & (t < (iblock+1)*t_st)]
        nx = len(ix)
        block_maxima[iblock] = np.max(ix)
    return block_maxima


def extremeDistribution_Weibull(x, x_e, t_x, t_st, locFlag=0):
    '''Approximates the short-term extreme distribution using the all peaks
    Weibull method.

    Parameters
    ----------
        x : np.array
            Independent random variable (global peaks)
        x_e : np.array
            Array of x values at which to evaluate the short-term extreme CDF
        t_x : float
            Time length of the x array
        t_st : float
            Short-term period
        locFlag : boolean
            locFlag = 0: Location parameter of Weibull distribution is forced to zero
            locFlag = 1: Location parameter of Weibull distribution is calculated in fit


    Returns
    -------
        stextreme_dist: ecmDist object
            Probability distribution of the short-term extreme.
        stextreme_dist: ecmDist object
            Probability distribution of the short-term extreme.
        peaks_dist: scipy.stats rv_frozen
            Probability distribution of the peaks.
        peaks_params: np.array length 4
            Parameters of peak's distribution (Weibull)
            [shape_a, shape_c, loc, scale].
    '''
    # peaks distribution
    if locFlag == 0:
        peaks_params = stats.exponweib.fit(x, f0=1, floc=0)
    elif locFlag == 1:
        peaks_params = stats.exponweib.fit(x, f0=1)
    peaks_dist = stats.exponweib(a=peaks_params[0],
                                 c=peaks_params[1],
                                 loc=peaks_params[2],
                                 scale=peaks_params[3])
    # short-term extreme distribution
    ratio = t_st / t_x
    N = len(x)
    N_st = N * ratio
    weib_cdf = peaks_dist.cdf(x_e)
    ste_cdf = weib_cdf ** N_st
    stextreme_dist = ecmDist(x_e, cdf=ste_cdf)
    # return
    return stextreme_dist, peaks_dist, peaks_params


def extremeDistribution_WeibullTailFit(x, x_e, t_x, t_st, avg=0, p0=None):
    '''Approximates the short-term extreme distribution using the Weibull tail
    fit method.

    Parameters
    ----------
        x : np.array
            Independent random variable (global peaks)
        x_e : np.array
            Array of x values at which to evaluate the short-term extreme CDF
        t_x : float
            Time length of the x array
        t_st : float
            Short-term period
        avg : float
            The average of the time response, if this was substracted before
            identifying global peaks. Else it is assumed that the average is
            zero.
        p0 : list length 2: [float, float]
            Initial guess for the Weibull parameters [shape, scale]

    Returns
    -------
        stextreme_dist: ecmDist object
            Probability distribution of the short-term extreme.
        stextreme_dist : ecmDist object
            Probability distribution of the short-term extreme.
        peaks_dist : scipy.stats rv_frozen
            Probability distribution of the peaks.
        subset_shape_params : np.array length 7
            Shape parameter for each of the seven Weibull fits for the
            subsets of data corresponding to F>[0.65,0.7,...,0.95].
        subset_scale_params : np.array length 7
            Scale parameter for each of the seven Weibull fits for the
            subsets of data corresponding to F>[0.65,0.7,...,0.95].
        peaks_params: np.array length 4
            Parameters of peak's distribution (Weibull)
            [shape_a, shape_c, loc, scale].
    '''
    # Two-parameter weibull distribution def
    def weibCDF(yy, shape, scale):
        loc = 0
        return 1. - np.exp(-1. * ((yy - loc) / scale)**shape)
    # Initial guess for Weibull parameters
    if p0 is None:
        p0_tmp = stats.exponweib.fit(x, f0=1, floc=0)
        p0 = np.zeros(2)
        p0[0] = p0_tmp[1]
        p0[1] = p0_tmp[3]
    # Approximate CDF
    x = np.sort(x)
    N = len(x)
    F = np.zeros(N)
    for i in range(N):
        F[i] = i / (N + 1.0)
    # Divide into seven sets
    subset_shape_params = np.zeros(7)
    subset_scale_params = np.zeros(7)
    setLim = np.arange(0.60, 0.91, 0.05)
    for set in range(7):
        xset = x[(F > setLim[set])]
        Fset = F[(F > setLim[set])]
        popt, _ = optim.curve_fit(weibCDF, xset, Fset, p0=p0)
        subset_shape_params[set] = popt[0]
        subset_scale_params[set] = popt[1]
    # peaks distribution
    peaks_params = [1, np.mean(subset_shape_params), avg,
                    np.mean(subset_scale_params)]
    peaks_dist = stats.exponweib(a=peaks_params[0],
                                 c=peaks_params[1],
                                 loc=peaks_params[2],
                                 scale=peaks_params[3])
    # short-term extreme
    ratio = t_st / t_x
    N_st = N * ratio
    weib_cdf = weibCDF(x_e, peaks_params[1], peaks_params[3])
    ste_cdf = weib_cdf ** N_st
    stextreme_dist = ecmDist(x_e, cdf=ste_cdf)
    return stextreme_dist, peaks_dist, subset_shape_params, \
               subset_scale_params, peaks_params


def extremeDistribution_peaksOverThreshold(x, x_e, t_x, t_st, u):
    '''Approximates the short-term extreme distribution using the peaks over
    threshold method.

    Parameters
    ----------
        x : np.array
            Independent random variable (global peaks)
        x_e : np.array
            Array of x values at which to evaluate the short-term extreme CDF
        t_x : float
            Time length of the x array
        t_st : float
            Short-term period
        u : float
            Threshold below which peaks (x) are ignored

    Returns
    -------
        stextreme_dist: ecmDist object
            Probability distribution of the short-term extreme.
        stextreme_dist : ecmDist object
            Probability distribution of the short-term extreme.
        peaks_dist : ecmDist object
            Probability distribution of the peaks.
        peaksOverThreshold_dist: scipy.stats rv_frozen
            Probaility distribution of the peaks over threshold.
        pot_params: np.array length 3
            Parameters of peak over threshold's distribution using Generalized
            Pareto[shape_c, loc, scale].
    '''
    # peaks over threshold
    pot = np.sort(x)
    pot = pot[(pot > u)] - u
    N = len(x)
    Npot = len(pot)
    # Fit a generalized Pareto
    pot_params = stats.genpareto.fit(pot, floc=0.)
    peaksOverThreshold_dist = stats.genpareto(c=pot_params[0],
                                              loc=pot_params[1],
                                              scale=pot_params[2])
    # peaks
    x_e_pot = x_e[x_e>=u]
    genpareto_cdf = peaksOverThreshold_dist.cdf(x_e_pot-u)
    A = 1. - genpareto_cdf
    k = 1.*Npot / (1.*N)
    peaks_cdf = 1. - (k * A)
    peaks_dist = ecmDist(x_e_pot, cdf=peaks_cdf)
    # short-term extreme
    ratio = t_st / t_x
    N_st = 1.*N * ratio
    ste_cdf = peaks_cdf ** N_st
    stextreme_dist = ecmDist(x_e_pot, cdf=ste_cdf)
    # return
    return stextreme_dist, peaks_dist, peaksOverThreshold_dist, pot_params


def extremeDistribution_blockMaximaGEV(x, t, t_st):
    '''Approximates the short-term extreme distribution using the block maxima
    method and the Generalized Extreme Value distribution.

    Parameters
    ----------
        x : np.array
            Independent random variable (global peaks)
        t : np.array
            Time vector corresponding to x
        t_st : float
            Short-term period

    Returns
    -------
        stextreme_dist: scipy.stats rv_frozen
            Probability distribution of the short-term extreme.
        stextreme_dist : scipy.stats rv_frozen
            Probability distribution of the short-term extreme.
        ste_params: np.array length 3
            Parameters of the short term extreme distribution (Generalized
            Extreme Value) [shape_c, loc, scale].
        block_maxima: np.array
            Block maxima (i.e. largest peak in each block).
    '''
    block_maxima = blockMaxima(x, t, t_st)
    ste_parameters = stats.genextreme.fit(block_maxima)
    stextreme_dist = stats.genextreme(c=ste_parameters[0],
                                      loc=ste_parameters[1],
                                      scale=ste_parameters[2])
    return stextreme_dist, ste_parameters, block_maxima


def extremeDistribution_blockMaximaGumb(x, t, t_st):
    '''Approximates the short-term extreme distribution using the block maxima
    method and the Gumbel (right) distribution.

    Parameters
    ----------
        x : np.array
            Independent random variable (global peaks)
        t : np.array
            Time vector corresponding to x
        t_st : float
            Short-term period

    Returns
    -------
        stextreme_dist: scipy.stats rv_frozen
            Probability distribution of the short-term extreme.
        stextreme_dist : scipy.stats rv_frozen
            Probability distribution of the short-term extreme.
        ste_params: np.array length 2
            Parameters of the short term extreme distribution (Gumbel_r)
            [loc, scale].
        block_maxima: np.array
            Block maxima (i.e. largest peak in each block).
    '''
    block_maxima = blockMaxima(x, t, t_st)
    ste_parameters = stats.gumbel_r.fit(block_maxima)
    stextreme_dist = stats.gumbel_r(loc=ste_parameters[0],
                                    scale=ste_parameters[1])

    return stextreme_dist, ste_parameters, block_maxima


def goodnessOfFitPlots(data, prob_func, x_pdf, bins_pdf=20, np_return=100001, m_prob=0., response_name='Response', response_name_2='Peaks', response_units='Response Units'):
    '''Creates plots showing the goodness of fit of a probability model to the
    actual data. The four subplots are: the probability plot, the quantile plot,
    the return level plot, and probability density function (PDF) plot.

    Parameters
    ----------
        data : np.array
            Data which the probability model represents
        prob_func : {scipy.stats rv_frozen} or {ecmDist object}
            Probability function that models the data.
        x_pdf : np.array
            Array of x values at which to evaluate the PDF.
        bins_pdf : int
            Number of bins to use for the PDF plot.
        np_return: int
            Number of points to use in the return level plot.
        m_prob: float
            Number of points below min(data).
            Use if data consists of only values above a certain minimum and
            goodness of fit plots for entire population wanted.

    Returns
    -------
        fgof : matplotlib.pyplot figure
            Figure containing the four goodness of fit subplots.
    '''
    # nan instead of error if out-of-range
    if str(prob_func.__class__)=='WDRT.shortTermExtreme.ecmDist':
        prob_func.pdf.bounds_error=False
        prob_func.cdf.bounds_error=False
        prob_func.ppf.bounds_error=False
    # figure
    fgof = plt.figure()
    ax1 = fgof.add_subplot(2,2,1)
    plt.title('Probability Plot')
    plt.xlabel('Fit CDF')
    plt.ylabel('Empirical CDF')
    plt.ticklabel_format(style='sci', axis='y', scilimits=(0, 0))
    plt.ticklabel_format(style='sci', axis='x', scilimits=(0, 0))
    ax2 = fgof.add_subplot(2,2,2)
    plt.title('Quantile Plot')
    plt.xlabel('Fit Inverse CDF')
    ylabel = response_name + ' ' + response_name_2 + ' [' + response_units + ']'
    plt.ylabel(ylabel)
    plt.ticklabel_format(style='sci', axis='y', scilimits=(0, 0))
    plt.ticklabel_format(style='sci', axis='x', scilimits=(0, 0))
    ax3 = fgof.add_subplot(2,2,3)
    ax3.set_xscale('log')
    plt.title('Return Level Plot')
    plt.xlabel('Return Period [' + response_name_2 + ']')
    ylabel = 'Return level [' + response_units + ']'
    plt.ylabel(ylabel)
    plt.ticklabel_format(style='sci', axis='y', scilimits=(0, 0))
    # plt.ticklabel_format(style='sci', axis='x', scilimits=(0, 0))
    ax4 = fgof.add_subplot(2,2,4)
    plt.title('PDF')
    xlabel = response_name + ' ' + response_name_2 + ' [' + response_units + ']'
    plt.xlabel(xlabel)
    plt.ticklabel_format(style='sci', axis='y', scilimits=(0, 0))
    plt.ticklabel_format(style='sci', axis='x', scilimits=(0, 0))
    plt.ylabel('PDF')
    plt.tight_layout()
    # probability plot
    data = np.sort(data)
    N = len(data)
    _,cdf_emp = empiricalCdf(data,m=m_prob)
    cdf_model = prob_func.cdf(data)
    cdf_emp_p = cdf_emp[~np.isnan(cdf_model)]
    cdf_model = cdf_model[~np.isnan(cdf_model)]
    unit = np.arange(2)
    ax1.plot(cdf_emp_p, cdf_model,'o')
    ax1.plot(unit,unit,'k--')
    # quantile plot
    inv_cdf_model = prob_func.ppf(cdf_emp)
    data_q = data[~np.isnan(inv_cdf_model)]
    inv_cdf_model = inv_cdf_model[~np.isnan(inv_cdf_model)]
    ax2.plot(inv_cdf_model,data_q,'o')
    unit2 = np.arange(int(np.max(data)+1))
    ax2.plot(unit2,unit2,'k--')
    # return level plot
    p = np.linspace(0,1,np_return)
    p = p[1:-1]
    p1 = 1. - p
    xp = 1./p
    zp = prob_func.ppf(p1)
    ax3.plot(xp,zp,'k--')
    ax3.plot(1./(1.-cdf_emp),data,'o')
    # pdf plot
    xpdf_emp,pdf_emp = empiricalPdf(data,bins_pdf,m=m_prob)
    pdf_model = prob_func.pdf(x_pdf)
    x_pdf = x_pdf[~np.isnan(pdf_model)]
    pdf_model = pdf_model[~np.isnan(pdf_model)]
    plt.bar(xpdf_emp,pdf_emp, align='center',
        width=0.75*((data[-1]-data[0])/bins_pdf), color='k',alpha=0.2)
    plt.plot(data,np.zeros(N),'bo')
    plt.plot(x_pdf, pdf_model, 'r-')
    xlim = ax4.get_xlim()
    ylim = ax4.get_ylim()
    plt.xlim([0,xlim[1]])
    plt.ylim([0,ylim[1]])
    # revert to error instead of nan if out-of-range
    if str(prob_func.__class__)=='WDRT.shortTermExtreme.ecmDist':
        prob_func.pdf.bounds_error=True
        prob_func.cdf.bounds_error=True
        prob_func.ppf.bounds_error=True
    # return
    return fgof


def compare_methods(x, t, t_st, methods=[1,2,3,4,5],colors=['g','b','r','k','k'],lines=['-','-','-','-','--']):
    '''Compares the results obtained using different methods to approximate the
    short-term extreme distribution. The methods are:
    1 - All peaks Weibull,
    2 - Weibull tail fit,
    3 - Peaks over threshold,
    4 - Block maxima GEV,
    5 - Block maxima Gumbel,

    Parameters
    ----------
        x : np.array
            Independent random variable (global peaks)
        t : np.array
            Time vector corresponding to x
        t_st : float
            Short-term period
        methods : list
            List of methods to be compared. Options are any combination of: 1, 2, 3, 4,and 5
        colors : list
            Strings defining the color used to plot each method
        lines : list
            Strings defining the lines used to plot each method

    Returns
    -------
        fig1 : matplotlib.pyplot figure
            Figure containing the comparison of peaks distribution.
        fig2 : matplotlib.pyplot figure
            Figure containing the comparison of short-term extreme distribution
        expected_value_of_short_term_extreme : dictionary
            Expected value of the short-term extreme distribution from each method
    '''
    # get the 1-hour extreme distribution using the different methods
    x_e = np.linspace(0, 2 * np.max(x), 10000)
    t_x = (t[-1]-t[0]) + ((t[-1]-t[0])/(1.*len(x)))
    expected_value_of_short_term_extreme = {}
    # 1 - All peaks Weibull
    if 1 in methods:
        m1 = {}
        m1['stextreme_dist'], m1['peaks_dist'], _ = extremeDistribution_Weibull(x=x, x_e=x_e, t_x=t_x, t_st=t_st)
        m1['ev'] = m1['stextreme_dist'].getExpVal()
        expected_value_of_short_term_extreme['all_peaks_weibull'] = m1['ev']
    # 2 - Weibull tail fit
    if 2 in methods:
        m2 = {}
        m2['stextreme_dist'], m2['peaks_dist'], _, _, _ = extremeDistribution_WeibullTailFit(x=x, x_e=x_e, t_x=t_x, t_st=t_st)
        m2['ev'] = m2['stextreme_dist'].getExpVal()
        expected_value_of_short_term_extreme['weibull_tail_fit'] = m2['ev']
    # 3 - Peaks over threshold
    if 3 in methods:
        m3 = {}
        thresh = np.mean(x) + 1.4*np.std(x)
        thresh_x = np.min(x_e[x_e>thresh])
        m3['stextreme_dist'], m3['peaks_dist'], m3['pot_dist'], _ = extremeDistribution_peaksOverThreshold(x=x, x_e=x_e, t_x=t_x, t_st=t_st, u=thresh)
        m3['ev'] = m3['stextreme_dist'].getExpVal()
        expected_value_of_short_term_extreme['peaks_over_threshhold'] = m3['ev']
    # 4 - Block maxima GEV
    if 4 in methods:
        m4 = {}
        m4['stextreme_dist'],_,_ = extremeDistribution_blockMaximaGEV(x=x, t=t, t_st=t_st)
        m4['ev'] = m4['stextreme_dist'].mean()
        expected_value_of_short_term_extreme['block_maxima_gev'] = m4['ev']
    # 5 - Block maxima Gumbel
    if 5 in methods:
        m5 = {}
        m5['stextreme_dist'],_,_ = extremeDistribution_blockMaximaGumb(x=x, t=t, t_st=t_st)
        m5['ev'] = m5['stextreme_dist'].mean()
        expected_value_of_short_term_extreme['block_maxima_gumbel'] = m5['ev']
    # plot peaks distribution
    fig1 = plt.figure()
    ax = plt.subplot(2, 1, 1)
    plt.hold(True)
    if 1 in methods:
        plt.plot(x_e, m1['peaks_dist'].pdf(x_e), colors[0]+lines[0], label='All Peaks Weibull')
    if 2 in methods:
        plt.plot(x_e, m2['peaks_dist'].pdf(x_e), colors[1]+lines[1], label='Weibull Tail Fit')
    if 3 in methods:
        plt.plot(x_e[x_e>thresh_x], m3['peaks_dist'].pdf(x_e[x_e>thresh_x]), colors[2]+lines[2], label='Peaks Over Threshhold')
    xlim = ax.get_xlim()
    ylim = ax.get_ylim()
    plt.ylim([0,ylim[1]])
    plt.xlim([0,xlim[1]])
    plt.ylabel('$PDF(x)$')
    plt.grid(True)
    plt.ticklabel_format(style='sci', axis='x', scilimits=(0, 0))
    plt.ticklabel_format(style='sci', axis='y', scilimits=(0, 0))
    plt.title('Peaks Distribution')
    plt.legend()
    ax = plt.subplot(2, 1, 2)
    plt.hold(True)
    if 1 in methods:
        plt.plot(x_e, m1['peaks_dist'].cdf(x_e), colors[0]+lines[0], label='All Peaks Weibull')
    if 2 in methods:
        plt.plot(x_e, m2['peaks_dist'].cdf(x_e), colors[1]+lines[1], label='Weibull Tail Fit')
    if 3 in methods:
        plt.plot(x_e[x_e>thresh_x], m3['peaks_dist'].cdf(x_e[x_e>thresh_x]), colors[2]+lines[2], label='Peaks Over Threshhold')
    xlim = ax.get_xlim()
    ylim = ax.get_ylim()
    plt.ylim([0,ylim[1]])
    plt.xlim([0,xlim[1]])
    plt.xlabel('Response, $x$')
    plt.ylabel('$CDF(x)$')
    plt.grid(True)
    plt.ticklabel_format(style='sci', axis='x', scilimits=(0, 0))
    plt.ticklabel_format(style='sci', axis='y', scilimits=(0, 0))
    # plot short-term extreme distribution
    fig2 = plt.figure()
    ax = plt.subplot(2, 1, 1)
    plt.hold(True)
    if 1 in methods:
        plt.plot(x_e, m1['stextreme_dist'].pdf(x_e), colors[0]+lines[0], label='All Peaks Weibull')
        plt.plot(m1['ev'],m1['stextreme_dist'].pdf(m1['ev']),colors[0]+'o')
    if 2 in methods:
        plt.plot(x_e, m2['stextreme_dist'].pdf(x_e), colors[1]+lines[1], label='Weibull Tail Fit')
        plt.plot(m2['ev'],m2['stextreme_dist'].pdf(m2['ev']),colors[1]+'o')
    if 3 in methods:
        plt.plot(x_e[x_e>thresh_x], m3['stextreme_dist'].pdf(x_e[x_e>thresh_x]), colors[2]+lines[2], label='Peaks Over Threshhold')
        plt.plot(m3['ev'],m3['stextreme_dist'].pdf(m3['ev']),colors[2]+'o')
    if 4 in methods:
        plt.plot(x_e, m4['stextreme_dist'].pdf(x_e), colors[3]+lines[3], label='Block Maxima (GEV)')
        plt.plot(m4['ev'],m4['stextreme_dist'].pdf(m4['ev']),colors[3]+'o')
    if 5 in methods:
        plt.plot(x_e, m5['stextreme_dist'].pdf(x_e), colors[4]+lines[4], label='Block Maxima (Gumb)')
        plt.plot(m5['ev'],m5['stextreme_dist'].pdf(m5['ev']),colors[4]+'o')
    xlim = ax.get_xlim()
    ylim = ax.get_ylim()
    plt.ylim([0,ylim[1]])
    plt.xlim([0,xlim[1]])
    plt.ylabel('$PDF(x)$')
    plt.grid(True)
    plt.ticklabel_format(style='sci', axis='x', scilimits=(0, 0))
    plt.ticklabel_format(style='sci', axis='y', scilimits=(0, 0))
    plt.title('Short-Term Extreme Distribution')
    plt.legend()
    ax = plt.subplot(2, 1, 2)
    plt.hold(True)
    if 1 in methods:
        plt.plot(x_e, m1['stextreme_dist'].cdf(x_e), colors[0]+lines[0], label='All Peaks Weibull')
        plt.plot(m1['ev'],m1['stextreme_dist'].cdf(m1['ev']),colors[0]+'o')
    if 2 in methods:
        plt.plot(x_e, m2['stextreme_dist'].cdf(x_e), colors[1]+lines[1], label='Weibull Tail Fit')
        plt.plot(m2['ev'],m2['stextreme_dist'].cdf(m2['ev']),colors[1]+'o')
    if 3 in methods:
        plt.plot(x_e[x_e>thresh_x], m3['stextreme_dist'].cdf(x_e[x_e>thresh_x]), colors[2]+lines[2], label='Peaks Over Threshhold')
        plt.plot(m3['ev'],m3['stextreme_dist'].cdf(m3['ev']),colors[2]+'o')
    if 4 in methods:
        plt.plot(x_e, m4['stextreme_dist'].cdf(x_e), colors[3]+lines[3], label='Block Maxima (GEV)')
        plt.plot(m4['ev'],m4['stextreme_dist'].cdf(m4['ev']),colors[3]+'o')
    if 5 in methods:
        plt.plot(x_e, m5['stextreme_dist'].cdf(x_e), colors[4]+lines[4], label='Block Maxima (Gumb)')
        plt.plot(m5['ev'],m5['stextreme_dist'].cdf(m5['ev']),colors[4]+'o')
    xlim = ax.get_xlim()
    ylim = ax.get_ylim()
    plt.ylim([0,ylim[1]])
    plt.xlim([0,xlim[1]])
    plt.xlabel('Response, $x$')
    plt.ylabel('$CDF(x)$')
    plt.grid(True)
    plt.ticklabel_format(style='sci', axis='x', scilimits=(0, 0))
    plt.ticklabel_format(style='sci', axis='y', scilimits=(0, 0))
    # return
    print 'expected_value_of_short_term_extreme:'
    print expected_value_of_short_term_extreme
    return fig1, fig2, expected_value_of_short_term_extreme

