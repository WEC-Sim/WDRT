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

import numpy as np
import scipy.stats as stats
import scipy.optimize as optim
import scipy.interpolate as interp
from sklearn.decomposition import PCA


def getContours(Hs, T, depth, size_bin, nb_steps, Time_SS, Time_r):
    '''WDRT Extreme Sea State Contour (ESSC) function

    This function calculates environmental contours of extreme sea states using
    principal component analysis and the inverse first-order reliability
    method.

    Parameters
    ----------
        Hs : np.array
            Vector of Hs values for each measurement in the input.
        T : np.array
            Vector of Te or Tp values for each measurement in the input.
        depth : float
            Water depth at buoy under analysis.
        size_bin : float
            Bin size that will be used to split Component 2 into bins
            after ranking according to values of Component 1.
        nb_steps : float
            Discretization of the circle in the normal space used for
            inverse FORM calculation.
        Time_SS : float
            Sea state duration (hours) of measurements in input.
        Time_r : np.array
            Desired return period (years) for calculation of environmental
            contour, can be a scalar or a vector.

    Returns
    -------
    Hs_Return : np.array
        Calculated Hs values along the contour boundary following
        return to original input orientation.
    T_Return : np.array
       Calculated T values along the contour boundary following
       return to original input orientation.

    Example
    -------
    To obtain the contours for a NDBC buoy::

        import numpy as np
        import WDRT.NDBCdata as NDBCdata
        import WDRT.ESSC as ESSC


        # Pull spectral data from NDBC website
        url = "http://www.ndbc.noaa.gov/station_history.php?station=46022"
        swdList, freqList, dateVals = NDBCdata.fetchFromWeb(46089, savePath='data')

        # Find relevant stats (Hs and Te)
        n = len(swdList)
        Hs = []
        T = []
        DateNum = []
        for ii in range(n):
            tmp1, tmp2 = NDBCdata.getStats(swdList[ii], freqList[ii])
            Hs.extend(tmp1)
            T.extend(tmp2)
            DateNum.extend(NDBCdata.getDateNums(dateVals[ii]))
        Hs = np.array(Hs, dtype=np.float)
        T = np.array(T, dtype=np.float)
        DateNum = np.array(DateNum, dtype=np.float)

        # Removing NaN data, assigning T label depending on input (Te or Tp)
        Nanrem = np.logical_not(np.isnan(T) | np.isnan(Hs))  # Find NaN data in Hs or T
        DateNum = DateNum[Nanrem]  # Remove any NaN data from DateNum
        Hs = Hs[Nanrem]  # Remove any NaN data from Hs
        T = T[Nanrem]  # Remove any NaN data from T

        # Declare required parameters
        depth = 391.4  # Depth at measurement point (m)
        size_bin = 250.  # Enter chosen bin size
        nb_steps = 1000.  # Enter discretization of the circle in the normal space
        # used for inverse FORM calculation
        Time_SS = 1.  # Sea state duration (hrs)
        Time_r = np.array([100])  # Return periods (yrs) of interest
        SteepMax = 0.07  # Optional: enter estimate of breaking steepness

        # Contour generation example
        Hs_Return, T_Return, _, _, _, _, _ = ESSC.getContours(Hs, T, depth, size_bin, nb_steps, Time_SS,
                                                       Time_r)
        # Sample Generation Example
        num_contour_points = 20  # Number of points to be sampled for each
        # contour interval.
        contour_probs = 10 ** (-1 * np.array([1, 2, 2.5, 3, 3.5, 4, 4.5, 5, 5.5, 6]))
        # Probabilities defining sampling contour bounds.
        random_seed = 2  # Random seed for sample generation

        Hs_sample, T_sample, Weight_sample = ESSC.getSamples(Hs, T, num_contour_points, contour_probs, random_seed, depth, size_bin, nb_steps, Time_SS, Time_r)
    '''

    n_data = len(Hs)  # Number of observations

    # Principal Component Analysis (PCA)
    # Perform PCA on significant wave height and energy period after first
    # centering the data.
    pca = PCA(n_components=2)
    pca.fit(np.array((Hs - Hs.mean(axis=0), T - T.mean(axis=0))).T)
    coeff = abs(pca.components_)  # Apply correct/expected sign convention
    coeff[1, 1] = -1.0 * coeff[1, 1]  # Apply correct/expected sign convention
    Comp1_Comp2 = np.dot(np.array((Hs, T)).T, coeff)

    shift = abs(min(Comp1_Comp2[:, 1])) + 0.1  # Calculate shift
    # Apply shift to Component 2 to make all values positive
    Comp1_Comp2[:, 1] = Comp1_Comp2[:, 1] + shift

    Comp1_Comp2_sort = Comp1_Comp2[Comp1_Comp2[:, 0].argsort(), :]

    # Fitting distribution of component 1
    Comp1_params = stats.invgauss.fit(Comp1_Comp2_sort[:, 0], floc=0)

    # Binning Component 2 by Component 1 and
    # application of normal distribution fits for Component 2 in each bin
    # Calculate location of bin edges
    edges = np.hstack((np.arange(0, size_bin * np.ceil(n_data / size_bin),
                                 size_bin), n_data + 1))
    ranks = np.arange(n_data)
    hist_count, _ = np.histogram(ranks, bins=edges)
    bin_inds = np.digitize(ranks, bins=edges) - 1
    Comp2_bins_params = np.zeros((2, int(max(bin_inds) + 1)))
    Comp1_mean = np.array([])

    for bin_loop in range(np.max(bin_inds) + 1):
        mask_bins = bin_inds == bin_loop  # Find location of bin values
        Comp2_bin = np.sort(Comp1_Comp2_sort[mask_bins, 1])
        Comp1_mean = np.append(Comp1_mean,
                               np.mean(Comp1_Comp2_sort[mask_bins, 0]))
        # Calcualte normal distribution parameters for C2 in each bin
        Comp2_bins_params[:, bin_loop] = np.array(stats.norm.fit(Comp2_bin))

    # mu function mu_fcn defined outside of ESSC function

    # Calculate mu parameters from linear fit function

    mu_param, pcov = optim.curve_fit(mu_fcn,
                                     Comp1_mean.T, Comp2_bins_params[0, :])

    # Sigma fitting fuction
    def betafcn(sig_p, rho):
        '''Penalty calculation for sigma parameter fitting function to impose
        positive value constraint.

        Parameters
        ----------
        sig_p: np.array
               Array of sigma fitting function parameters.
        rho: float
             Penalty function variable that drives the solution towards
             required constraint.

        Returns
        -------
        Beta1: float
               Penalty function variable that applies the constraint requiring
               the y-intercept of the sigma fitting function to be greater than
               or equal to 0.
        Beta2: float
               Penalty function variable that applies the constraint requiring
               the minimum of the sigma fitting function to be greater than or
               equal to 0.
        '''
        if -sig_p[2] <= 0:
            Beta1 = 0.0
        else:
            Beta1 = rho
        if -sig_p[2] + (sig_p[1]**2) / (4 * sig_p[0]) <= 0:
            Beta2 = 0.0
        else:
            Beta2 = rho
        return Beta1, Beta2

    # Sigma function sigma_fcn defined outside of ESSC function

    def objfun(sig_p, x, y_actual):
        '''Sum of least square error objective function used in sigma
        minimization.

        Parameters
        ----------
        sig_p: np.array
               Array of sigma fitting function parameters.
        x: np.array
           Array of values (Component 1 mean for each bin) at which to evaluate
           the sigma fitting function.
        y_actual: np.array
                  Array of actual sigma values for each bin to use in least
                  square error calculation with fitted values.

        Returns
        -------
        obj_fun_result: float
                        Sum of least square error objective function for fitted
                        and actual values.
        '''
        obj_fun_result = np.sum((sigma_fcn(sig_p, x) - y_actual)**2)
        return obj_fun_result  # Sum of least square error

    def objfun_penalty(sig_p, x, y_actual, Beta1, Beta2):
        '''Penalty function used for sigma function constrained optimization.

        Parameters
        ----------
        sig_p: np.array
               Array of sigma fitting function parameters.
        x: np.array
           Array of values (Component 1 mean for each bin) at which to evaluate
           the sigma fitting function.
        y_actual: np.array
                  Array of actual sigma values for each bin to use in least
                  square error calculation with fitted values.
        Beta1: float
               Penalty function variable that applies the constraint requiring
               the y-intercept of the sigma fitting function to be greater than
               or equal to 0.
        Beta2: float
               Penalty function variable that applies the constraint requiring
               the minimum of the sigma fitting function to be greater than or
               equal to 0.

        Returns
        -------
        penalty_fcn: float
                     Objective function result with constraint penalties
                     applied for out of bound solutions.
        '''
        penalty_fcn = (objfun(sig_p, x, y_actual) + Beta1 * (-sig_p[2])**2 +
                       Beta2 * (-sig_p[2] + (sig_p[1]**2) / (4 * sig_p[0]))**2)
        return penalty_fcn

    def sigma_fits(Comp1_mean, sigma_vals):
        '''Sigma parameter fitting function using penalty optimization.

        Parameters
        ----------
        Comp1_mean: np.array
                    Mean value of Component 1 for each bin of Component 2.
        sigma_vals: np.array
                    Value of Component 2 sigma for each bin derived from normal
                    distribution fit.

        Returns
        -------
        sig_final: np.array
                   Final sigma parameter values after constrained optimization.
        '''
        sig_0 = np.array((0.1, 0.1, 0.1))  # Set initial guess
        rho = 1.0  # Set initial penalty value
        # Set tolerance, very small values (i.e.,smaller than 10^-5) may cause
        # instabilities
        epsilon = 10**-5
        # Set inital beta values using beta function
        Beta1, Beta2 = betafcn(sig_0, rho)
        # Initial search for minimum value using initial guess
        sig_1 = optim.fmin(func=objfun_penalty, x0=sig_0,
                           args=(Comp1_mean, sigma_vals, Beta1, Beta2), disp=False)
        # While either the difference between iterations or the difference in
        # objective function evaluation is greater than the tolerance, continue
        # iterating
        while (np.amin(abs(sig_1 - sig_0)) > epsilon and
               abs(objfun(sig_1, Comp1_mean, sigma_vals) -
                   objfun(sig_0, Comp1_mean, sigma_vals)) > epsilon):
            sig_0 = sig_1
            # Calculate penalties for this iteration
            Beta1, Beta2 = betafcn(sig_0, rho)
            # Find a new minimum
            sig_1 = optim.fmin(func=objfun_penalty, x0=sig_0,
                               args=(Comp1_mean, sigma_vals, Beta1, Beta2), disp=False)
            rho = 10 * rho  # Increase penalization
        sig_final = sig_1
        return sig_final

    sigma_param = sigma_fits(Comp1_mean, Comp2_bins_params[1, :])

    # IFORM
    # Failure probability for the desired return period (Time_r) given the
    # duration of the measurements (Time_SS)
    p_f = 1 / (365 * (24 / Time_SS) * Time_r)
    beta = stats.norm.ppf((1 - p_f), loc=0, scale=1)  # Reliability
    theta = np.linspace(0, 2 * np.pi, num=nb_steps)
    # Vary U1, U2 along circle sqrt(U1^2+U2^2)=beta
    U1 = beta * np.cos(theta)
    U2 = beta * np.sin(theta)
    # Calculate C1 values along the contour
    Comp1_R = stats.invgauss.ppf(stats.norm.cdf(U1, loc=0, scale=1),
                                 mu=Comp1_params[0], loc=0,
                                 scale=Comp1_params[2])
    # Calculate mu values at each point on the circle
    mu_R = mu_fcn(Comp1_R, mu_param[0], mu_param[1])
    # Calculate sigma values at each point on the circle
    sigma_R = sigma_fcn(sigma_param, Comp1_R)
    # Use calculated mu and sigma values to calculate C2 along the contour
    Comp2_R = stats.norm.ppf(stats.norm.cdf(U2, loc=0, scale=1),
                             loc=mu_R, scale=sigma_R)

    # Re-rotate

    # Principal component rotation formula princomp_inv defined outside of ESSC
    # function

    # Calculate Hs and T along the contour
    Hs_Return, T_Return = princomp_inv(Comp1_R, Comp2_R, coeff, shift)
    Hs_Return = np.maximum(0, Hs_Return)  # Remove negative values
    return Hs_Return, T_Return, coeff, shift, Comp1_params, mu_param, sigma_param


def mu_fcn(x, mu_p_1, mu_p_2):
    ''' Linear fitting function for the mean(mu) of Component 2 normal
    distribution as a function of the Component 1 mean for each bin.
    Used in the ESSC and getSamples functions.

    Parameters
    ----------
    mu_p: np.array
           Array of mu fitting function parameters.
    x: np.array
       Array of values (Component 1 mean for each bin) at which to evaluate
       the mu fitting function.

    Returns
    -------
    mu_fit: np.array
            Array of fitted mu values.
    '''
    mu_fit = mu_p_1 * x + mu_p_2
    return mu_fit


def sigma_fcn(sig_p, x):
    '''Quadratic fitting formula for the standard deviation(sigma) of Component
    2 normal distribution as a function of the Component 1 mean for each bin.
    Used in the ESSC and getSamples functions.

    Parameters
    ----------
    sig_p: np.array
           Array of sigma fitting function parameters.
    x: np.array
       Array of values (Component 1 mean for each bin) at which to evaluate
       the sigma fitting function.

    Returns
    -------
    sigma_fit: np.array
               Array of fitted sigma values.
    '''
    sigma_fit = sig_p[0] * x**2 + sig_p[1] * x + sig_p[2]
    return sigma_fit


def princomp_inv(princip_data1, princip_data2, coeff, shift):
    '''Takes the inverse of the principal component rotation given data,
    coefficients, and shift. Used in the ESSC and getSamples functions.

    Parameters
    ----------
    princip_data1: np.array
                   Array of Component 1 values.
    princip_data2: np.array
                   Array of Component 2 values.
    coeff: np.array
           Array of principal component coefficients.
    shift: float
           Shift applied to Component 2 to make all values positive.

    Returns
    -------
    original1: np.array
               Hs values following rotation from principal component space.
    original2: np.array
               T values following rotation from principal component space.
    '''
    original1 = np.zeros(len(princip_data1))
    original2 = np.zeros(len(princip_data1))
    for i in range(len(princip_data2)):
        original1[i] = (((coeff[0, 1] * (princip_data2[i] - shift)) +
                         (coeff[0, 0] * princip_data1[i])) / (coeff[0, 1]**2 +
                                                              coeff[0, 0]**2))
        original2[i] = (((coeff[0, 1] * princip_data1[i]) -
                         (coeff[0, 0] * (princip_data2[i] -
                                         shift))) / (coeff[0, 1]**2 + coeff[0, 0]**2))
    return original1, original2


def getSamples(Hs, T, num_contour_points, contour_probs, random_seed, depth,
               size_bin, nb_steps, Time_SS, Time_r):
    '''WDRT Extreme Sea State Contour (ESSC) Sampling function

    This function calculates samples of Hs and T using the ESSC function to
    sample between contours of user-defined probabilities.

    Parameters
    ----------
    Hs : np.array
        Vector of Hs values for each measurement in the input.
    T : np.array
        Vector of Te or Tp values for each measurement in the input.
    num_contour_points : int
        Number of sample points to be calculated per contour interval.
    contour_probs: np.array
        Vector of probabilities that define the contour intervals in
        which samples will be taken. Values must be greater than zero and
        less than 1.
    random_seed: int
        Random seed for sample generation, required for sample
        repeatability.
    depth : float
        Water depth at buoy under analysis.
    size_bin : float
        Bin size that will be used to split Component 2 into bins
        after ranking according to values of Component 1.
    nb_steps : float
        Discretization of the circle in the normal space used for
        inverse FORM calculation.
    Time_SS : float
        Sea state duration (hours) of measurements in input.
    Time_r : np.array
        Desired return period (years) for calculation of environmental
        contour, can be a scalar or a vector.

    Returns
    -------
    Hs_samples: np.array
        Vector of Hs values for each sample point.
    Te_samples: np.array
        Vector of Te values for each sample point.
    Weight_points: np.array
        Vector of probabilistic weights for each sampling point
        to be used in risk calculations.

    Example
    -------
    To get weighted samples from a set of contours::

        import numpy as np
        import ESSC

        # Load data from existing text files
        # swdList, freqList, dateVals = NDBCdata.loadFromText(
        #     os.path.join('data', 'NDBC46022'))

        # Find relevant stats (Hs and Te)
        n = len(swdList)
        Hs = []
        T = []
        DateNum = []
        for ii in range(n):
            tmp1, tmp2 = NDBCdata.getStats(swdList[ii], freqList[ii])
            Hs.extend(tmp1)
            T.extend(tmp2)
            DateNum.extend(NDBCdata.getDateNums(dateVals[ii]))
        Hs = np.array(Hs, dtype=np.float)
        T = np.array(T, dtype=np.float)
        DateNum = np.array(DateNum, dtype=np.float)

        # Removing NaN data, assigning T label depending on input (Te or Tp)
        Nanrem = np.logical_not(np.isnan(T) | np.isnan(Hs))  # Find NaN data in Hs or T
        DateNum = DateNum[Nanrem]  # Remove any NaN data from DateNum
        Hs = Hs[Nanrem]  # Remove any NaN data from Hs
        T = T[Nanrem]  # Remove any NaN data from T

        depth = float(675) # Depth at measurement point (m)
        size_bin = float(250) # Enter chosen bin size
        nb_steps = float(1000) # Enter discretization of the circle in the
        # normal space. Used for inverse FORM calculation.
        Time_SS = float(1) # Sea state duration (hrs)
        Time_r = np.array([100]) # Return periods (yrs) of interest

        # Removing NaN data, assigning T label depending on input (Te or Tp)
        Nanrem = np.logical_not(np.isnan(T) | np.isnan(Hs)) # Find NaN data
        DateNum = DateNum[Nanrem] #Remove any NaN data from DateNum
        Hs = Hs[Nanrem] #Remove any NaN data from Hs
        T = T[Nanrem] #Remove any NaN data from T

        num_contour_points = 10 # Number of points to be sampled for each
        # contour interval.
        contour_probs = 10**(-1*np.array([1,2,2.5,3,3.5,4,4.5,5,5.5,6]))
        # Probabilities defining sampling contour bounds.
        random_seed = 2 # Random seed for sample generation

        Hs_sample,T_sample,Weight_points = ESSC.getSamples(Hs,T,
        num_contour_points,contour_probs,random_seed,depth,size_bin,nb_steps,
        Time_SS,Time_r)
    '''
    # Call getContours function to calculate required parameters for sampling
    Hs_Return, T_Return, coeff, shift, Comp1_params, mu_param, sigma_param = getContours(Hs, T, depth, size_bin, nb_steps, Time_SS, Time_r)

    # Calculate line where Hs = 0 to avoid sampling Hs in negative space
    Te_zeroline = np.linspace(2.5, 30, 1000)
    Te_zeroline = np.transpose(Te_zeroline)
    Hs_zeroline = np.zeros(len(Te_zeroline))

    # Transform zero line into principal component space
    Comp_zeroline = np.dot(np.transpose(np.vstack([Hs_zeroline, Te_zeroline])),
                           coeff)
    Comp_zeroline[:, 1] = Comp_zeroline[:, 1] + shift

    # Find quantiles along zero line
    C1_zeroline_prob = stats.invgauss.cdf(Comp_zeroline[:, 0],
                                          mu=Comp1_params[0], loc=0,
                                          scale=Comp1_params[2])
    mu_zeroline = mu_fcn(Comp_zeroline[:, 0], mu_param[0], mu_param[1])
    sigma_zeroline = sigma_fcn(sigma_param, Comp_zeroline[:, 0])
    C2_zeroline_prob = stats.norm.cdf(Comp_zeroline[:, 1],
                                      loc=mu_zeroline, scale=sigma_zeroline)
    C1_normzeroline = stats.norm.ppf(C1_zeroline_prob, 0, 1)
    C2_normzeroline = stats.norm.ppf(C2_zeroline_prob, 0, 1)

    # Reliability contour generation
    beta_lines = stats.norm.ppf(
        (1 - contour_probs), 0, 1)  # Calculate reliability
    beta_lines = np.hstack((0, beta_lines))  # Add zero as lower bound to first
    # contour
    theta_lines = np.linspace(0, 2 * np.pi, 1000)  # Discretize the circle

    contour_probs = np.hstack((1, contour_probs))  # Add probablity of 1 to the
    # reliability set, corresponding to probability of the center point of the
    # normal space

    # Vary U1,U2 along circle sqrt(U1^2+U2^2) = beta
    U1_lines = np.dot(np.cos(theta_lines[:, None]), beta_lines[None, :])
    U2_lines = np.dot(np.sin(theta_lines[:, None]), beta_lines[None, :])

    # Removing values on the H_s = 0 line that are far from the circles in the
    # normal space that will be evaluated to speed up calculations
    minval = np.amin(U1_lines) - 0.5
    mask = C1_normzeroline > minval
    C1_normzeroline = C1_normzeroline[mask]
    C2_normzeroline = C2_normzeroline[mask]

    # Transform to polar coordinates
    Theta_zeroline = np.arctan2(C2_normzeroline, C1_normzeroline)
    Rho_zeroline = np.sqrt(C1_normzeroline**2 + C2_normzeroline**2)
    Theta_zeroline[Theta_zeroline < 0] = Theta_zeroline[
        Theta_zeroline < 0] + 2 * np.pi

    num_samples = (len(beta_lines) - 1) * num_contour_points
    np.random.seed(random_seed)

    Alpha_bounds = np.zeros((len(beta_lines) - 1, 2))
    Angular_dist = np.zeros(len(beta_lines) - 1)
    Angular_ratio = np.zeros(len(beta_lines) - 1)
    Alpha = np.zeros((len(beta_lines) - 1, num_contour_points + 1))
    Weight = np.zeros(len(beta_lines) - 1)
    Sample_beta = np.zeros(num_samples)
    Sample_alpha = np.zeros(num_samples)
    Weight_points = np.zeros(num_samples)

    for i in range(len(beta_lines) - 1):  # Loop over contour intervals
        # Check if any of the radii for the
        r = Rho_zeroline - beta_lines[i + 1]
        # Hs=0, line are smaller than the radii of the contour, meaning
        # that these lines intersect
        if any(r < 0):
            left = np.amin(np.where(r < -0.01))
            right = np.amax(np.where(r < -0.01))
            Alpha_bounds[i, :] = (Theta_zeroline[left], Theta_zeroline[right] -
                                  2 * np.pi)  # Save sampling bounds
        else:
            Alpha_bounds[i, :] = np.array((0, 2 * np.pi))
        # Find the angular distance that will be covered by sampling the disc
        Angular_dist[i] = sum(abs(Alpha_bounds[i]))
        # Calculate ratio of area covered for each contour
        Angular_ratio[i] = Angular_dist[i] / (2 * np.pi)
        # Discretize the remaining portion of the disc into 10 equally spaced
        # areas to be sampled
        Alpha[i, :] = np.arange(min(Alpha_bounds[i]),
                                max(Alpha_bounds[i]) + 0.1, Angular_dist[i] / num_contour_points)
        # Calculate the weight of each point sampled per contour
        Weight[i] = ((contour_probs[i] - contour_probs[i + 1]) *
                     Angular_ratio[i] / num_contour_points)
        for j in range(num_contour_points):
            # Generate sample radius by adding a randomly sampled distance to
            # the 'disc' lower bound
            Sample_beta[(i) * num_contour_points + j] = (beta_lines[i] +
                                                         np.random.random_sample() * (beta_lines[i + 1] - beta_lines[i]))
            # Generate sample angle by adding a randomly sampled distance to
            # the lower bound of the angle defining a discrete portion of the
            # 'disc'
            Sample_alpha[(i) * num_contour_points + j] = (Alpha[i, j] +
                                                          np.random.random_sample() * (Alpha[i, j + 1] - Alpha[i, j]))
            # Save the weight for each sample point
            Weight_points[(i) * num_contour_points + j] = Weight[i]

    # Transform to cartesian coordinates
    Sample_U1 = Sample_beta * np.cos(Sample_alpha)
    Sample_U2 = Sample_beta * np.sin(Sample_alpha)

    # Sample transformation to principal component space
    Comp1_sample = stats.invgauss.ppf(stats.norm.cdf(Sample_U1, loc=0, scale=1),
                                      mu=Comp1_params[0], loc=0,
                                      scale=Comp1_params[2])
    mu_sample = mu_fcn(Comp1_sample, mu_param[0], mu_param[1])
    # Calculate sigma values at each point on the circle
    sigma_sample = sigma_fcn(sigma_param, Comp1_sample)
    # Use calculated mu and sigma values to calculate C2 along the contour
    Comp2_sample = stats.norm.ppf(stats.norm.cdf(Sample_U2, loc=0, scale=1),
                                  loc=mu_sample, scale=sigma_sample)
    # Sample transformation into Hs-T space
    Hs_Sample, T_Sample = princomp_inv(
        Comp1_sample, Comp2_sample, coeff, shift)
    return Hs_Sample, T_Sample, Weight_points


def getContourPoints(T_Return, Hs_Return, T_Sample):
    '''Get points along a specified environmental contour.

    Parameters
    ----------
        T_Return : nparray
            points defining period of return contour
        Hs_Return : nparray
            points defining sig. wave height of return contour
        T_Sample : nparray
            points for sampling along return contour

    Returns
    -------
        Hs_Sample : nparray
            points sampled along return contour
    '''
    amin = np.argmin(T_Return)
    amax = np.argmax(T_Return)

    w1 = Hs_Return[amin:amax]
    w2 = np.concatenate((Hs_Return[amax:], Hs_Return[:amin]))
    if (np.max(w1) > np.max(w2)):
        x1 = T_Return[amin:amax]
        y = Hs_Return[amin:amax]
    else:
        x1 = np.concatenate((T_Return[amax:], T_Return[:amin]))
        y1 = np.concatenate((Hs_Return[amax:], Hs_Return[:amin]))

    ms = np.argsort(x1)
    x = x1[ms]
    y = y1[ms]

    si = interp.interp1d(x, y)
    Hs_Sample = si(T_Sample)
    return Hs_Sample


def steepness(depth, SteepMax, T_vals):
    '''This function calculates a steepness curve to be plotted on an H vs. T
    diagram.  First, the function calculates the wavelength based on the
    depth and T. The T vector can be the input data vector, or will be
    created below to cover the span of possible T values.

    The function solves the dispersion relation for water waves
    using the Newton-Raphson method. All outputs are solved for exactly
    using: (w^2*h/g=kh*tanh(khG)

    Approximations that could be used in place of this code for deep
    and shallow water, as appropriate:
    deep water:h/lambda >= 1/2, tanh(kh)~1, lambda = (g.*T.^2)./(2*.pi)
    shallow water:h/lambda <= 1/20, tanh(kh)~kh, lambda = T.*(g.*h)^0.5

    Parameters
    ----------
    depth: float
        Water depth at site [m].
    SteepMax: float
        Wave breaking steepness estimate (e.g., 0.07).
    T_vals :np.array
        Array of T values [sec] at which to calculate the breaking height.

    Returns
    -------
    SteepH: np.array
        H values [m] that correspond to the T_mesh values creating the
        steepness curve.
    T_steep: np.array
        T values [sec] over which the steepness curve is defined.

    Example
    -------
    To find limit the steepness of waves on a contour by breaking::

        import numpy as np
        import BreakingSteepness

        data_file = np.load('NDBC46022_1996_2012_Hs_Te.npz')
        DateNum = data_file['DateNum']
        T = data_file['T']
        Hs = data_file['Hs']

        depth = 391.1  # Depth at measurement point (m)
        SteepMax = 0.07  # Optional: enter estimate of breaking steepness

        # Removing NaN data, assigning T label depending on input (Te or Tp)
        Nanrem = np.logical_not(np.isnan(T) | np.isnan(Hs)) #Find NaN data in Hs or T
        DateNum = DateNum[Nanrem] #Remove any NaN data from DateNum
        Hs = Hs[Nanrem] #Remove any NaN data from Hs
        T = T[Nanrem] #Remove any NaN data from T

        T_vals = np.arange(0.1,np.amax(T),0.1)

        SteepH = BreakingSteepness.steepness(depth,SteepMax,T_vals)
    '''
    # Calculate the wavelength at a given depth at each value of T
    lambdaT = []

    g = 9.81  # [m/s^2]
    omega = ((2 * np.pi) / T_vals)
    lambdaT = []

    for i in range(len(T_vals)):
        # Initialize kh using Eckert 1952 (mentioned in Holthuijsen pg. 124)
        kh = (omega[i]**2) * depth / \
            (g * (np.tanh((omega[i]**2) * depth / g)**0.5))
        # Find solution using the Newton-Raphson Method
        for j in range(1000):
            kh0 = kh
            f0 = (omega[i]**2) * depth / g - kh0 * np.tanh(kh0)
            df0 = -np.tanh(kh) - kh * (1 - np.tanh(kh)**2)
            kh = -f0 / df0 + kh0
            f = (omega[i]**2) * depth / g - kh * np.tanh(kh)
            if abs(f0 - f) < 10**(-6):
                break

        lambdaT.append((2 * np.pi) / (kh / depth))
        del kh, kh0

    lambdaT = np.array(lambdaT, dtype=np.float)
    SteepH = lambdaT * SteepMax
    return SteepH
