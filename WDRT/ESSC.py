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
The Extreme Sea State Contour (ESSC) module contains the tools necessary to 
calculate environmental contours of extreme sea states for buoy data. 
'''

import numpy as np
import scipy.stats as stats
import scipy.optimize as optim
import scipy.interpolate as interp
import matplotlib.pyplot as plt
import h5py
from sklearn.decomposition import PCA as skPCA
import requests
import bs4
import urllib.request
import re
from datetime import datetime, date
import os
import glob
import copy
import statsmodels.api as sm
from statsmodels import robust
import urllib
import matplotlib


class EA:
    '''The Environmental Assessment (EA) class points to functions for 
    various contour methods (including getContours and getSamples) and allows 
    the user to plot results (plotData), sample along the contour 
    (getContourPoints), calculate the wave breaking steepness curve (steepness)
    and/or use the bootstrap method to calculate 95% confidence bounds about  
    the contours (bootStrap).'''
    def __init__():
        return
    def getContours():
        '''Points to the getContours function in whatever contouring method is used'''
        return
    def getSamples():
        '''Points to the getSamples function in whatever contouring method is
        used, currently only implemented for PCA contours. Implementation for 
        additional contour methods planned for future release.'''
        return

    def saveContour(self, fileName=None):
        '''Saves all available contour data obtained via the EA module to
        a .h5 file

        Parameters
        ----------
        fileName : string
            relevent path and filename where the .h5 file will be created and
            saved. If no filename, the h5 file will be named NDBC(buoyNum).h5
        '''
        if (fileName is None):
            fileName = 'NDBC' + str(self.buoy.buoyNum) + '.h5'
        else:
            _, file_extension = os.path.splitext(fileName)
            if not file_extension:
                fileName = fileName + '.h5'

        print(fileName);
        with h5py.File(fileName, 'a') as f:
            if('method' in f):
                f['method'][...] = self.method
            else:
                f.create_dataset('method', data=self.method)
            
            if('parameters' in f):
                gp = f['parameters']
            else:
                gp = f.create_group('parameters')
            
            self._saveParams(gp)

            if(self.Hs_ReturnContours is not None):
                if('ReturnContours' in f):                
                    grc = f['ReturnContours']
                else:                   
                    grc = f.create_group('ReturnContours')
                
                if('T_Return' in grc):
                    f_T_Return = grc['T_Return']
                    f_T_Return[...] = self.T_ReturnContours
                else:
                    f_T_Return = grc.create_dataset('T_Return', data=self.T_ReturnContours)                
                
                f_T_Return.attrs['units'] = 's'
                f_T_Return.attrs['description'] = 'contour, energy period'
                
                if('Hs_Return' in grc):
                    f_Hs_Return = grc['Hs_Return']
                    f_Hs_Return[...] = self.Hs_ReturnContours
                else:
                    f_Hs_Return = grc.create_dataset('Hs_Return', data=self.Hs_ReturnContours)
                
                f_Hs_Return.attrs['units'] = 'm'
                f_Hs_Return.attrs['description'] = 'contours, significant wave height'

            # Samples for full sea state long term analysis
            if(hasattr(self, 'Hs_SampleFSS') and self.Hs_SampleFSS is not None):
                if('Samples_FullSeaState' in f):
                    gfss = f['Samples_FullSeaState']
                else:    
                    gfss = f.create_group('Samples_FullSeaState')
                
                if('Hs_SampleFSS' in gfss):
                    f_Hs_SampleFSS = gfss['Hs_SampleFSS']
                    f_Hs_SampleFSS[...] = self.Hs_SampleFSS
                else:
                    f_Hs_SampleFSS = gfss.create_dataset('Hs_SampleFSS', data=self.Hs_SampleFSS)
                
                f_Hs_SampleFSS.attrs['units'] = 'm'
                f_Hs_SampleFSS.attrs['description'] = 'full sea state significant wave height samples'
                
                if('T_SampleFSS' in gfss):
                    f_T_SampleFSS = gfss['T_SampleFSS']
                    f_T_SampleFSS[...] = self.T_SampleFSS
                else:   
                    f_T_SampleFSS = gfss.create_dataset('T_SampleFSS', data=self.T_SampleFSS)                
                
                f_T_SampleFSS.attrs['units'] = 's'
                f_T_SampleFSS.attrs['description'] = 'full sea state energy period samples'                
                
                if('Weight_SampleFSS' in gfss):
                    f_Weight_SampleFSS = gfss['Weight_SampleFSS']
                    f_Weight_SampleFSS[...] = self.Weight_SampleFSS
                else:
                    f_Weight_SampleFSS = gfss.create_dataset('Weight_SampleFSS', data = self.Weight_SampleFSS)
                
                f_Weight_SampleFSS.attrs['description'] = 'full sea state relative weighting samples'

            # Samples for contour approach long term analysis
            if(hasattr(self, 'Hs_SampleCA') and self.Hs_SampleCA is not None):
                if('Samples_ContourApproach' in f):
                    gca = f['Samples_ContourApproach']
                else:
                    gca = f.create_group('Samples_ContourApproach')
                
                if('Hs_SampleCA' in gca):
                    f_Hs_sampleCA = gca['Hs_SampleCA']
                    f_Hs_sampleCA[...] = self.Hs_SampleCA
                else:
                    f_Hs_sampleCA = gca.create_dataset('Hs_SampleCA', data=self.Hs_SampleCA)
                
                f_Hs_sampleCA.attrs['units'] = 'm'
                f_Hs_sampleCA.attrs['description'] = 'contour approach significant wave height samples'
                
                if('T_SampleCA' in gca):
                    f_T_sampleCA = gca['T_SampleCA']
                    f_T_sampleCA[...] = self.T_SampleCA
                else:
                    f_T_sampleCA = gca.create_dataset('T_SampleCA', data=self.T_SampleCA)
                
                f_T_sampleCA.attrs['units'] = 's'
                f_T_sampleCA.attrs['description'] = 'contour approach energy period samples'

    def plotData(self):
        """
        Display a plot of the 100-year return contour, full sea state samples
        and contour samples
        """
        plt.figure()
        plt.plot(self.buoy.T, self.buoy.Hs, 'bo', alpha=0.1, label='NDBC data')
        plt.plot(self.T_ReturnContours, self.Hs_ReturnContours, 'k-', label='100 year contour')
        #plt.plot(self.T_SampleFSS, self.Hs_SampleFSS, 'ro', label='full sea state samples')
        #plt.plot(self.T_SampleCA, self.Hs_SampleCA, 'y^', label='contour approach samples')
        plt.legend(loc='lower right', fontsize='small')
        plt.grid(True)
        plt.xlabel('Energy period, $T_e$ [s]')
        plt.ylabel('Sig. wave height, $H_s$ [m]')
        plt.show()

    def getContourPoints(self, T_Sample):
        '''Get Hs points along a specified environmental contour using 
        user-defined T values.

        Parameters
        ----------
            T_Sample : nparray
                points for sampling along return contour

        Returns
        -------
            Hs_SampleCA : nparray
                points sampled along return contour
        
        Example
        -------
            To calculate Hs values along the contour at specific 
            user-defined T values:
            
                import WDRT.ESSC as ESSC
                import numpy as np
                
                # Pull spectral data from NDBC website
                buoy46022 = ESSC.Buoy('46022','NDBC')
                buoy46022.fetchFromWeb()
                
                # Create PCA EA object for buoy
                pca46022 = ESSC.PCA(buoy46022)
                
                # Declare required parameters
                Time_SS = 1.  # Sea state duration (hrs)
                Time_r = 100  # Return periods (yrs) of interest
                nb_steps = 1000  # Enter discretization of the circle in the normal space (optional)
                
                # Generate contour
                Hs_Return, T_Return = pca46022.getContours(Time_SS, Time_r,nb_steps)
                
                # Use getContourPoints to find specific points along the contour
                T_sampleCA = np.arange(12, 26, 2)
                Hs_sampleCA = pca46022.getContourPoints(T_sampleCA)
        '''
        #finds minimum and maximum energy period values
        amin = np.argmin(self.T_ReturnContours)
        amax = np.argmax(self.T_ReturnContours)
        #finds points along the contour
        w1 = self.Hs_ReturnContours[amin:amax]
        w2 = np.concatenate((self.Hs_ReturnContours[amax:], self.Hs_ReturnContours[:amin]))
        if (np.max(w1) > np.max(w2)):
            x1 = self.T_ReturnContours[amin:amax]
            y = self.Hs_ReturnContours[amin:amax]
        else:
            x1 = np.concatenate((self.T_ReturnContours[amax:], self.T_ReturnContours[:amin]))
            y1 = np.concatenate((self.Hs_ReturnContours[amax:], self.Hs_ReturnContours[:amin]))
        #sorts data based on the max and min energy period values
        ms = np.argsort(x1)
        x = x1[ms]
        y = y1[ms]
        #interpolates the sorted data
        si = interp.interp1d(x, y)
        #finds the wave height based on the user specified energy period values
        Hs_SampleCA = si(T_Sample)

        self.T_SampleCA = T_Sample
        self.Hs_SampleCA = Hs_SampleCA
        return Hs_SampleCA

    def steepness(self, SteepMax, T_vals, depth = None):
        '''This function calculates a steepness curve to be plotted on an H vs. T
        diagram.  First, the function calculates the wavelength based on the
        depth and T. The T vector can be the input data vector, or will be
        created below to cover the span of possible T values.
        The function solves the dispersion relation for water waves
        using the Newton-Raphson method. All outputs are solved for exactly
        using: :math:`hw^2/g = kh*tanh(khG)` 
        
        Approximations that could be used in place of this code for deep
        and shallow water, as appropriate:
            
        deep water: :math:`h/\lambda \geq 1/2, tanh(kh) \sim 1, \lambda = (gT^2)/(2\pi)`
        
        shallow water: :math:`h/\lambda \leq 1/20, tanh(kh) \sim kh, \lambda = \sqrt{T(gh)}`

        Parameters
        ----------
        SteepMax: float
            Wave breaking steepness estimate (e.g., 0.07).
        T_vals :np.array
            Array of T values [sec] at which to calculate the breaking height.
        depth: float
            Depth at site
            Note: if not inputted, the depth will tried to be grabbed from the respective
            buoy type's website.

        Returns
        -------
        SteepH: np.array
            H values [m] that correspond to the T_mesh values creating the
            steepness curve.
        T_steep: np.array
            T values [sec] over which the steepness curve is defined.

        Example
        -------

        To find limit the steepness of waves on a contour by breaking:
            
            import numpy as np
            import WDRT.ESSC as ESSC
            
            # Pull spectral data from NDBC website
            buoy46022 = ESSC.Buoy('46022','NDBC')
            buoy46022.fetchFromWeb()
            
            # Create PCA EA object for buoy
            pca46022 = ESSC.PCA(buoy46022)
            
            T_vals = np.arange(0.1, np.amax(buoy46022.T), 0.1)
            # Enter estimate of breaking steepness
            SteepMax = 0.07  # Reference DNV-RP-C205
            
            # Declare required parameters
            depth = 391.4  # Depth at measurement point (m)
            
            SteepH = pca46022.steepness(depth,SteepMax,T_vals)
        '''

        # Calculate the wavelength at a given depth at each value of T

        if depth == None:
            depth = self.__fetchDepth()
        lambdaT = []

        g = 9.81  # [m/s^2]
        omega = ((2 * np.pi) / T_vals)
        lambdaT = []

        for i in range(len(T_vals)):
            # Initialize kh using Eckart 1952 (mentioned in Holthuijsen pg. 124)
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
    
    def __fetchDepth(self):
        '''Obtains the depth from the website for a buoy (either NDBC or CDIP)'''
        if self.buoy.buoyType == "NDBC":
            url = "https://www.ndbc.noaa.gov/station_page.php?station=%s" % (46022)
            ndbcURL = requests.get(url)
            ndbcURL.raise_for_status()
            ndbcHTML = bs4.BeautifulSoup(ndbcURL.text, "lxml")
            header = ndbcHTML.find("b", text="Water depth:")
            return float(str(header.nextSibling).split()[0])
        elif self.buoy.buoyType == "CDIP":
            url = "http://cdip.ucsd.edu/cgi-bin/wnc_metadata?ARCHIVE/%sp1/%sp1_historic" % (self.buoy.buoyNum, self.buoy.buoyNum)
            cdipURL = requests.get(url)
            cdipURL.raise_for_status()
            cdipHTML = bs4.BeautifulSoup(cdipURL.text, "lxml")
            #Parse the table for the depth value
            depthString = str(cdipHTML.findChildren("td", {"class" : "plus"})[0])
            depthString = depthString.split("<br/>")[2]
            return float(re.findall(r"[-+]?\d*\.\d+|\d+", depthString)[0])



    def bootStrap(self, boot_size=1000, plotResults=True):
        '''Get 95% confidence bounds about a contour using the bootstrap
        method. Warning - this function is time consuming. Computation 
        time depends on selected boot_size.

        Parameters
        ----------
            boot_size: int (optional)
                Number of bootstrap samples that will be used to calculate 95%
                confidence interval. Should be large enough to calculate stable
                statistics. If left blank will be set to 1000.
            plotResults: boolean (optional)
                Option for showing plot of bootstrap confidence bounds. If left
                blank will be set to True and plot will be shown.

        Returns
        -------
            contourmean_Hs : nparray
                Hs values for mean contour calculated as the average over all
                bootstrap contours.
            contourmean_T : nparray
                T values for mean contour calculated as the average over all
                bootstrap contours.
        
        Example
        -------
        To generate 95% boostrap contours for a given contour method:
            
            import WDRT.ESSC as ESSC
            
            # Pull spectral data from NDBC website
            buoy46022 = ESSC.Buoy('46022','NDBC')
            buoy46022.fetchFromWeb()
            
            # Create PCA EA object for buoy
            pca46022 = ESSC.PCA(buoy46022)
            
            # Declare required parameters
            Time_SS = 1.  # Sea state duration (hrs)
            Time_r = 100  # Return periods (yrs) of interest
            nb_steps = 1000  # Enter discretization of the circle in the normal space (optional)
            
            # Contour generation
            Hs_Return, T_Return = pca46022.getContours(Time_SS, Time_r,nb_steps)
            
            # Calculate boostrap confidence interval
            contourmean_Hs, contourmean_T = pca46022.bootStrap(boot_size=10)
        '''
        #preallocates arrays
        n = len(self.buoy.Hs)
        Hs_Return_Boot = np.zeros([self.nb_steps,boot_size])
        T_Return_Boot = np.zeros([self.nb_steps,boot_size])
        buoycopy = copy.deepcopy(self.buoy);
        #creates copies of the data based on how it was modeled.
        for i in range(boot_size):
            boot_inds = np.random.randint(0, high=n, size=n)
            buoycopy.Hs = copy.deepcopy(self.buoy.Hs[boot_inds])
            buoycopy.T = copy.deepcopy(self.buoy.T[boot_inds])
            essccopy=None            
            if self.method == "Principle component analysis":
                essccopy = PCA(buoycopy, self.size_bin)
            elif self.method == "Gaussian Copula":
                essccopy = GaussianCopula(buoycopy, self.n_size, self.bin_1_limit, self.bin_step)
            elif self.method == "Rosenblatt":
                essccopy = Rosenblatt(buoycopy, self.n_size, self.bin_1_limit, self.bin_step)
            elif self.method == "Clayton Copula":
                essccopy = ClaytonCopula(buoycopy, self.n_size, self.bin_1_limit, self.bin_step)
            elif self.method == "Gumbel Copula":
                essccopy = GumbelCopula(buoycopy, self.n_size, self.bin_1_limit, self.bin_step, self.Ndata)
            elif self.method == "Non-parametric Gaussian Copula":
                essccopy = NonParaGaussianCopula(buoycopy, self.Ndata, self.max_T, self.max_Hs)
            elif self.method == "Non-parametric Clayton Copula":
                essccopy = NonParaClaytonCopula(buoycopy, self.Ndata, self.max_T, self.max_Hs)
            elif self.method == "Non-parametric Gumbel Copula":
                essccopy = NonParaGumbelCopula(buoycopy, self.Ndata, self.max_T, self.max_Hs)
            Hs_Return_Boot[:,i],T_Return_Boot[:,i] = essccopy.getContours(self.time_ss, self.time_r, self.nb_steps)

        #finds 95% CI values for wave height and energy
        contour97_5_Hs = np.percentile(Hs_Return_Boot,97.5,axis=1)
        contour2_5_Hs = np.percentile(Hs_Return_Boot,2.5,axis=1)
        contourmean_Hs = np.mean(Hs_Return_Boot, axis=1)

        contour97_5_T = np.percentile(T_Return_Boot,97.5,axis=1)
        contour2_5_T = np.percentile(T_Return_Boot,2.5,axis=1)
        contourmean_T = np.mean(T_Return_Boot, axis=1)

        self.contourMean_Hs = contourmean_Hs
        self.contourMean_T = contourmean_T
        #plotting function
        def plotResults():
            plt.figure()
            plt.plot(self.buoy.T, self.buoy.Hs, 'bo', alpha=0.1, label='NDBC data')
            plt.plot(self.T_ReturnContours, self.Hs_ReturnContours, 'k-', label='100 year contour')
            plt.plot(contour97_5_T, contour97_5_Hs, 'r--', label='95% bootstrap confidence interval')
            plt.plot(contour2_5_T, contour2_5_Hs, 'r--')
            plt.plot(contourmean_T, contourmean_Hs, 'r-', label='Mean bootstrap contour')
            plt.legend(loc='lower right', fontsize='small')
            plt.grid(True)
            plt.xlabel('Energy period, $T_e$ [s]')
            plt.ylabel('Sig. wave height, $H_s$ [m]')
            plt.show()
        if plotResults:
            plotResults()

        return contourmean_Hs, contourmean_T

    def outsidePoints(self):
        
        '''Determines which buoy observations are outside of a given contour.
        
        Parameters
        ----------
            None
            
        Returns
        -------
            outsideHs : nparray
                The Hs values of the observations that are outside of the contour
        
            outsideT : nparray
                The T values of the observations that are outside of the contour
        
        Example
        -------
        
            To get correseponding T and Hs arrays of observations that are outside
            of a given contour:
                
                import WDRT.ESSC as ESSC
                import numpy as np
                
                # Pull spectral data from NDBC website
                buoy46022 = ESSC.Buoy('46022','NDBC')
                buoy46022.fetchFromWeb()
                
                # Create PCA EA object for buoy
                rosen46022 = ESSC.Rosenblatt(buoy46022)
                
                # Declare required parameters
                Time_SS = 1.  # Sea state duration (hrs)
                Time_r = 100  # Return periods (yrs) of interest
                
                # Generate contour
                Hs_Return, T_Return = rosen46022.getContours(Time_SS, Time_r)
                
                # Return the outside point Hs/T combinations
                outsideT, outsideHs = rosen46022.outsidePoints()
        
        
        '''
        #checks if the contour type is a KDE contour - if so, finds the outside points for the KDE contour.
        if isinstance(self.T_ReturnContours,list): 
            contains_test = np.zeros(len(self.buoy.T),dtype=bool)
            for t,hs in zip(self.T_ReturnContours,self.Hs_ReturnContours):        
                path_contour = []        
                path_contour = matplotlib.path.Path(np.column_stack((t,hs)))
                contains_test = contains_test+path_contour.contains_points(np.column_stack((self.buoy.T,self.buoy.Hs)))
            out_inds = np.where(~contains_test)
        else: # For non-KDE methods (copulas, etc.)
            path_contour = matplotlib.path.Path(np.column_stack((self.T_ReturnContours,self.Hs_ReturnContours)))
            contains_test = path_contour.contains_points(np.column_stack((self.buoy.T,self.buoy.Hs)))
            out_inds = np.where(~contains_test)
        outsideHs =self.buoy.Hs[out_inds]
        outsideT = self.buoy.T[out_inds]
        
        return(outsideT, outsideHs)
    
    def contourIntegrator(self):    
             
        '''Calculates the area of the contour over the two-dimensional input
        space of interest. 
        
        Parameters
        ----------
            None
            
        Returns
        -------
            area : float
                The area of the contour in TxHs units. 
        
        Example
        -------
        
        To obtain the area of the contour:
            
            import WDRT.ESSC as ESSC
            
            # Pull spectral data from NDBC website
            buoy46022 = ESSC.Buoy('46022','NDBC')
            buoy46022.fetchFromWeb()
            
            # Create PCA EA object for buoy
            rosen46022 = ESSC.Rosenblatt(buoy46022)
            
            # Declare required parameters
            Time_SS = 1.  # Sea state duration (hrs)
            Time_r = 100  # Return periods (yrs) of interest
            
            # Generate contour
            Hs_Return, T_Return = rosen46022.getContours(Time_SS, Time_r)
            
            # Return the area of the contour
            rosenArea = rosen46022.contourIntegrator()
        
        
        '''        
        
        contourTs = self.T_ReturnContours
        contourHs = self.Hs_ReturnContours
    
        area = 0.5*np.abs(np.dot(contourTs,np.roll(contourHs,1))-np.dot(contourHs,np.roll(contourTs,1))) 
    
        return area
    
    def dataContour(self, tStepSize = 1, hsStepSize = .5):
    
        '''Creates a contour around the ordered pairs of buoy observations. How tightly
           the contour fits around the data will be determined by step size parameters.
           Please note that this function currently is in beta; it needs further work to be
           optimized for use.
        
        Parameters
        ----------
            tStepSize : float
                Determines how far to search for the next point in the T direction. 
                Smaller values will produce contours that follow the data more closely.
    
            hsStepSize : float
                Determines how far to search for the next point in the Hs direction. 
                Smaller values will produce contours that follow the data more closely.
                
        Returns
        -------
            dataBoundryHs : nparray
                The Hs values of the boundry observations
        
            dataBoundryT : nparray
                The Hs values of the boundry observations
        
        Example
        -------
        
        To get the corresponding data contour:
            
            import WDRT.ESSC as ESSC
            
            # Pull spectral data from NDBC website
            buoy46022 = ESSC.Buoy('46022','NDBC')
            buoy46022.fetchFromWeb()
            
            # Create PCA EA object for buoy
            rosen46022 = ESSC.Rosenblatt(buoy46022)
            
            # Calculate the data contour
            dataHs, dataT = rosen46022.dataContour(tStepSize = 1, hsStepSize = .5)
        
        
        '''    
    
        
        
        maxHs = max(self.buoy.Hs)
        minHs = min(self.buoy.Hs)
        
        sortedHsBuoy = copy.deepcopy(self.buoy)
        sortedTBuoy = copy.deepcopy(self.buoy)
        
        sortedTIndex = sorted(range(len(self.buoy.T)),key=lambda x:self.buoy.T[x])
        sortedHsIndex = sorted(range(len(self.buoy.Hs)),key=lambda x:self.buoy.Hs[x])
        
        sortedHsBuoy.Hs = self.buoy.Hs[sortedHsIndex]
        sortedHsBuoy.T = self.buoy.T[sortedHsIndex]
        sortedTBuoy.Hs = self.buoy.Hs[sortedTIndex]
        sortedTBuoy.T = self.buoy.T[sortedTIndex]
        
        hsBin1 = []
        hsBin2 = []
        hsBin3 = []
        hsBin4 = []
        tBin1 = []
        tBin2 = []
        tBin3 = []
        tBin4 = []
    
        
        startingPoint = sortedTBuoy.T[0]
        hsBin4.append(sortedTBuoy.Hs[0])
        tBin4.append(sortedTBuoy.T[0])
        while True:
            tempNextBinTs = sortedTBuoy.T[sortedTBuoy.T < startingPoint + tStepSize]    
            tempNextBinHs = sortedTBuoy.Hs[sortedTBuoy.T < startingPoint + tStepSize]
            nextBinTs = tempNextBinTs[tempNextBinTs > startingPoint]
            nextBinHs = tempNextBinHs[tempNextBinTs > startingPoint]
            
            try:
                nextHs = max(nextBinHs)
                nextT = nextBinTs[nextBinHs.argmax(axis=0)]
                
                hsBin4.append(nextHs)
                tBin4.append(nextT)
            
                startingPoint = nextT
            except ValueError:
                startingPoint += tStepSize
                break
     
            if nextHs == maxHs:
                break
        
        startingPoint = sortedTBuoy.T[0]
        hsBin1.append(sortedTBuoy.Hs[0])
        tBin1.append(sortedTBuoy.T[0])
        while True:
            tempNextBinTs = sortedTBuoy.T[sortedTBuoy.T < startingPoint + tStepSize]    
            tempNextBinHs = sortedTBuoy.Hs[sortedTBuoy.T < startingPoint + tStepSize]
            nextBinTs = tempNextBinTs[tempNextBinTs > startingPoint]
            nextBinHs = tempNextBinHs[tempNextBinTs > startingPoint]
            
            try:
                nextHs = min(nextBinHs)
                nextT = nextBinTs[nextBinHs.argmin(axis=0)]
                
                hsBin1.append(nextHs)
                tBin1.append(nextT)
            
                startingPoint = nextT
            except ValueError:
                startingPoint += tStepSize
                break
            
            
            if nextHs == minHs:
                break
        
        
        startingPoint = sortedHsBuoy.Hs[sortedHsBuoy.T.argmax(axis=0)]
        hsBin3.append(sortedHsBuoy.Hs[sortedHsBuoy.T.argmax(axis=0)])
        tBin3.append(sortedHsBuoy.T[sortedHsBuoy.T.argmax(axis=0)])
        while True:
            tempNextBinTs = sortedHsBuoy.T[sortedHsBuoy.Hs < startingPoint + hsStepSize]    
            tempNextBinHs = sortedHsBuoy.Hs[sortedHsBuoy.Hs < startingPoint + hsStepSize]
            nextBinTs = tempNextBinTs[tempNextBinHs > startingPoint]
            nextBinHs = tempNextBinHs[tempNextBinHs > startingPoint]
        
            try:
                nextT = max(nextBinTs)
                nextHs = nextBinHs[nextBinTs.argmax(axis=0)]
                
                if nextHs not in hsBin4 and nextHs not in hsBin1:
                
                    hsBin3.append(nextHs)
                    tBin3.append(nextT)
        
                startingPoint = nextHs
            
            except ValueError:
                startingPoint += hsStepSize
                break
    
            if nextHs == maxHs:
                break
     
        
        startingPoint = sortedHsBuoy.Hs[sortedHsBuoy.T.argmax(axis=0)]
        while True:
            tempNextBinTs = sortedHsBuoy.T[sortedHsBuoy.Hs > startingPoint - hsStepSize]    
            tempNextBinHs = sortedHsBuoy.Hs[sortedHsBuoy.Hs > startingPoint - hsStepSize]
            nextBinTs = tempNextBinTs[tempNextBinHs < startingPoint]
            nextBinHs = tempNextBinHs[tempNextBinHs < startingPoint]
        
            try:    
                nextT = max(nextBinTs)
                nextHs = nextBinHs[nextBinTs.argmax(axis=0)]
            
                if nextHs not in hsBin1 and nextHs not in hsBin4:
            
                    hsBin2.append(nextHs)
                    tBin2.append(nextT)
        
                startingPoint = nextHs
            except ValueError:
                startingPoint = startingPoint - hsStepSize
                break
    
            
            if nextHs == minHs:
                break
        
        hsBin2 = hsBin2[::-1] # Reverses the order of the array
        tBin2 = tBin2[::-1]
        hsBin4 = hsBin4[::-1] # Reverses the order of the array
        tBin4 = tBin4[::-1]
        
        dataBoundryHs = np.concatenate((hsBin1,hsBin2,hsBin3,hsBin4),axis = 0)
        dataBoundryT = np.concatenate((tBin1,tBin2,tBin3,tBin4),axis = 0)
        dataBoundryHs = dataBoundryHs[::-1]
        dataBoundryT = dataBoundryT[::-1]
        
        return(dataBoundryHs, dataBoundryT)
    


    def __getCopulaParams(self,n_size,bin_1_limit,bin_step):
        sorted_idx = sorted(range(len(self.buoy.Hs)),key=lambda x:self.buoy.Hs[x])
        Hs = self.buoy.Hs[sorted_idx]
        T = self.buoy.T[sorted_idx]

        # Estimate parameters for Weibull distribution for component 1 (Hs) using MLE
        # Estimate parameters for Lognormal distribution for component 2 (T) using MLE
        para_dist_1=stats.exponweib.fit(Hs,floc=0,fa=1)
        para_dist_2=stats.norm.fit(np.log(T))

        # Binning
        ind = np.array([])
        ind = np.append(ind,sum(Hs_val <= bin_1_limit for Hs_val in Hs))
        # Make sure first bin isn't empty or too small to avoid errors        
        while ind == 0 or ind < n_size:         
            ind = np.array([])    
            bin_1_limit = bin_1_limit + bin_step
            ind = np.append(ind,sum(Hs_val <= bin_1_limit for Hs_val in Hs))        
        for i in range(1,200):
            bin_i_limit = bin_1_limit+bin_step*(i)
            ind = np.append(ind,sum(Hs_val <= bin_i_limit for Hs_val in Hs))
            if (ind[i-0]-ind[i-1]) < n_size:
                break
    
        # Parameters for conditional distribution of T|Hs for each bin
        num=len(ind) # num+1: number of bins
        para_dist_cond = []
        hss = []

        para_dist_cond.append(stats.norm.fit(np.log(T[range(0,int(ind[0]))])))  # parameters for first bin
        hss.append(np.mean(Hs[range(0,int(ind[0])-1)])) # mean of Hs (component 1 for first bin)
        para_dist_cond.append(stats.norm.fit(np.log(T[range(0,int(ind[1]))]))) # parameters for second bin
        hss.append(np.mean(Hs[range(0,int(ind[1])-1)])) # mean of Hs (component 1 for second bin)

        for i in range(2,num):
            para_dist_cond.append(stats.norm.fit(np.log(T[range(int(ind[i-2]),int(ind[i]))])));
            hss.append(np.mean(Hs[range(int(ind[i-2]),int(ind[i]))]))

        # Estimate coefficient using least square solution (mean: third order, sigma: 2nd order)
        para_dist_cond.append(stats.norm.fit(np.log(T[range(int(ind[num-2]),int(len(Hs)))])));  # parameters for last bin
        hss.append(np.mean(Hs[range(int(ind[num-2]),int(len(Hs)))])) # mean of Hs (component 1 for last bin)

        para_dist_cond = np.array(para_dist_cond)
        hss = np.array(hss)

        phi_mean = np.column_stack((np.ones(num+1),hss[:],hss[:]**2,hss[:]**3))
        phi_std = np.column_stack((np.ones(num+1),hss[:],hss[:]**2))

        # Estimate coefficients of mean of Ln(T|Hs)(vector 4x1) (cubic in Hs)
        mean_cond = np.linalg.lstsq(phi_mean,para_dist_cond[:,0])[0]
        # Estimate coefficients of standard deviation of Ln(T|Hs) (vector 3x1) (quadratic in Hs)
        std_cond = np.linalg.lstsq(phi_std,para_dist_cond[:,1])[0]

        return para_dist_1, para_dist_2, mean_cond, std_cond
        
    def __getNonParaCopulaParams(self,Ndata, max_T, max_Hs):
        sorted_idx = sorted(range(len(self.buoy.Hs)),key=lambda x:self.buoy.Hs[x])
        Hs = self.buoy.Hs[sorted_idx]
        T = self.buoy.T[sorted_idx]
        
        # Calcualte KDE bounds (this may be added as an input later)
        min_limit_1 = 0
        max_limit_1 = max_Hs
        min_limit_2 = 0
        max_limit_2 = max_T   
        
        # Discretize for KDE
        pts_hs = np.linspace(min_limit_1, max_limit_1, self.Ndata) 
        pts_t = np.linspace(min_limit_2, max_limit_2, self.Ndata)
        
        # Calculate optimal bandwidth for T and Hs
        sig = robust.scale.mad(T)
        num = float(len(T))
        bwT = sig*(4.0/(3.0*num))**(1.0/5.0)
        
        sig = robust.scale.mad(Hs)
        num = float(len(Hs))
        bwHs = sig*(4.0/(3.0*num))**(1.0/5.0)
        
        # Nonparametric PDF for T
        temp = sm.nonparametric.KDEUnivariate(T)
        temp.fit(bw = bwT)
        f_t = temp.evaluate(pts_t)
        
        # Nonparametric CDF for Hs
        temp = sm.nonparametric.KDEUnivariate(Hs)
        temp.fit(bw = bwHs)
        tempPDF = temp.evaluate(pts_hs)
        F_hs = tempPDF/sum(tempPDF)
        F_hs = np.cumsum(F_hs)
        
        # Nonparametric CDF for T
        F_t = f_t/sum(f_t)
        F_t = np.cumsum(F_t)
        
        nonpara_dist_1 = np.transpose(np.array([pts_hs, F_hs]))
        nonpara_dist_2 = np.transpose(np.array([pts_t, F_t]))
        nonpara_pdf_2 = np.transpose(np.array([pts_t, f_t]))
        
        return nonpara_dist_1, nonpara_dist_2, nonpara_pdf_2
        
    def __gumbelCopula(self, u, alpha):
        ''' Calculates the Gumbel copula density
        Parameters
        ----------
        u: np.array
                    Vector of equally spaced points between 0 and twice the
                    maximum value of T.
        alpha: float
                    Copula parameter. Must be greater than or equal to 1.
        Returns
        -------
        y: np.array
                   Copula density function.
        '''
        #Ignore divide by 0 warnings and resulting NaN warnings
        np.seterr(all='ignore')        
        v = -np.log(u)
        v = np.sort(v, axis=0)
        vmin = v[0, :]
        vmax = v[1, :]
        nlogC = vmax * (1 + (vmin / vmax) ** alpha) ** (1 / alpha)
        y = (alpha - 1 +nlogC)*np.exp(-nlogC+np.sum((alpha-1)*np.log(v)+v, axis =0) +(1-2*alpha)*np.log(nlogC))
        np.seterr(all='warn')

        return(y)        

class PCA(EA):
    def __init__(self, buoy, size_bin=250.):
        '''
        Create a PCA EA class for a buoy object. Contours generated under this
        class will use principal component analysis (PCA) with improved 
        distribution fitting (Eckert et. al 2015) and the I-FORM.
        
        Parameters
        ___________
            size_bin : float
                chosen bin size
            buoy : NDBCData
                ESSC.Buoy Object
        '''
        self.method = "Principle component analysis"
        self.buoy = buoy
        if size_bin > len(buoy.Hs)*0.25:
            self.size_bin = len(buoy.Hs)*0.25
            print(round(len(buoy.Hs)*0.25,2),'is the max bin size for this buoy. The bin size has been set to this amount.')
        else:
            self.size_bin = size_bin
            
        self.Hs_ReturnContours = None
        self.Hs_SampleCA = None
        self.Hs_SampleFSS = None

        self.T_ReturnContours = None
        self.T_SampleCA = None
        self.T_SampleFSS = None

        self.Weight_points = None

        self.coeff, self.shift, self.comp1_params, self.sigma_param, self.mu_param = self.__generateParams(self.size_bin)

    def __generateParams(self, size_bin=250.0):
        pca = skPCA(n_components=2)
        pca.fit(np.array((self.buoy.Hs - self.buoy.Hs.mean(axis=0), self.buoy.T - self.buoy.T.mean(axis=0))).T)
        coeff = abs(pca.components_)  # Apply correct/expected sign convention
        coeff[1, 1] = -1.0 * coeff[1, 1]  # Apply correct/expected sign convention

        Comp1_Comp2 = np.dot (np.array((self.buoy.Hs, self.buoy.T)).T, coeff)

        shift = abs(min(Comp1_Comp2[:, 1])) + 0.1  # Calculate shift


        shift = abs(min(Comp1_Comp2[:, 1])) + 0.1  # Calculate shift
        # Apply shift to Component 2 to make all values positive
        Comp1_Comp2[:, 1] = Comp1_Comp2[:, 1] + shift

        Comp1_Comp2_sort = Comp1_Comp2[Comp1_Comp2[:, 0].argsort(), :]

        # Fitting distribution of component 1
        comp1_params = stats.invgauss.fit(Comp1_Comp2_sort[:, 0], floc=0)

        n_data = len(self.buoy.Hs)  # Number of observations

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

        mu_param, pcov = optim.curve_fit(self.__mu_fcn,
                                                 Comp1_mean.T, Comp2_bins_params[0, :])

        sigma_param = self.__sigma_fits(Comp1_mean, Comp2_bins_params[1, :])

        return coeff, shift, comp1_params, sigma_param, mu_param

    def _saveParams(self, groupObj):
        if('nb_steps' in groupObj):
            groupObj['nb_steps'][...] = self.nb_steps
        else:
            groupObj.create_dataset('nb_steps', data=self.nb_steps)
        
        if('time_r' in groupObj):
            groupObj['time_r'][...] = self.time_r
        else:
            groupObj.create_dataset('time_r', data=self.time_r)
        
        if('time_ss' in groupObj):
            groupObj['time_ss'][...] = self.time_ss
        else:       
            groupObj.create_dataset('time_ss', data=self.time_ss)

        if('coeff' in groupObj):
            groupObj['coeff'][...] = self.coeff
        else:
            groupObj.create_dataset('coeff', data=self.coeff)
        
        if('shift' in groupObj):
            groupObj['shift'][...] = self.shift
        else:
            groupObj.create_dataset('shift', data=self.shift)
        
        if('comp1_params' in groupObj):
            groupObj['comp1_params'][...] = self.comp1_params
        else:
            groupObj.create_dataset('comp1_params', data=self.comp1_params)
        
        if('sigma_param' in groupObj):
            groupObj['sigma_param'][...] = self.sigma_param
        else:
            groupObj.create_dataset('sigma_param', data=self.sigma_param)
        
        if('mu_param' in groupObj):
            groupObj['mu_param'][...] = self.mu_param
        else:
            groupObj.create_dataset('mu_param', data=self.mu_param)

    def getContours(self, time_ss, time_r, nb_steps=1000):
        '''WDRT Extreme Sea State PCA Contour function
        This function calculates environmental contours of extreme sea states using
        principal component analysis and the inverse first-order reliability
        method.

        Parameters
        ___________
        time_ss : float
            Sea state duration (hours) of measurements in input.
        time_r : np.array
            Desired return period (years) for calculation of environmental
            contour, can be a scalar or a vector.
        nb_steps : int
            Discretization of the circle in the normal space used for
            inverse FORM calculation.

        Returns
        -------
        Hs_Return : np.array
            Calculated Hs values along the contour boundary following
            return to original input orientation.
        T_Return : np.array
           Calculated T values along the contour boundary following
           return to original input orientation.
        nb_steps : float
            Discretization of the circle in the normal space

        Example
        -------
        To obtain the contours for a NDBC buoy::
            
            import WDRT.ESSC as ESSC
            
            # Pull spectral data from NDBC website
            buoy46022 = ESSC.Buoy('46022','NDBC')
            buoy46022.fetchFromWeb()
            
            # Create PCA EA object for buoy
            pca46022 = ESSC.PCA(buoy46022)
            
            # Declare required parameters
            Time_SS = 1.  # Sea state duration (hrs)
            Time_r = 100  # Return periods (yrs) of interest
            nb_steps = 1000  # Enter discretization of the circle in the normal space (optional)
            
            # Contour generation example
            Hs_Return, T_Return = pca46022.getContours(Time_SS, Time_r,nb_steps)
        '''

        self.time_ss = time_ss
        self.time_r = time_r
        self.nb_steps = nb_steps

        # IFORM
        # Failure probability for the desired return period (time_R) given the
        # duration of the measurements (time_ss)
        p_f = 1 / (365 * (24 / time_ss) * time_r)
        beta = stats.norm.ppf((1 - p_f), loc=0, scale=1)  # Reliability
        theta = np.linspace(0, 2 * np.pi, num = nb_steps)
        # Vary U1, U2 along circle sqrt(U1^2+U2^2)=beta
        U1 = beta * np.cos(theta)
        U2 = beta * np.sin(theta)
        # Calculate C1 values along the contour
        Comp1_R = stats.invgauss.ppf(stats.norm.cdf(U1, loc=0, scale=1),
                                     mu= self.comp1_params[0], loc=0,
                                     scale= self.comp1_params[2])
        # Calculate mu values at each point on the circle
        mu_R = self.__mu_fcn(Comp1_R, self.mu_param[0], self.mu_param[1])
        # Calculate sigma values at each point on the circle
        sigma_R = self.__sigma_fcn(self.sigma_param, Comp1_R)
        # Use calculated mu and sigma values to calculate C2 along the contour
        Comp2_R = stats.norm.ppf(stats.norm.cdf(U2, loc=0, scale=1),
                                 loc=mu_R, scale=sigma_R)

        # Calculate Hs and T along the contour
        Hs_Return, T_Return = self.__princomp_inv(Comp1_R, Comp2_R, self.coeff, self.shift)
        Hs_Return = np.maximum(0, Hs_Return)  # Remove negative values
        self.Hs_ReturnContours = Hs_Return
        self.T_ReturnContours = T_Return
        return Hs_Return, T_Return

    def getSamples(self, num_contour_points, contour_returns, random_seed=None):
        '''WDRT Extreme Sea State Contour Sampling function.
        This function calculates samples of Hs and T using the EA function to
        sample between contours of user-defined return periods.

        Parameters
        ----------
        num_contour_points : int
            Number of sample points to be calculated per contour interval.
        contour_returns: np.array
            Vector of return periods that define the contour intervals in
            which samples will be taken. Values must be greater than zero and
            must be in increasing order.
        random_seed: int (optional)
            Random seed for sample generation, required for sample
            repeatability. If left blank, a seed will automatically be
            generated.

        Returns
        -------
        Hs_Samples: np.array
            Vector of Hs values for each sample point.
        Te_Samples: np.array
            Vector of Te values for each sample point.
        Weight_points: np.array
            Vector of probabilistic weights for each sampling point
            to be used in risk calculations.

        Example
        -------
        To get weighted samples from a set of contours::

                import numpy as np
                import WDRT.ESSC as ESSC
                
                # Pull spectral data from NDBC website
                buoy46022 = ESSC.Buoy('46022','NDBC')
                buoy46022.fetchFromWeb()
                
                # Create PCA EA object for buoy
                pca46022 = ESSC.PCA(buoy46022)
                
                # Declare required parameters
                Time_SS = 1.  # Sea state duration (hrs)
                Time_r = 100  # Return periods (yrs) of interest
                num_contour_points = 10 # Number of points to be sampled for each contour interval
                contour_returns = np.array([0.001, 0.01, 0.05, 0.1, 0.5, 1, 5, 10, 50, 100])
                
                # Calculate contour to save required variables to PCA EA object
                pca46022.getContours(Time_SS, Time_r,nb_steps)
                
                # Probabilities defining sampling contour bounds.
                random_seed = 2  # Random seed for sample generation
                
                # Get samples for a full sea state long term analysis
                Hs_sampleFSS, T_sampleFSS, Weight_sampleFSS = pca46022.getSamples(num_contour_points,
                                                             contour_returns, random_seed)
        '''

        # Calculate line where Hs = 0 to avoid sampling Hs in negative space
        Te_zeroline = np.linspace(2.5, 30, 1000)
        Te_zeroline = np.transpose(Te_zeroline)
        Hs_zeroline = np.zeros(len(Te_zeroline))

        # Transform zero line into principal component space
        Comp_zeroline = np.dot(np.transpose(np.vstack([Hs_zeroline, Te_zeroline])),
                               self.coeff)
        Comp_zeroline[:, 1] = Comp_zeroline[:, 1] + self.shift

        # Find quantiles along zero line
        C1_zeroline_prob = stats.invgauss.cdf(Comp_zeroline[:, 0],
                                              mu = self.comp1_params[0], loc=0,
                                              scale = self.comp1_params[2])
        mu_zeroline = self.__mu_fcn(Comp_zeroline[:, 0], self.mu_param[0], self.mu_param[1])
        sigma_zeroline = self.__sigma_fcn(self.sigma_param, Comp_zeroline[:, 0])
        C2_zeroline_prob = stats.norm.cdf(Comp_zeroline[:, 1],
                                          loc=mu_zeroline, scale=sigma_zeroline)
        C1_normzeroline = stats.norm.ppf(C1_zeroline_prob, 0, 1)
        C2_normzeroline = stats.norm.ppf(C2_zeroline_prob, 0, 1)

        contour_probs = 1 / (365 * (24 / self.time_ss) * contour_returns)
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


        Sample_alpha, Sample_beta, Weight_points = self.__generateData(beta_lines,
            Rho_zeroline, Theta_zeroline, num_contour_points,contour_probs,random_seed)

        Hs_Sample, T_Sample = self.__transformSamples(Sample_alpha, Sample_beta)

        self.Hs_SampleFSS = Hs_Sample
        self.T_SampleFSS = T_Sample
        self.Weight_SampleFSS = Weight_points

        return Hs_Sample, T_Sample, Weight_points
    
    def plotSampleData(self):
        """
        Display a plot of the 100-year return contour, full sea state samples
        and contour samples
        """
        plt.figure()
        plt.plot(self.buoy.T, self.buoy.Hs, 'bo', alpha=0.1, label='NDBC data')
        plt.plot(self.T_ReturnContours, self.Hs_ReturnContours, 'k-', label='100 year contour')
        plt.plot(self.T_SampleFSS, self.Hs_SampleFSS, 'ro', label='full sea state samples')
        plt.plot(self.T_SampleCA, self.Hs_SampleCA, 'y^', label='contour approach samples')
        plt.legend(loc='lower right', fontsize='small')
        plt.grid(True)
        plt.xlabel('Energy period, $T_e$ [s]')
        plt.ylabel('Sig. wave height, $H_s$ [m]')
        plt.show() 

    def __generateData(self, beta_lines, Rho_zeroline, Theta_zeroline, num_contour_points, contour_probs, random_seed):
        """
        Calculates radius, angle, and weight for each sample point
        """
        
        ''' Data generating function that calculates the radius, angle, and
        weight for each sample point.
        Parameters
        ----------
        beta_lines: np.array
               Array of mu fitting function parameters.
        Rho_zeroline: np.array
               array of radii
        Theta_zeroline: np.array
        num_contour_points: np.array
        contour_probs: np.array
        random_seed: int
            seed for generating random data.
        Returns
        -------
        Sample_alpha: np.array
                Array of fitted sample angle values.
        Sample_beta: np.array
                Array of fitted sample radius values.
        Weight_points: np.array
                Array of weights for each point.
        '''
        np.random.seed(random_seed)

        num_samples = (len(beta_lines) - 1) * num_contour_points
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

        return Sample_alpha, Sample_beta, Weight_points

    def __transformSamples(self, Sample_alpha, Sample_beta):
        Sample_U1 = Sample_beta * np.cos(Sample_alpha)
        Sample_U2 = Sample_beta * np.sin(Sample_alpha)

        # Sample transformation to principal component space
        Comp1_sample = stats.invgauss.ppf(stats.norm.cdf(Sample_U1, loc=0, scale=1),
                                          mu=self.comp1_params[0], loc=0,
                                          scale=self.comp1_params[2])
        mu_sample = self.__mu_fcn(Comp1_sample, self.mu_param[0], self.mu_param[1])
        # Calculate sigma values at each point on the circle
        sigma_sample = self.__sigma_fcn(self.sigma_param, Comp1_sample)
        # Use calculated mu and sigma values to calculate C2 along the contour
        Comp2_sample = stats.norm.ppf(stats.norm.cdf(Sample_U2, loc=0, scale=1),
                                      loc=mu_sample, scale=sigma_sample)
        # Sample transformation into Hs-T space
        Hs_Sample, T_Sample = self.__princomp_inv(
            Comp1_sample, Comp2_sample, self.coeff, self.shift)

        return Hs_Sample, T_Sample


    def __mu_fcn(self, x, mu_p_1, mu_p_2):
        ''' Linear fitting function for the mean(mu) of Component 2 normal
        distribution as a function of the Component 1 mean for each bin.
        Used in the EA and getSamples functions.
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


    def __sigma_fcn(self,sig_p, x):
        '''Quadratic fitting formula for the standard deviation(sigma) of Component
        2 normal distribution as a function of the Component 1 mean for each bin.
        Used in the EA and getSamples functions.
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


    def __princomp_inv(self, princip_data1, princip_data2, coeff, shift):
        '''Takes the inverse of the principal component rotation given data,
        coefficients, and shift. Used in the EA and getSamples functions.
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

    def __betafcn(self, sig_p, rho):
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

    # Sigma function sigma_fcn defined outside of EA function

    def __objfun(self, sig_p, x, y_actual):
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
        obj_fun_result = np.sum((self.__sigma_fcn(sig_p, x) - y_actual)**2)
        return obj_fun_result  # Sum of least square error

    def __objfun_penalty(self, sig_p, x, y_actual, Beta1, Beta2):
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
        penalty_fcn = (self.__objfun(sig_p, x, y_actual) + Beta1 * (-sig_p[2])**2 +
                       Beta2 * (-sig_p[2] + (sig_p[1]**2) / (4 * sig_p[0]))**2)
        return penalty_fcn

    def __sigma_fits(self, Comp1_mean, sigma_vals):
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
        Beta1, Beta2 = self.__betafcn(sig_0, rho)
        # Initial search for minimum value using initial guess
        sig_1 = optim.fmin(func=self.__objfun_penalty, x0=sig_0,
                           args=(Comp1_mean, sigma_vals, Beta1, Beta2), disp=False)
        # While either the difference between iterations or the difference in
        # objective function evaluation is greater than the tolerance, continue
        # iterating
        while (np.amin(abs(sig_1 - sig_0)) > epsilon and
               abs(self.__objfun(sig_1, Comp1_mean, sigma_vals) -
                   self.__objfun(sig_0, Comp1_mean, sigma_vals)) > epsilon):
            sig_0 = sig_1
            # Calculate penalties for this iteration
            Beta1, Beta2 = self.__betafcn(sig_0, rho)
            # Find a new minimum
            sig_1 = optim.fmin(func=self.__objfun_penalty, x0=sig_0,
                               args=(Comp1_mean, sigma_vals, Beta1, Beta2), disp=False)
            rho = 10 * rho  # Increase penalization
        sig_final = sig_1
        return sig_final


class GaussianCopula(EA):
    '''Create a GaussianCopula EA class for a buoy object. Contours generated 
    under this class will use a Gaussian copula.'''
    def __init__(self, buoy, n_size=40., bin_1_limit=1., bin_step=0.25):
        '''
        Parameters
        ----------
            buoy : NDBCData
                ESSC.Buoy Object
            n_size: float
                minimum bin size used for Copula contour methods
            bin_1_limit: float
                maximum value of Hs for the first bin
            bin_step: float
                overlap interval for each bin
        '''
        self.method = "Gaussian Copula"
        self.buoy = buoy
        self.n_size = n_size
        self.bin_1_limit = bin_1_limit
        self.bin_step = bin_step

        self.Hs_ReturnContours = None
#        self.Hs_SampleCA = None
#        self.Hs_SampleFSS = None

        self.T_ReturnContours = None
#        self.T_SampleCA = None
#        self.T_SampleFSS = None

#        self.Weight_points = None

#        self.coeff, self.shift, self.comp1_params, self.sigma_param, self.mu_param = self.__generateParams(size_bin)
        self.para_dist_1,self.para_dist_2,self.mean_cond,self.std_cond = self._EA__getCopulaParams(n_size,bin_1_limit,bin_step)

    def getContours(self, time_ss, time_r, nb_steps = 1000):
        '''WDRT Extreme Sea State Gaussian Copula Contour function.
        This function calculates environmental contours of extreme sea states using
        a Gaussian copula and the inverse first-order reliability
        method.

        Parameters
        ___________
        time_ss : float
            Sea state duration (hours) of measurements in input.
        time_r : np.array
            Desired return period (years) for calculation of environmental
            contour, can be a scalar or a vector.
        nb_steps : float
            Discretization of the circle in the normal space used for
            inverse FORM calculation.

        Returns
        -------
        Hs_Return : np.array
            Calculated Hs values along the contour boundary following
            return to original input orientation.
        T_Return : np.array
           Calculated T values along the contour boundary following
           return to original input orientation.
        nb_steps : float
            Discretization of the circle in the normal space

        Example
        -------
        To obtain the contours for a NDBC buoy::
            
            import WDRT.ESSC as ESSC
            
            # Pull spectral data from NDBC website
            buoy46022 = ESSC.Buoy('46022','NDBC')
            buoy46022.fetchFromWeb()
            
            # Create Environtmal Analysis object using above parameters
            Gauss46022 = ESSC.GaussianCopula(buoy46022)
            
            # Declare required parameters
            Time_SS = 1.  # Sea state duration (hrs)
            Time_r = 100  # Return periods (yrs) of interest
            nb_steps = 1000  # Enter discretization of the circle in the normal space (optional)
            
            # Gaussian copula contour generation example
            Hs_Return, T_Return = Gauss46022.getContours(Time_SS, Time_r,nb_steps)
        '''
        self.time_ss = time_ss
        self.time_r = time_r
        self.nb_steps = nb_steps

        p_f = 1 / (365 * (24 / time_ss) * time_r)
        beta = stats.norm.ppf((1 - p_f), loc=0, scale=1)  # Reliability
        theta = np.linspace(0, 2 * np.pi, num = nb_steps)
        # Vary U1, U2 along circle sqrt(U1^2+U2^2)=beta
        U1 = beta * np.cos(theta)
        U2 = beta * np.sin(theta)

        comp_1 = stats.exponweib.ppf(stats.norm.cdf(U1),a=self.para_dist_1[0],c=self.para_dist_1[1],loc=self.para_dist_1[2],scale=self.para_dist_1[3])

        tau = stats.kendalltau(self.buoy.T,self.buoy.Hs)[0] # Calculate Kendall's tau
        rho_gau=np.sin(tau*np.pi/2.)

        z2_Gau=stats.norm.cdf(U2*np.sqrt(1.-rho_gau**2.)+rho_gau*U1);
        comp_2_Gaussian = stats.lognorm.ppf(z2_Gau,s=self.para_dist_2[1],loc=0,scale=np.exp(self.para_dist_2[0])) #lognormalinverse

        Hs_Return = comp_1
        T_Return = comp_2_Gaussian

        self.Hs_ReturnContours = Hs_Return
        self.T_ReturnContours = T_Return
        return Hs_Return, T_Return

    def getSamples(self):
        '''Currently not implemented in this version.'''
        raise NotImplementedError

    def _saveParams(self, groupObj):
        groupObj.create_dataset('n_size', data=self.n_size)
        groupObj.create_dataset('bin_1_limit', data=self.bin_1_limit)
        groupObj.create_dataset('bin_step', data=self.bin_step)
        groupObj.create_dataset('para_dist_1', data=self.para_dist_1)
        groupObj.create_dataset('para_dist_2', data=self.para_dist_2)
        groupObj.create_dataset('mean_cond', data=self.mean_cond)
        groupObj.create_dataset('std_cond', data=self.std_cond)


class Rosenblatt(EA):
    '''Create a Rosenblatt EA class for a buoy object. Contours generated 
    under this class will use a Rosenblatt transformation and the I-FORM.'''    
    def __init__(self, buoy, n_size=50., bin_1_limit= .5, bin_step=0.25):
        '''
        Parameters
        ----------
            buoy : NDBCData
                ESSC.Buoy Object
            n_size: float
                minimum bin size used for Copula contour methods
            bin_1_limit: float
                maximum value of Hs for the first bin
            bin_step: float
                overlap interval for each bin
        '''
        self.method = "Rosenblatt"
        self.buoy = buoy
        
        if n_size > 100:
            self.n_size = 100
            print(100,'is the maximum "minimum bin size" for this buoy. The minimum bin size has been set to this amount.')
        else:
            self.n_size = n_size
        
        if bin_step > max(buoy.Hs)*.1:
            self.bin_step = max(buoy.Hs)*.1
            print(round(max(buoy.Hs)*.1,2),'is the maximum bin overlap for this buoy. The bin overlap has been set to this amount.')
        else:
            self.bin_step = bin_step

        if bin_1_limit  > max(buoy.Hs)*.25:
            self.bin_1_limit = max(buoy.Hs)*.25
            print(round(max(buoy.Hs)*.25,2),'is the maximum limit for the first for this buoy. The first bin limit has been set to this amount.')
        else:
            self.bin_1_limit = bin_1_limit           
        
        

        self.Hs_ReturnContours = None
#        self.Hs_SampleCA = None
#        self.Hs_SampleFSS = None

        self.T_ReturnContours = None
#        self.T_SampleCA = None
#        self.T_SampleFSS = None

#        self.Weight_points = None

#        self.coeff, self.shift, self.comp1_params, self.sigma_param, self.mu_param = self.__generateParams(size_bin)
        self.para_dist_1,self.para_dist_2,self.mean_cond,self.std_cond = self._EA__getCopulaParams(self.n_size,self.bin_1_limit,self.bin_step)

    def getContours(self, time_ss, time_r, nb_steps = 1000):
        '''WDRT Extreme Sea State Rosenblatt Copula Contour function.
        This function calculates environmental contours of extreme sea states using
        a Rosenblatt transformation and the inverse first-order reliability
        method.

        Parameters
        ___________
        time_ss : float
            Sea state duration (hours) of measurements in input.
        time_r : np.array
            Desired return period (years) for calculation of environmental
            contour, can be a scalar or a vector.
        nb_steps : float
            Discretization of the circle in the normal space used for
            inverse FORM calculation.

        Returns
        -------
        Hs_Return : np.array
            Calculated Hs values along the contour boundary following
            return to original input orientation.
        T_Return : np.array
           Calculated T values along the contour boundary following
           return to original input orientation.
        nb_steps : float
            Discretization of the circle in the normal space

        Example
        -------
        To obtain the contours for a NDBC buoy::
            
            import WDRT.ESSC as ESSC
            
            # Pull spectral data from NDBC website
            buoy46022 = ESSC.Buoy('46022','NDBC')
            buoy46022.fetchFromWeb()
            
            # Create Environtmal Analysis object using above parameters
            Rosen46022 = ESSC.Rosenblatt(buoy46022)
            
            # Declare required parameters
            Time_SS = 1.  # Sea state duration (hrs)
            Time_r = 100  # Return periods (yrs) of interest
            nb_steps = 1000  # Enter discretization of the circle in the normal space (optional)
            
            # Rosenblatt contour generation example
            Hs_Return, T_Return = Rosen46022.getContours(Time_SS, Time_r,nb_steps)
        '''
        self.time_ss = time_ss
        self.time_r = time_r
        self.nb_steps = nb_steps

        p_f = 1 / (365 * (24 / time_ss) * time_r)
        beta = stats.norm.ppf((1 - p_f), loc=0, scale=1)  # Reliability
        theta = np.linspace(0, 2 * np.pi, num = nb_steps)
        # Vary U1, U2 along circle sqrt(U1^2+U2^2)=beta
        U1 = beta * np.cos(theta)
        U2 = beta * np.sin(theta)

        comp_1 = stats.exponweib.ppf(stats.norm.cdf(U1),a=self.para_dist_1[0],c=self.para_dist_1[1],loc=self.para_dist_1[2],scale=self.para_dist_1[3])

        lamda_cond=self.mean_cond[0]+self.mean_cond[1]*comp_1+self.mean_cond[2]*comp_1**2+self.mean_cond[3]*comp_1**3      # mean of Ln(T) as a function of Hs
        sigma_cond=self.std_cond[0]+self.std_cond[1]*comp_1+self.std_cond[2]*comp_1**2                                # Standard deviation of Ln(T) as a function of Hs

        comp_2_Rosenblatt = stats.lognorm.ppf(stats.norm.cdf(U2),s=sigma_cond,loc=0,scale=np.exp(lamda_cond))  # lognormal inverse

        Hs_Return = comp_1
        T_Return = comp_2_Rosenblatt

        self.Hs_ReturnContours = Hs_Return
        self.T_ReturnContours = T_Return
        return Hs_Return, T_Return

    def getSamples(self):
        '''Currently not implemented in this version'''
        raise NotImplementedError

    def _saveParams(self, groupObj):
        groupObj.create_dataset('n_size', data=self.n_size)
        groupObj.create_dataset('bin_1_limit', data=self.bin_1_limit)
        groupObj.create_dataset('bin_step', data=self.bin_step)
        groupObj.create_dataset('para_dist_1', data=self.para_dist_1)
        groupObj.create_dataset('para_dist_2', data=self.para_dist_2)
        groupObj.create_dataset('mean_cond', data=self.mean_cond)
        groupObj.create_dataset('std_cond', data=self.std_cond)


class ClaytonCopula(EA):
    '''Create a ClaytonCopula EA class for a buoy object. Contours generated 
    under this class will use a Clayton copula.'''    
    def __init__(self, buoy, n_size=40., bin_1_limit=1., bin_step=0.25):
        '''
        Parameters
        ----------
            buoy : NDBCData
                ESSC.Buoy Object
            n_size: float
                minimum bin size used for Copula contour methods
            bin_1_limit: float
                maximum value of Hs for the first bin
            bin_step: float
                overlap interval for each bin
        '''
        self.method = "Clayton Copula"
        self.buoy = buoy
        self.n_size = n_size
        self.bin_1_limit = bin_1_limit
        self.bin_step = bin_step

        self.Hs_ReturnContours = None
#        self.Hs_SampleCA = None
#        self.Hs_SampleFSS = None

        self.T_ReturnContours = None
#        self.T_SampleCA = None
#        self.T_SampleFSS = None

#        self.Weight_points = None

#        self.coeff, self.shift, self.comp1_params, self.sigma_param, self.mu_param = self.__generateParams(size_bin)
        self.para_dist_1,self.para_dist_2,self.mean_cond,self.std_cond = self._EA__getCopulaParams(n_size,bin_1_limit,bin_step)

    def getContours(self, time_ss, time_r, nb_steps = 1000):
        '''WDRT Extreme Sea State Clayton Copula Contour function.
        This function calculates environmental contours of extreme sea states using
        a Clayton copula and the inverse first-order reliability
        method.

        Parameters
        ----------
        time_ss : float
            Sea state duration (hours) of measurements in input.
        time_r : np.array
            Desired return period (years) for calculation of environmental
            contour, can be a scalar or a vector.
        nb_steps : float
            Discretization of the circle in the normal space used for
            inverse FORM calculation.

        Returns
        -------
        Hs_Return : np.array
            Calculated Hs values along the contour boundary following
            return to original input orientation.
        T_Return : np.array
           Calculated T values along the contour boundary following
           return to original input orientation.
        nb_steps : float
            Discretization of the circle in the normal space

        Example
        -------
        To obtain the contours for a NDBC buoy::
            
            import WDRT.ESSC as ESSC
            
            # Pull spectral data from NDBC website
            buoy46022 = ESSC.Buoy('46022','NDBC')
            buoy46022.fetchFromWeb()
            
            # Create Environtmal Analysis object using above parameters
            Clayton46022 = ESSC.ClaytonCopula(buoy46022)
            
            # Declare required parameters
            Time_SS = 1.  # Sea state duration (hrs)
            Time_r = 100  # Return periods (yrs) of interest
            nb_steps = 1000  # Enter discretization of the circle in the normal space (optional)
            
            # Clayton copula contour generation example
            Hs_Return, T_Return = Clayton46022.getContours(Time_SS, Time_r,nb_steps)
        '''
        self.time_ss = time_ss
        self.time_r = time_r
        self.nb_steps = nb_steps

        p_f = 1 / (365 * (24 / time_ss) * time_r)
        beta = stats.norm.ppf((1 - p_f), loc=0, scale=1)  # Reliability
        theta = np.linspace(0, 2 * np.pi, num = nb_steps)
        # Vary U1, U2 along circle sqrt(U1^2+U2^2)=beta
        U1 = beta * np.cos(theta)
        U2 = beta * np.sin(theta)

        comp_1 = stats.exponweib.ppf(stats.norm.cdf(U1),a=self.para_dist_1[0],c=self.para_dist_1[1],loc=self.para_dist_1[2],scale=self.para_dist_1[3])

        tau = stats.kendalltau(self.buoy.T,self.buoy.Hs)[0] # Calculate Kendall's tau
        theta_clay = (2.*tau)/(1.-tau)

        z2_Clay=((1.-stats.norm.cdf(U1)**(-theta_clay)+stats.norm.cdf(U1)**(-theta_clay)/stats.norm.cdf(U2))**(theta_clay/(1.+theta_clay)))**(-1./theta_clay)
        comp_2_Clayton = stats.lognorm.ppf(z2_Clay,s=self.para_dist_2[1],loc=0,scale=np.exp(self.para_dist_2[0])) #lognormalinverse

        Hs_Return = comp_1
        T_Return = comp_2_Clayton

        self.Hs_ReturnContours = Hs_Return
        self.T_ReturnContours = T_Return
        return Hs_Return, T_Return

    def getSamples(self):
        '''Currently not implemented in this version'''
        raise NotImplementedError

    def _saveParams(self, groupObj):
        groupObj.create_dataset('n_size', data=self.n_size)
        groupObj.create_dataset('bin_1_limit', data=self.bin_1_limit)
        groupObj.create_dataset('bin_step', data=self.bin_step)
        groupObj.create_dataset('para_dist_1', data=self.para_dist_1)
        groupObj.create_dataset('para_dist_2', data=self.para_dist_2)
        groupObj.create_dataset('mean_cond', data=self.mean_cond)
        groupObj.create_dataset('std_cond', data=self.std_cond)


class GumbelCopula(EA):
    '''Create a GumbelCopula EA class for a buoy object. Contours generated 
    under this class will use a Gumbel copula.'''    
    def __init__(self, buoy, n_size=40., bin_1_limit=1., bin_step=0.25,Ndata = 1000):
        '''
        Parameters
        ----------
            buoy : NDBCData
                ESSC.Buoy Object
            n_size: float
                minimum bin size used for Copula contour methods
            bin_1_limit: float
                maximum value of Hs for the first bin
            bin_step: float
                overlap interval for each bin
            Ndata: int
                discretization used in the Gumbel copula density estimation, 
                must be less than the number of contour points used in 
                getContours
        '''
        self.method = "Gumbel Copula"
        self.buoy = buoy
        self.n_size = n_size
        self.bin_1_limit = bin_1_limit
        self.bin_step = bin_step

        self.Hs_ReturnContours = None
#        self.Hs_SampleCA = None
#        self.Hs_SampleFSS = None
        self.T_ReturnContours = None
#        self.T_SampleCA = None
#        self.T_SampleFSS = None
#        self.Weight_points = None

#        self.coeff, self.shift, self.comp1_params, self.sigma_param, self.mu_param = self.__generateParams(size_bin)
        self.Ndata = Ndata
        self.min_limit_2 = 0.
        self.max_limit_2 = np.ceil(np.amax(self.buoy.T)*2)
        self.para_dist_1,self.para_dist_2,self.mean_cond,self.std_cond = self._EA__getCopulaParams(n_size,bin_1_limit,bin_step)

    def getContours(self, time_ss, time_r, nb_steps = 1000):
        '''WDRT Extreme Sea State Gumbel Copula Contour function
        This function calculates environmental contours of extreme sea states using
        a Gumbel copula and the inverse first-order reliability
        method.

        Parameters
        ___________
        time_ss : float
            Sea state duration (hours) of measurements in input.
        time_r : np.array
            Desired return period (years) for calculation of environmental
            contour, can be a scalar or a vector.
        nb_steps : float
            Discretization of the circle in the normal space used for
            inverse FORM calculation.

        Returns
        -------
        Hs_Return : np.array
            Calculated Hs values along the contour boundary following
            return to original input orientation.
        T_Return : np.array
           Calculated T values along the contour boundary following
           return to original input orientation.
        nb_steps : float
            Discretization of the circle in the normal space

        Example
        -------
        To obtain the contours for a NDBC buoy::
            
            import WDRT.ESSC as ESSC
            
            # Pull spectral data from NDBC website
            buoy46022 = ESSC.Buoy('46022','NDBC')
            buoy46022.fetchFromWeb()
            
            # Create Environtmal Analysis object using above parameters
            Gumbel46022 = ESSC.GumbelCopula(buoy46022)
            
            # Declare required parameters
            Time_SS = 1.  # Sea state duration (hrs)
            Time_r = 100  # Return periods (yrs) of interest
            nb_steps = 1000  # Enter discretization of the circle in the normal space (optional)
            
            # Gumbel copula contour generation example
            Hs_Return, T_Return = Gumbel46022.getContours(Time_SS,Time_r,nb_steps)
        '''
        self.time_ss = time_ss
        self.time_r = time_r
        self.nb_steps = nb_steps

        p_f = 1 / (365 * (24 / time_ss) * time_r)
        beta = stats.norm.ppf((1 - p_f), loc=0, scale=1)  # Reliability
        theta = np.linspace(0, 2 * np.pi, num = nb_steps)
        # Vary U1, U2 along circle sqrt(U1^2+U2^2)=beta
        U1 = beta * np.cos(theta)
        U2 = beta * np.sin(theta)

        comp_1 = stats.exponweib.ppf(stats.norm.cdf(U1),a=self.para_dist_1[0],c=self.para_dist_1[1],loc=self.para_dist_1[2],scale=self.para_dist_1[3])

        tau = stats.kendalltau(self.buoy.T,self.buoy.Hs)[0] # Calculate Kendall's tau
        theta_gum = 1./(1.-tau)

        fi_u1=stats.norm.cdf(U1);
        fi_u2=stats.norm.cdf(U2);
        x2 = np.linspace(self.min_limit_2,self.max_limit_2,self.Ndata)
        z2 = stats.lognorm.cdf(x2,s=self.para_dist_2[1],loc=0,scale=np.exp(self.para_dist_2[0]))

        comp_2_Gumb = np.zeros(nb_steps)
        for k in range(0,int(nb_steps)):
            z1 = np.linspace(fi_u1[k],fi_u1[k],self.Ndata)
            Z = np.array((z1,z2))
            Y = self._EA__gumbelCopula(Z, theta_gum) # Copula density function
            Y =np.nan_to_num(Y)
            p_x2_x1 = Y*(stats.lognorm.pdf(x2, s = self.para_dist_2[1], loc=0, scale = np.exp(self.para_dist_2[0]))) # pdf 2|1, f(comp_2|comp_1)=c(z1,z2)*f(comp_2)
            dum = np.cumsum(p_x2_x1)
            cdf = dum/(dum[self.Ndata-1]) # Estimate CDF from PDF
            table = np.array((x2, cdf)) # Result of conditional CDF derived based on Gumbel copula
            table = table.T
            for j in range(self.Ndata):
                if fi_u2[k] <= table[0,1]:
                    comp_2_Gumb[k] = min(table[:,0])
                    break
                elif fi_u2[k] <= table[j,1]:
                    comp_2_Gumb[k] = (table[j,0]+table[j-1,0])/2
                    break
                else:
                    comp_2_Gumb[k] = table[:,0].max()

        Hs_Return = comp_1
        T_Return = comp_2_Gumb

        self.Hs_ReturnContours = Hs_Return
        self.T_ReturnContours = T_Return
        return Hs_Return, T_Return

    def getSamples(self):
        '''Currently not implemented in this version.'''
        raise NotImplementedError

    def _saveParams(self, groupObj):
        groupObj.create_dataset('Ndata', data=self.Ndata)
        groupObj.create_dataset('min_limit_2', data=self.min_limit_2)
        groupObj.create_dataset('max_limit_2', data=self.max_limit_2)
        groupObj.create_dataset('n_size', data=self.n_size)
        groupObj.create_dataset('bin_1_limit', data=self.bin_1_limit)
        groupObj.create_dataset('bin_step', data=self.bin_step)
        groupObj.create_dataset('para_dist_1', data=self.para_dist_1)
        groupObj.create_dataset('para_dist_2', data=self.para_dist_2)
        groupObj.create_dataset('mean_cond', data=self.mean_cond)
        groupObj.create_dataset('std_cond', data=self.std_cond)

class NonParaGaussianCopula(EA):
    '''Create a NonParaGaussianCopula EA class for a buoy object. Contours
    generated under this class will use a Gaussian copula with non-parametric
    marginal distribution fits.'''
    def __init__(self, buoy, Ndata = 1000, max_T=None, max_Hs=None):
        '''
        Parameters
        ----------
            buoy : NDBCData
                ESSC.Buoy Object
            NData: int
                discretization resolution used in KDE construction
            max_T:float
                Maximum T value for KDE contstruction, must include possible 
                range of contour. Default value is 2*max(T)
            max_Hs:float
                Maximum Hs value for KDE contstruction, must include possible 
                range of contour. Default value is 2*max(Hs)    
        '''
        self.method = "Non-parametric Gaussian Copula"
        self.buoy = buoy
        self.Ndata = Ndata

        self.Hs_ReturnContours = None
#        self.Hs_SampleCA = None
#        self.Hs_SampleFSS = None

        self.T_ReturnContours = None
#        self.T_SampleCA = None
#        self.T_SampleFSS = None

#        self.Weight_points = None

#        self.coeff, self.shift, self.comp1_params, self.sigma_param, self.mu_param = self.__generateParams(size_bin)
        if max_T == None:
            max_T = max(self.buoy.T)*2.
        if max_Hs == None:
            max_Hs = max(self.buoy.Hs)*2.
        
        self.max_T = max_T
        self.max_Hs = max_Hs
        self.nonpara_dist_1,self.nonpara_dist_2,self.nonpara_pdf_2 = self._EA__getNonParaCopulaParams(Ndata,max_T,max_Hs)

    def getContours(self, time_ss, time_r, nb_steps = 1000):
        '''WDRT Extreme Sea State Gaussian Copula Contour function.
        This function calculates environmental contours of extreme sea states 
        using a Gaussian copula with non-parametric marginal distribution fits
        and the inverse first-order reliability method.

        Parameters
        ___________
        time_ss : float
            Sea state duration (hours) of measurements in input.
        time_r : np.array
            Desired return period (years) for calculation of environmental
            contour, can be a scalar or a vector.
        nb_steps : float
            Discretization of the circle in the normal space used for
            inverse FORM calculation.

        Returns
        -------
        Hs_Return : np.array
            Calculated Hs values along the contour boundary following
            return to original input orientation.
        T_Return : np.array
           Calculated T values along the contour boundary following
           return to original input orientation.
        nb_steps : float
            Discretization of the circle in the normal space

        Example
        -------
        To obtain the contours for a NDBC buoy::
            
            import WDRT.ESSC as ESSC
            
            # Pull spectral data from NDBC website
            buoy46022 = ESSC.Buoy('46022','NDBC')
            buoy46022.fetchFromWeb()
            
            # Create Environtmal Analysis object using above parameters
            NonParaGauss46022 = ESSC.NonParaGaussianCopula(buoy46022)
            
            # Declare required parameters
            Time_SS = 1.  # Sea state duration (hrs)
            Time_r = 100  # Return periods (yrs) of interest
            nb_steps = 1000  # Enter discretization of the circle in the normal space (optional)
            
            # Non-Parametric Gaussian copula contour generation example
            Hs_Return, T_Return = NonParaGauss46022.getContours(Time_SS, Time_r,nb_steps)
        '''
        self.time_ss = time_ss
        self.time_r = time_r
        self.nb_steps = nb_steps
        
        comp_1 = np.zeros(nb_steps)
        comp_2_Gau = np.zeros(nb_steps)
        
        # Inverse FORM
        p_f = 1 / (365 * (24 / time_ss) * time_r)
        beta = stats.norm.ppf((1 - p_f), loc=0, scale=1)  # Reliability
        
        # Normal Space
        theta = np.linspace(0, 2 * np.pi, num = nb_steps)
        U1 = beta * np.cos(theta)
        U2 = beta * np.sin(theta)
        
        # Copula parameters
        tau = stats.kendalltau(self.buoy.T,self.buoy.Hs)[0]# Calculate Kendall's tau    
        rho_gau=np.sin(tau*np.pi/2.);
        
        # Component 1 (Hs)
        z1_Hs = stats.norm.cdf(U1)
        for k in range(0,nb_steps):
            for j in range(0,np.size(self.nonpara_dist_1,0)):
                if z1_Hs[k] <= self.nonpara_dist_1[0,1]: 
                    comp_1[k] = min(self.nonpara_dist_1[:,0]) 
                    break
                elif z1_Hs[k] <= self.nonpara_dist_1[j,1]: 
                    comp_1[k] = (self.nonpara_dist_1[j,0] + self.nonpara_dist_1[j-1,0])/2
                    break
                else:
                    comp_1[k]= max(self.nonpara_dist_1[:,0])
        
        # Component 2 (T)
        
        z2_Gau=stats.norm.cdf(U2*np.sqrt(1.-rho_gau**2.)+rho_gau*U1);
        for k in range(0,nb_steps):
            for j in range(0,np.size(self.nonpara_dist_2,0)):
                if z2_Gau[k] <= self.nonpara_dist_2[0,1]: 
                    comp_2_Gau[k] = min(self.nonpara_dist_2[:,0]) 
                    break
                elif z2_Gau[k] <= self.nonpara_dist_2[j,1]: 
                    comp_2_Gau[k] = (self.nonpara_dist_2[j,0] + self.nonpara_dist_2[j-1,0])/2
                    break
                else:
                    comp_2_Gau[k]= max(self.nonpara_dist_2[:,0])

        Hs_Return = comp_1
        T_Return = comp_2_Gau

        self.Hs_ReturnContours = Hs_Return
        self.T_ReturnContours = T_Return
        return Hs_Return, T_Return

    def getSamples(self):
        '''Currently not implemented in this version.'''
        raise NotImplementedError

    def _saveParams(self, groupObj):
        groupObj.create_dataset('nonpara_dist_1', data=self.nonpara_dist_1)
        groupObj.create_dataset('nonpara_dist_2', data=self.nonpara_dist_2)
        
class NonParaClaytonCopula(EA):
    '''Create a NonParaClaytonCopula EA class for a buoy object. Contours
    generated under this class will use a Clayton copula with non-parametric
    marginal distribution fits.'''
    def __init__(self, buoy, Ndata = 1000, max_T=None, max_Hs=None):
        '''
        Parameters
        ----------
            buoy : NDBCData
                ESSC.Buoy Object
            NData: int
                discretization resolution used in KDE construction
            max_T:float
                Maximum T value for KDE contstruction, must include possible 
                range of contour. Default value is 2*max(T)
            max_Hs:float
                Maximum Hs value for KDE contstruction, must include possible 
                range of contour. Default value is 2*max(Hs)    
        '''
        self.method = "Non-parametric Clayton Copula"
        self.buoy = buoy
        self.Ndata = Ndata

        self.Hs_ReturnContours = None
#        self.Hs_SampleCA = None
#        self.Hs_SampleFSS = None

        self.T_ReturnContours = None
#        self.T_SampleCA = None
#        self.T_SampleFSS = None

#        self.Weight_points = None

#        self.coeff, self.shift, self.comp1_params, self.sigma_param, self.mu_param = self.__generateParams(size_bin)
        if max_T == None:
            max_T = max(self.buoy.T)*2.
        if max_Hs == None:
            max_Hs = max(self.buoy.Hs)*2.
        
        self.max_T = max_T
        self.max_Hs = max_Hs
        self.nonpara_dist_1,self.nonpara_dist_2,self.nonpara_pdf_2 = self._EA__getNonParaCopulaParams(Ndata,max_T,max_Hs)

    def getContours(self, time_ss, time_r, nb_steps = 1000):
        '''WDRT Extreme Sea State non-parameteric Clayton Copula Contour
        function. This function calculates environmental contours of extreme
        sea states using a Clayton copula with non-parametric marginal
        distribution fits and the inverse first-order reliability method.

        Parameters
        ___________
        time_ss : float
            Sea state duration (hours) of measurements in input.
        time_r : np.array
            Desired return period (years) for calculation of environmental
            contour, can be a scalar or a vector.
        nb_steps : float
            Discretization of the circle in the normal space used for
            inverse FORM calculation.

        Returns
        -------
        Hs_Return : np.array
            Calculated Hs values along the contour boundary following
            return to original input orientation.
        T_Return : np.array
           Calculated T values along the contour boundary following
           return to original input orientation.
        nb_steps : float
            Discretization of the circle in the normal space

        Example
        -------
        To obtain the contours for a NDBC buoy::
            
            import WDRT.ESSC as ESSC
            
            # Pull spectral data from NDBC website
            buoy46022 = ESSC.Buoy('46022','NDBC')
            buoy46022.fetchFromWeb()
            
            # Create Environtmal Analysis object using above parameters
            NonParaClayton46022 = ESSC.NonParaClaytonCopula(buoy46022)
            
            # Declare required parameters
            Time_SS = 1.  # Sea state duration (hrs)
            Time_r = 100  # Return periods (yrs) of interest
            nb_steps = 1000  # Enter discretization of the circle in the normal space (optional)
            
            # Non-Parametric Clayton copula contour generation example
            Hs_Return, T_Return = NonParaClayton46022.getContours(Time_SS, Time_r,nb_steps)
        '''
        self.time_ss = time_ss
        self.time_r = time_r
        self.nb_steps = nb_steps
        
        comp_1 = np.zeros(nb_steps)
        comp_2_Clay = np.zeros(nb_steps)
        
        # Inverse FORM
        p_f = 1 / (365 * (24 / time_ss) * time_r)
        beta = stats.norm.ppf((1 - p_f), loc=0, scale=1)  # Reliability
        
        # Normal Space
        theta = np.linspace(0, 2 * np.pi, num = nb_steps)
        U1 = beta * np.cos(theta)
        U2 = beta * np.sin(theta)
        
        # Copula parameters
        tau = stats.kendalltau(self.buoy.T,self.buoy.Hs)[0]# Calculate Kendall's tau    
        theta_clay = (2.*tau)/(1.-tau);
        
        # Component 1 (Hs)
        z1_Hs = stats.norm.cdf(U1)
        for k in range(0,nb_steps):
            for j in range(0,np.size(self.nonpara_dist_1,0)):
                if z1_Hs[k] <= self.nonpara_dist_1[0,1]: 
                    comp_1[k] = min(self.nonpara_dist_1[:,0]) 
                    break
                elif z1_Hs[k] <= self.nonpara_dist_1[j,1]: 
                    comp_1[k] = (self.nonpara_dist_1[j,0] + self.nonpara_dist_1[j-1,0])/2
                    break
                else:
                    comp_1[k]= max(self.nonpara_dist_1[:,0])
        
        # Component 2 (T)
        
        z2_Clay=((1.-stats.norm.cdf(U1)**(-theta_clay)+stats.norm.cdf(U1)**(-theta_clay)/stats.norm.cdf(U2))**(theta_clay/(1.+theta_clay)))**(-1./theta_clay)
        for k in range(0,nb_steps):
            for j in range(0,np.size(self.nonpara_dist_2,0)):
                if z2_Clay[k] <= self.nonpara_dist_2[0,1]: 
                    comp_2_Clay[k,0] = min(self.nonpara_dist_2[:,0]) 
                    break
                elif z2_Clay[k] <= self.nonpara_dist_2[j,1]: 
                    comp_2_Clay[k] = (self.nonpara_dist_2[j,0] + self.nonpara_dist_2[j-1,0])/2
                    break
                else:
                    comp_2_Clay[k]= max(self.nonpara_dist_2[:,0])

        Hs_Return = comp_1
        T_Return = comp_2_Clay

        self.Hs_ReturnContours = Hs_Return
        self.T_ReturnContours = T_Return
        return Hs_Return, T_Return

    def getSamples(self):
        '''Currently not implemented in this version.'''
        raise NotImplementedError

    def _saveParams(self, groupObj):
        groupObj.create_dataset('nonpara_dist_1', data=self.nonpara_dist_1)
        groupObj.create_dataset('nonpara_dist_2', data=self.nonpara_dist_2)

class NonParaGumbelCopula(EA):
    '''Create a NonParaGumbelCopula EA class for a buoy object. Contours
    generated under this class will use a Gumbel copula with non-parametric
    marginal distribution fits.'''
    def __init__(self, buoy, Ndata = 1000, max_T=None, max_Hs=None):
        '''
        Parameters
        ----------
            buoy : NDBCData
                ESSC.Buoy Object
            NData: int
                discretization resolution used in KDE construction
            max_T:float
                Maximum T value for KDE contstruction, must include possible 
                range of contour. Default value is 2*max(T)
            max_Hs:float
                Maximum Hs value for KDE contstruction, must include possible 
                range of contour. Default value is 2*max(Hs)    
        '''
        self.method = "Non-parametric Gumbel Copula"
        self.buoy = buoy
        self.Ndata = Ndata

        self.Hs_ReturnContours = None
#        self.Hs_SampleCA = None
#        self.Hs_SampleFSS = None

        self.T_ReturnContours = None
#        self.T_SampleCA = None
#        self.T_SampleFSS = None

#        self.Weight_points = None

#        self.coeff, self.shift, self.comp1_params, self.sigma_param, self.mu_param = self.__generateParams(size_bin)
        if max_T == None:
            max_T = max(self.buoy.T)*2.
        if max_Hs == None:
            max_Hs = max(self.buoy.Hs)*2.
        
        self.max_T = max_T
        self.max_Hs = max_Hs
        self.nonpara_dist_1,self.nonpara_dist_2,self.nonpara_pdf_2 = self._EA__getNonParaCopulaParams(Ndata,max_T,max_Hs)

    def getContours(self, time_ss, time_r, nb_steps = 1000):
        '''WDRT Extreme Sea State non-parameteric Gumbel Copula Contour
        function. This function calculates environmental contours of extreme
        sea states using a Gumbel copula with non-parametric marginal
        distribution fits and the inverse first-order reliability method.

        Parameters
        ___________
        time_ss : float
            Sea state duration (hours) of measurements in input.
        time_r : np.array
            Desired return period (years) for calculation of environmental
            contour, can be a scalar or a vector.
        nb_steps : float
            Discretization of the circle in the normal space used for
            inverse FORM calculation.

        Returns
        -------
        Hs_Return : np.array
            Calculated Hs values along the contour boundary following
            return to original input orientation.
        T_Return : np.array
           Calculated T values along the contour boundary following
           return to original input orientation.
        nb_steps : float
            Discretization of the circle in the normal space

        Example
        -------
        To obtain the contours for a NDBC buoy::
            
            import WDRT.ESSC as ESSC
            
            # Pull spectral data from NDBC website
            buoy46022 = ESSC.Buoy('46022','NDBC')
            buoy46022.fetchFromWeb()
            
            # Create Environtmal Analysis object using above parameters
            NonParaGumbel46022 = ESSC.NonParaGumbelCopula(buoy46022)
            
            # Declare required parameters
            Time_SS = 1.  # Sea state duration (hrs)
            Time_r = 100  # Return periods (yrs) of interest
            nb_steps = 1000  # Enter discretization of the circle in the normal space (optional)
            
            # Non-Parametric Gumbel copula contour generation example
            Hs_Return, T_Return = NonParaGumbel46022.getContours(Time_SS, Time_r,nb_steps)
        '''
        self.time_ss = time_ss
        self.time_r = time_r
        self.nb_steps = nb_steps
        
        comp_1 = np.zeros(nb_steps)
        comp_2_Gumb = np.zeros(nb_steps)
        
        # Inverse FORM
        p_f = 1 / (365 * (24 / time_ss) * time_r)
        beta = stats.norm.ppf((1 - p_f), loc=0, scale=1)  # Reliability
        
        # Normal Space
        theta = np.linspace(0, 2 * np.pi, num = nb_steps)
        U1 = beta * np.cos(theta)
        U2 = beta * np.sin(theta)
        
        # Copula parameters
        tau = stats.kendalltau(self.buoy.T,self.buoy.Hs)[0]# Calculate Kendall's tau    
        theta_gum = 1./(1.-tau);
        
        # Component 1 (Hs)
        z1_Hs = stats.norm.cdf(U1)
        for k in range(0,nb_steps):
            for j in range(0,np.size(self.nonpara_dist_1,0)):
                if z1_Hs[k] <= self.nonpara_dist_1[0,1]: 
                    comp_1[k] = min(self.nonpara_dist_1[:,0]) 
                    break
                elif z1_Hs[k] <= self.nonpara_dist_1[j,1]: 
                    comp_1[k] = (self.nonpara_dist_1[j,0] + self.nonpara_dist_1[j-1,0])/2
                    break
                else:
                    comp_1[k]= max(self.nonpara_dist_1[:,0])
        
        # Component 2 (T)
        
        fi_u1=stats.norm.cdf(U1);
        fi_u2=stats.norm.cdf(U2);
        
        for k in range(0,nb_steps):
            z1 = np.linspace(fi_u1[k],fi_u1[k],self.Ndata)
            Z = np.array((np.transpose(z1),self.nonpara_dist_2[:,1]))
            Y = self._EA__gumbelCopula(Z, theta_gum)
            Y =np.nan_to_num(Y) # Need to look into this
            p_x2_x1 = Y*self.nonpara_pdf_2[:,1]
            dum = np.cumsum(p_x2_x1)
            cdf = dum/(dum[self.Ndata-1])
            table = np.array((self.nonpara_pdf_2[:,0], cdf))
            table = table.T
            for j in range(self.Ndata):
                if fi_u2[k] <= table[0,1]:
                    comp_2_Gumb[k] = min(table[:,0])
                    break
                elif fi_u2[k] <= table[j,1]:
                    comp_2_Gumb[k] = (table[j,0]+table[j-1,0])/2
                    break
                else: 
                    comp_2_Gumb[k] = max(table[:,0])


        Hs_Return = comp_1
        T_Return = comp_2_Gumb

        self.Hs_ReturnContours = Hs_Return
        self.T_ReturnContours = T_Return
        return Hs_Return, T_Return

    def getSamples(self):
        '''Currently not implemented in this version.'''
        raise NotImplementedError

    def _saveParams(self, groupObj):
        groupObj.create_dataset('nonpara_dist_1', data=self.nonpara_dist_1)
        groupObj.create_dataset('nonpara_dist_2', data=self.nonpara_dist_2)

class BivariateKDE(EA):
    '''Create a BivariateKDE EA class for a buoy object. Contours
    generated under this class will use a non-parametric KDE to fit the joint distribution.'''
    def __init__(self, buoy, bw, NData = 100, logTransform = False, max_T=None, max_Hs=None):
        '''
        Parameters
        ----------
            buoy : NDBCData
                ESSC.Buoy Object
            bw: np.array
                Array containing KDE bandwidth for Hs and T
            NData: int
                Discretization resolution used in KDE construction
            logTransform: Boolean
                Logical. True if log transformation should be taken prior to 
                KDE construction. Default value is False. 
            max_T:float
                Maximum T value for KDE contstruction, must include possible 
                range of contour. Default value is 2*max(T)
            max_Hs:float
                Maximum Hs value for KDE contstruction, must include possible 
                range of contour. Default value is 2*max(Hs)    
        '''
        
        if logTransform:        
            self.method = "Bivariate KDE, Log Transform"
        else:
            self.method = "Bivariate KDE"
        self.buoy = buoy
        
        if max_T == None:
            max_T = max(self.buoy.T)*2.
        if max_Hs == None:
            max_Hs = max(self.buoy.Hs)*2.

        self.max_T = max_T
        self.max_Hs = max_Hs
        self.Hs_ReturnContours = None
        self.T_ReturnContours = None
        self.NData = NData
        self.bw = bw
        self.logTransform = logTransform

    def getContours(self, time_ss, time_r):
        '''WDRT Extreme Sea State non-parameteric bivariate KDE Contour
        function. This function calculates environmental contours of extreme
        sea states using a bivariate KDE to estimate the joint distribution. 
        The contour is then calculcated directly from the joint distribution.

        Parameters
        ___________
        time_ss : float
            Sea state duration (hours) of measurements in input.
        time_r : np.array
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
            
            import WDRT.ESSC as ESSC
            
            # Pull spectral data from NDBC website
            buoy46022 = ESSC.Buoy('46022','NDBC')
            buoy46022.fetchFromWeb()
            
            # Create Environmental Analysis object using above parameters
            BivariateKDE46022 = ESSC.BivariateKDE(buoy46022, bw = [0.23,0.19], logTransform = False)
            
            # Declare required parameters
            Time_SS = 1.  # Sea state duration (hrs)
            Time_r = 100  # Return periods (yrs) of interest
            
            # KDE contour generation example
            Hs_Return, T_Return = BivariateKDE46022.getContours(Time_SS, Time_r)
        '''       
        p_f = 1 / (365 * (24 / time_ss) * time_r)

        if self.logTransform: 
            # Take log of both variables
            logTp = np.log(self.buoy.T)
            logHs = np.log(self.buoy.Hs)
            ty = [logTp, logHs]
        else: 
            ty = [self.buoy.T, self.buoy.Hs]
      

        # Create grid of points
        Ndata = self.NData
        min_limit_1 = 0.01
        max_limit_1 = self.max_T
        min_limit_2 = 0.01
        max_limit_2 = self.max_Hs
        pts_tp = np.linspace(min_limit_1, max_limit_1, Ndata) 
        pts_hs = np.linspace(min_limit_2, max_limit_2, Ndata)
        pt1,pt2 = np.meshgrid(pts_tp, pts_hs)
        pts_tp = pt1.flatten()
        pts_hs = pt2.flatten()

        # Transform gridded points using log
        xi = [pts_tp, pts_hs]
        if self.logTransform: 
            txi = [np.log(pts_tp), np.log(pts_hs)]
        else: 
            txi = xi

        m = len(txi[0])
        n = len(ty[0])
        d = 2

        # Create contour
        f = np.zeros((1,m))
        weight = np.ones((1,n))
        for i in range(0,m):
            ftemp = np.ones((n,1))
            for j in range(0,d):
                z = (txi[j][i] - ty[j])/self.bw[j]
                fk = stats.norm.pdf(z)
                if self.logTransform:     
                    fnew = fk*(1/np.transpose(xi[j][i]))
                else: 
                    fnew = fk
                fnew = np.reshape(fnew, (n,1))
                ftemp = np.multiply(ftemp,fnew)
            f[:,i] = np.dot(weight,ftemp)


        fhat = f.reshape(100,100)
        vals = plt.contour(pt1,pt2,fhat, levels = [p_f])
        plt.clf()
        self.Hs_ReturnContours = []
        self.T_ReturnContours = []
        for i,seg in enumerate(vals.allsegs[0]):
            self.Hs_ReturnContours.append(seg[:,1])
            self.T_ReturnContours.append(seg[:,0])
        
        self.Hs_ReturnContours = np.transpose(np.asarray(self.Hs_ReturnContours)[0])
        self.T_ReturnContours = np.transpose(np.asarray(self.T_ReturnContours)[0])
            
#        contourVals = np.empty((0,2))
#        for seg in vals.allsegs[0]:
#            contourVals = np.append(contourVals,seg, axis = 0)
#        self.Hs_ReturnContours = contourVals[:,1]
#        self.T_ReturnContours = contourVals[:,0]

        return self.Hs_ReturnContours, self.T_ReturnContours

        
class Buoy(object):
    '''
    This class creates a buoy object to store buoy data for use in the 
    environmental assessment functions available in the ESSC module.
    
    Attributes
    __________
    swdList : list
        List that contains numpy arrays of the spectral wave density data,
        separated by year.
    freqList: list
        List that contains numpy arrays that contain the frequency values
        for each year
    dateList : list
        List that contains numpy arrays of the date values for each line of
        spectral data, separated by year
    Hs : list
        Significant wave height.
    T : list
        Energy period.
    dateNum : list
        List of datetime objects.
    '''



    def __init__(self, buoyNum, buoyType):

        '''
        Parameters
        ___________
            buoyNum : string
                device number for desired buoy
            buoyType : string
                type of buoy device, available options are 'NDBC' or 'CDIP'
            savePath : string
                relative path where the data read from ndbc.noaa.gov will be stored


        '''
        self.swdList = []
        self.freqList = []
        self.dateList = []
        self.Hs = []
        self.T = []
        self.dateNum = []

        self.buoyNum = buoyNum
        self.buoyType = buoyType.upper()



    def fetchFromWeb(self, savePath = "./Data/",proxy=None):
        '''
        Calls either __fetchCDIP() or __fetchNDBC() depending on the given
        buoy's type and fetches the necessary data from its respective website.

        Parameters
        ----------
        saveType: string
            If set to to "h5", the data will be saved in a compressed .h5
            file
            If set to "txt", the data will be stored in a raw .txt file
            Otherwise, a file will not be created
            NOTE: Only applies 
        savePath : string
            Relative path to place directory with data files.
        proxy: dict
            Proxy server and port, i.e., {http":"http://proxyserver:port"}
        Example
        _________
        >>> import WDRT.ESSC as ESSC
        >>> buoy = ESSC.Buoy('46022','NDBC')
        >>> buoy.fetchFromWeb()
        '''
        if self.buoyType == "NDBC":
            self.__fetchNDBC(proxy)
        elif self.buoyType == "CDIP":
            self.__fetchCDIP(savePath,proxy)

    def __fetchNDBC(self, proxy):
        '''
        Searches ndbc.noaa.gov for the historical spectral wave density
        data of a given device and writes the annual files from the website
        to a single .txt file, and stores the values in the swdList, freqList,
        and dateList member variables.

        Parameters
        ----------
        saveType: string
            If set to to "h5", the data will be saved in a compressed .h5
            file
            If set to "txt", the data will be stored in a raw .txt file
            Otherwise, a file will not be created
            NOTE: Only applies 
        savePath : string
            Relative path to place directory with data files.
        '''
        maxRecordedDateValues = 4
        #preallocates data
        numLines = 0
        numCols = 0
        numDates = 0
        dateVals = []
        spectralVals = []
        #prepares to pull the data from the NDBC website
        url = "https://www.ndbc.noaa.gov/station_history.php?station=%s" % (self.buoyNum)
        if proxy == None:
            ndbcURL = requests.get(url)
        else:
            ndbcURL = requests.get(url,proxies=proxy)
        ndbcURL.raise_for_status()
        ndbcHTML = bs4.BeautifulSoup(ndbcURL.text, "lxml")
        headers = ndbcHTML.findAll("b", text="Spectral wave density data: ")
        #checks for headers in differently formatted webpages
        if len(headers) == 0:
            raise Exception("Spectral wave density data for buoy #%s not found" % self.buoyNum)


        if len(headers) == 2:
            headers = headers[1]
        else:
            headers = headers[0]

        links = [a["href"] for a in headers.find_next_siblings("a", href=True)]
        #downloads files
        for link in links:
            dataLink = "https://ndbc.noaa.gov" + link

            fileName = dataLink.replace('download_data', 'view_text_file')
            data = urllib.request.urlopen(fileName)
            print("Reading from:", data.geturl())

            #First Line of every file contains the frequency data
            frequency = data.readline()
            if frequency.split()[4] == b'mm':
                numDates = 5

            else:
                numDates = 4

            frequency = np.array(frequency.split()[numDates:], dtype = np.float)

            #splits and organizes data into arrays.
            for line in data:
                currentLine = line.split()
                numCols = len(currentLine)
                if numCols - numDates != len(frequency):
                    print("NDBC File is corrupted - Skipping and deleting data")
                    spectralVals = []
                    dateVals = []
                    break

                if float(currentLine[numDates+1]) < 999:
                    numLines += 1
                    for j in range(maxRecordedDateValues):
                        dateVals.append(currentLine[j])
                    for j in range(numCols - numDates):
                        spectralVals.append(currentLine[j + numDates])

            if len(spectralVals) != 0:
                dateValues = np.array(dateVals, dtype=np.int)
                spectralValues = np.array(spectralVals, dtype=np.float)

                dateValues = np.reshape(dateValues, (numLines, maxRecordedDateValues))
                spectralValues = np.reshape(spectralValues, (numLines,
                                                             (numCols - numDates)))
            numLines = 0
            numCols = 0

            if len(spectralVals) != 0:
                del dateVals[:]
                del spectralVals[:]
                self.swdList.append(spectralValues)
                self.freqList.append(frequency)
                self.dateList.append(dateValues)
        self._prepData()


    def loadFromTxt(self, dirPath = None):
        '''Loads NDBC data previously downloaded to a series of text files in the
        specified directory.

        Parameters
        ----------
            dirPath : string
                Relative path to directory containing NDBC text files (created by
                NBDCdata.fetchFromWeb). If left blank, the method will search
                all directories for the data using the current directory as
                the root.


        Example
        -------
        To load data from previously downloaded files 
        created using fetchFromWeb():

            import WDRT.ESSC as ESSC
            buoy46022 = ESSC.Buoy('46022','NDBC')
            buoy46022.loadFromText()
        '''
        #preallocates arrays
        dateVals = []
        spectralVals = []
        numLines = 0
        maxRecordedDateValues = 4
        #finds the text files (if they exist on the machine)
        if dirPath is None:
            for dirpath, subdirs, files in os.walk('.'):
                for dirs in subdirs:
                    if ("NDBC%s" % self.buoyNum) in dirs:
                        dirPath = os.path.join(dirpath,dirs)
                        break
        if dirPath is None:
            raise IOError("Could not find directory containing NDBC data")

        fileList = glob.glob(os.path.join(dirPath,'SWD*.txt'))

        if len(fileList) == 0:
            raise IOError("No NDBC data files found in " + dirPath)
        #reads in the files
        for fileName in fileList:
            print('Reading from: %s' % (fileName))
            f = open(fileName, 'r')
            frequency = f.readline().split()
            numCols = len(frequency)

            if frequency[4] == 'mm':
                frequency = np.array(frequency[5:], dtype=np.float)
                numTimeVals = 5

            else:
                frequency = np.array(frequency[4:], dtype=np.float)
                numTimeVals = 4

            for line in f:
                currentLine = line.split()
                if float(currentLine[numTimeVals + 1]) < 999:
                    numLines += 1
                    for i in range(maxRecordedDateValues):
                        dateVals.append(currentLine[i])
                    for i in range(numCols - numTimeVals):
                        spectralVals.append(currentLine[i + numTimeVals])

            dateValues = np.array(dateVals, dtype=np.int)
            spectralValues = np.array(spectralVals, dtype=np.double)
            dateValues = np.reshape(dateValues, (numLines, maxRecordedDateValues))
            spectralValues = np.reshape(
                spectralValues, (numLines, (numCols - numTimeVals)))

            del dateVals[:]
            del spectralVals[:]

            numLines = 0
            numCols = 0
            self.swdList.append(spectralValues)
            self.freqList.append(frequency)
            self.dateList.append(dateValues)
        self._prepData()

    def loadFile(self, dirPath = None):
        '''Loads file depending on whether it's NDBC or CDIP.'''
        if self.buoyType == "NDBC":
            self.loadFromText(dirPath)
        if self.buoyType == "CDIP":
            self.loadCDIP(dirPath)

    def loadFromH5(self, fileName = None):
        """
        Loads NDBC data previously saved in a .h5 file

        Parameters
        ----------
            fileName : string
                Name of the .h5 file to load data from.
        Example
        -------
        To load data from previously downloaded files:

            import WDRT.ESSC as ESSC
            buoy46022 = ESSC.Buoy('46022','NDBC')
            buoy46022.fetchFromWeb()
            buoy46022.saveData()
            buoy46022.loadFromH5('NDBC46022.h5')
        """
        if fileName == None:
            fileName = self.buoyType + self.buoyNum + ".h5"
        _, file_extension = os.path.splitext(fileName)
        if not file_extension:
            fileName = fileName + '.h5'
        print("Reading from: ", fileName)
        try:
            f = h5py.File(fileName, 'r')
        except IOError:
            raise IOError("Could not find file: " + fileName)
        self.Hs = np.array(f['buoy_Data/Hs'][:])
        self.T = np.array(f['buoy_Data/Te'][:])
        self.dateNum = np.array(f['buoy_Data/dateNum'][:])
        self.dateList = np.array(f['buoy_Data/dateList'][:])
        print("----> SUCCESS")

    def saveAsH5(self, fileName=None):
        '''
        Saves NDBC buoy data to h5 file after fetchFromWeb() or loadFromText(). 
        This data can later be used to create a buoy object using the 
        loadFromH5() function.

        Parameters
        ----------
            fileName : string
                relevent path and filename where the .h5 file will be created
                and saved. If no filename, the h5 file will be named 
                NDBC(buoyNum).h5 in location where code is running.
                
        Example
        -------
        To save data to h5 file after fetchFromWeb or loadFromText:
            
            import WDRT.ESSC as ESSC
            buoy46022 = ESSC.Buoy('46022','NDBC')
            buoy46022.fetchFromWeb()
            buoy46022.saveAsH5()
        '''
        if (fileName == None):
            fileName = 'NDBC' + str(self.buoyNum) + '.h5'
        else:
            _, file_extension = os.path.splitext(fileName)
            if not file_extension:
                fileName = fileName + '.h5'
        f = h5py.File(fileName, 'w')
        self._saveData(f)
        f.close()
        print("Saved buoy data");
     
    def saveAsTxt(self, savePath = "./Data/"):
        """
        Saves spectral wave density data to a .txt file in the same format as the files 
        found on NDBC's website.

        Parameters
        ----------
            savePath : string
                Relative file path where the .txt files will be saved. 
                
        Example
        -------
        To save data to h5 file after fetchFromWeb or loadFromText:
            
            import WDRT.ESSC as ESSC
            buoy46022 = ESSC.Buoy('46022','NDBC')
            buoy46022.fetchFromWeb()
            buoy46022.saveAsTxt()
        """
        curYear = self.dateList[0][0]
        dateIndexDiff = 0
        bFile = False #NDBC sometimes splits years into two files, the second one titled "YYYYb"
        saveDir = os.path.join(savePath, 'NDBC%s' % (self.buoyNum))
        print("Saving in :", saveDir)
        if not os.path.exists(saveDir):
            os.makedirs(saveDir)
        for i in range(len(self.swdList)):
            if not bFile:
                swdFile = open(os.path.join(saveDir, "SWD-%s-%d.txt" %
                                   (self.buoyNum, curYear)), 'w')
            else:
                swdFile = open(os.path.join(saveDir, "SWD-%s-%db.txt" %
                                   (self.buoyNum, curYear)), 'w')
                bFile = False
            
            freqLine = "YYYY MM DD hh"
            for j in range(len(self.freqList[i])):
                freqLine += ("   " + "%2.4f" % self.freqList[i][j])
            freqLine += "\n"
            swdFile.write(freqLine)
            for j in range(len(self.dateList)):
                if (j + dateIndexDiff + 1) > len(self.dateList):
                    break

                newYear = self.dateList[j + dateIndexDiff][0]
                if curYear != newYear:
                    dateIndexDiff += (j)
                    curYear = newYear
                    break
                if (j+1) > len(self.swdList[i]):
                    dateIndexDiff += (j)
                    bFile = True
                    break            
                swdLine = ' '.join("%0*d" % (2,dateVal) for dateVal in self.dateList[j + dateIndexDiff]) + "   "
                swdLine += "   ".join("%6s" % val for val in self.swdList[i][j]) + "\n"
                swdFile.write(swdLine)
                

        


    def createSubsetBuoy(self, trainingSize):
        
        '''Takes a given buoy and creates a subset buoy of a given length in years.
        
        Parameters
        ----------
            trainingSize : int
                The size in years of the subset buoy you would like to create
            
        Returns
        -------
            # subsetBuoy : ESSC.Buoy object
                A buoy (with Hs, T, and dateList values) that is a subset of the given buoy 
        
        Example
        -------
        
        To get a corresponding subset of a buoy with a given number of years:
            
            import WDRT.ESSC as ESSC
            
            # Pull spectral data from NDBC website
            buoy46022 = ESSC.Buoy('46022','NDBC')
            buoy46022.fetchFromWeb()
            
            # Create a subset of buoy 46022 consisting of the first 10 years
            subsetBuoy = buoy46022.createSubsetBuoy(10)
        
        
        '''  
        
        subsetBuoy = copy.deepcopy(self)
    
        sortedIndex = sorted(range(len(self.dateNum)),key=lambda x:self.dateNum[x])
        self.dateNum = self.dateNum[sortedIndex]        
        self.dateList = self.dateList[sortedIndex]
        self.Hs = self.Hs[sortedIndex]
        self.T = self.T[sortedIndex]
        
        years = [0] * len(self.dateList)
        for i in range(len(self.dateList)):
            years[i] = self.dateList[i][0]
            
        trainingYear = self.dateList[0][0] + trainingSize
        cond = years <= trainingYear  
        
        subsetBuoy.Hs = self.Hs[cond]
        subsetBuoy.T = self.T[cond]
        subsetBuoy.dateList = self.dateList[cond]
    
    
        return(subsetBuoy)    

    def _saveData(self, fileObj):
        '''Organizes and saves wave height, energy period, and date data.'''
        if(self.Hs is not None):
            gbd = fileObj.create_group('buoy_Data')
            f_Hs = gbd.create_dataset('Hs', data=self.Hs)
            f_Hs.attrs['units'] = 'm'
            f_Hs.attrs['description'] = 'significant wave height'
            f_T = gbd.create_dataset('Te', data=self.T)
            f_T.attrs['units'] = 'm'
            f_T.attrs['description'] = 'energy period'
            f_dateNum = gbd.create_dataset('dateNum', data=self.dateNum)
            f_dateNum.attrs['description'] = 'datenum'
            f_dateList = gbd.create_dataset('dateList', data=self.dateList)
            f_dateList.attrs['description'] = 'date list'
        else:
            RuntimeError('Buoy object contains no data')

    def __fetchCDIP(self,savePath,proxy):
        """
        Fetches the Hs and T values of a CDIP site by downloading the respective .nc file from
        http://cdip.ucsd.edu/

        Parameters
        ----------
        savePath : string
            Relative path to place directory with data files.
        """
        url = "http://thredds.cdip.ucsd.edu/thredds/fileServer/cdip/archive/" + str(self.buoyNum) + "p1/" + \
               str(self.buoyNum) +"p1_historic.nc"        
        print("Downloading data from: " + url)    
        filePath = savePath + "/" + str(self.buoyNum) + "-CDIP.nc"
        urllib.request.urlretrieve (url, filePath)

        self.__processCDIPData(filePath)

    def loadCDIP(self, filePath = None):
        """
        Loads the Hs and T values of the given site from the .nc file downloaded from 
        http://cdip.ucsd.edu/
        Parameters
        ----------
            filePath : string
                File path to the respective .nc file containing the Hs and T values
        """
        if filePath == None:
            filePath = "data/" + self.buoyNum + "-CDIP.nc"
        self.__processCDIPData(filePath)

    def __averageValues(self):
        """
        Averages the Hs and T values of the given buoy to get hour time-steps rather than
        half hour time-steps
        """
        self.Hs = np.mean(self.Hs.reshape(-1,2), axis = 1)
        self.T = np.mean(self.T.reshape(-1,2), axis = 1)


    def __processCDIPData(self,filePath):
        """
        Loads the Hs and T values from the .nc file downloaded from http://cdip.ucsd.edu/
        Parameters
        ----------
            filePath : string
                File path to the respective .nc file containing the Hs and T values
        """
        import netCDF4
        try:
            data = netCDF4.Dataset(filePath)
        except IOError:
            raise IOError("Could not find data for CDIP site: " + self.buoyNum)
            
        self.Hs = np.array(data["waveHs"][:], dtype = np.double)
        self.T = np.array(data["waveTa"][:],  dtype = np.double)
        data.close()

        #Some CDIP buoys record data every half hour rather than every hour
        if len(self.Hs)%2 == 0:
            self.__averageValues()



    def _prepData(self):
        '''Runs _getStats and _getDataNums for full set of data, then removes any
        NaNs. This cleans and prepares the data for use. Returns wave height, 
        energy period, and dates.
        '''
        n = len(self.swdList)
        Hs = []
        T = []
        dateNum = []
        for ii in range(n):
            tmp1, tmp2 = _getStats(self.swdList[ii], self.freqList[ii])
            Hs.extend(tmp1)
            T.extend(tmp2)
            dateNum.extend(_getDateNums(self.dateList[ii]))

        dateList = [date for year in self.dateList for date in year]
        Hs = np.array(Hs, dtype=np.float)
        T = np.array(T, dtype=np.float)
        dateNum = np.array(dateNum, dtype=np.float)
        dateList = np.array(dateList)

        # Removing NaN data, assigning T label depending on input (Te or Tp)
        Nanrem = np.logical_not(np.isnan(T) | np.isnan(Hs))
        # Find NaN data in Hs or T
        dateNum = dateNum[Nanrem]  # Remove any NaN data from DateNum
        Hs = Hs[Nanrem]  # Remove any NaN data from Hs
        T = T[Nanrem]  # Remove any NaN data from T
        self.Hs = Hs
        self.T = T
        self.dateNum = dateNum
        self.dateList = dateList

        return Hs, T, dateNum, dateList

def _getDateNums(dateArr):
    '''datetime objects

    Parameters
    ----------
        dateArr : np.array
            Array of a specific years date vals from NDBC.fetchFromWeb

    Returns
    -------
        dateNum : np.array
            Array of datetime objects.
    '''
    dateNum = []
    for times in dateArr:
        if  times[0] < 1900:
            times[0] = 1900 + times[0]
        dateNum.append(date.toordinal(datetime(times[0], times[1],
                                               times[2], times[3])))

    return dateNum

def _getStats(swdArr, freqArr):
        '''Significant wave height and energy period

        Parameters
        ----------
            swdArr : np.array
                Numpy array of the spectral wave density data for a specific year
            freqArr: np.array
                Numpy array that contains the frequency values for a specific year

        Returns
        -------
            Hm0 : list
                Significant wave height.
            Te : list
                Energy period.
        '''
        #Ignore divide by 0 warnings and resulting NaN warnings
        np.seterr(all='ignore')

        Hm0 = []
        T = []

        for line in swdArr:
            m_1 = np.trapz(line * freqArr ** (-1), freqArr)
            m0 = np.trapz(line, freqArr)
            Hm0.append(4.004 * m0 ** 0.5)
            T.append(m_1 / m0)
        np.seterr(all='warn')

        return Hm0, T
