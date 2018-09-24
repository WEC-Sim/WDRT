# -*- coding: utf-8 -*-
"""
Created on Fri Sep 21 13:38:12 2018

@author: nevmart
"""

# Create Test Files for ESSCTest.py 
import numpy as np
import WDRT.ESSC as ESSC

# Proxy information
import urllib2
proxy_support = urllib2.ProxyHandler({"https":"https://wwwproxy.sandia.gov:80"})
opener = urllib2.build_opener(proxy_support)
urllib2.install_opener(opener)



# Load in data from NDBC 
buoy46022 = ESSC.Buoy('46022','NDBC')
#buoy46022.fetchFromWeb(proxy = {"https":"https://wwwproxy.sandia.gov:80"})
     
# Save as h5 
#buoy46022.saveAsH5('TestFiles/')

# Load from h5
buoy46022.loadFromH5('NDBC46022.h5')

# Calculate contours 
# Declare required parameters
Time_SS = 1.  # Sea state duration (hrs)
Time_R = 100  # Return periods (yrs) of interest

# Create EA objects for contour methods
pca46022 = ESSC.PCA(buoy46022)
Gauss46022 = ESSC.GaussianCopula(buoy46022)
Gumbel46022 = ESSC.GumbelCopula(buoy46022)
Clayton46022 = ESSC.ClaytonCopula(buoy46022)
rosen46022 = ESSC.Rosenblatt(buoy46022)
NonParaGauss46022 = ESSC.NonParaGaussianCopula(buoy46022)
NonParaClay46022 = ESSC.NonParaClaytonCopula(buoy46022)
NonParaGum46022 = ESSC.NonParaGumbelCopula(buoy46022)

# Calculate contours for all contour methods
pca_Hs_Return, pca_T_Return = pca46022.getContours(Time_SS, Time_R)
Gauss_Hs_Return, Gauss_T_Return = Gauss46022.getContours(Time_SS, Time_R)
Gumbel_Hs_Return, Gumbel_T_Return = Gumbel46022.getContours(Time_SS, Time_R)
Clayton_Hs_Return, Clayton_T_Return = Clayton46022.getContours(Time_SS, Time_R)
rosen_Hs_Return, rosen_T_Return = rosen46022.getContours(Time_SS, Time_R)
NonParaGau_Hs_Return, NonParaGau_T_Return = NonParaGauss46022.getContours(Time_SS, Time_R)
NonParaClay_Hs_Return, NonParaClay_T_Return = NonParaClay46022.getContours(Time_SS, Time_R)
NonParaGum_Hs_Return, NonParaGum_T_Return = NonParaGum46022.getContours(Time_SS, Time_R)

# Save results 
np.save("TestFiles/pca46022-Hs.npy", pca_Hs_Return)
np.save("TestFiles/pca46022-T.npy", pca_T_Return)
np.save("TestFiles/gauss46022-Hs.npy", Gauss_Hs_Return)
np.save("TestFiles/gauss46022-T.npy", Gauss_T_Return)
np.save("TestFiles/gumbel46022-Hs.npy", Gumbel_Hs_Return)
np.save("TestFiles/gumbel46022-T.npy", Gumbel_T_Return)
np.save("TestFiles/cc46022-Hs.npy", Clayton_Hs_Return)
np.save("TestFiles/cc46022-T.npy", Clayton_T_Return)	
np.save("TestFiles/Gauss-NonPara46022-Hs.npy", NonParaGau_Hs_Return)
np.save("TestFiles/Gauss-NonPara46022-T.npy", NonParaGau_T_Return)
np.save("TestFiles/Gum-NonPara46022-Hs.npy", NonParaGum_Hs_Return)
np.save("TestFiles/Gum-NonPara46022-T.npy", NonParaGum_T_Return)
np.save("TestFiles/CC-NonPara46022-Hs.npy", NonParaClay_Hs_Return)
np.save("TestFiles/CC-NonPara46022-T.npy", NonParaClay_T_Return)
np.save("TestFiles/rosen46022-Hs.npy", rosen_Hs_Return)
np.save("TestFiles/rosen46022-T.npy", rosen_T_Return)
