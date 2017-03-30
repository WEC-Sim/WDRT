import numpy as np
import WDRT.ESSC as ESSC
import copy
import matplotlib.pyplot as plt

# Create buoy object, in this case for Station #46022
buoy46022 = ESSC.Buoy('46022')

# Read data from ndbc.noaa.gov
# buoy46022.fetchFromWeb()

# Load data from .txt file if avilable
buoy46022.loadFromText()

# Load data from .h5 file if available
#buoy46022.loadFromH5()

# Declare required parameters
depth = 391.4  # Depth at measurement point (m)
size_bin = 250.  # Enter chosen bin size

# # Create EA object using above parameters
pca46022 = ESSC.PCA(depth, size_bin, buoy46022)
Gauss46022 = ESSC.GaussianCopula(depth, buoy46022)
Gumbel46022 = ESSC.GumbelCopula(depth, buoy46022)
cc46022 = ESSC.ClaytonCopula(depth, buoy46022)
rosen46022 = ESSC.Rosenblatt(depth, buoy46022)

Time_SS = 1.  # Sea state duration (hrs)
Time_R = 100  # Return periods (yrs) of interest

pca_Hs_Return, pca_T_Return = pca46022.getContours(Time_SS, Time_R)
Gauss_Hs_Return, Gauss_T_Return = Gauss46022.getContours(Time_SS, Time_R)
Gumbel_Hs_Return, Gumbel_T_Return = Gumbel46022.getContours(Time_SS, Time_R)
cc_Hs_Return, cc_T_Return = cc46022.getContours(Time_SS, Time_R)
rosen_Hs_Return, rosen_T_Return = rosen46022.getContours(Time_SS, Time_R)

f = plt.figure()
f.canvas.set_window_title('NDBC%s, %i-year contours' % (buoy46022.buoyNum, Time_R))
plt.plot(buoy46022.T, buoy46022.Hs, 'bo', alpha=0.1, label='Data')
plt.plot(pca_T_Return, pca_Hs_Return, '-', label='PCA')
plt.plot(Gauss_T_Return, Gauss_Hs_Return, '-', label='Gaussian')
plt.plot(Gumbel_T_Return, Gumbel_Hs_Return, '-', label='Gumbel')
plt.plot(cc_T_Return, cc_Hs_Return, '-', label='Clayton')
plt.plot(rosen_T_Return, rosen_Hs_Return, '-', label='Rosenblatt')
plt.xlabel('Energy period, $T_e$ [s]')
plt.ylabel('Sig. wave height, $H_s$ [m]')
plt.grid(True)
plt.legend(loc='upper left', fontsize=10, fancybox=True)
plt.show()



# Sample Generation Example
num_contour_points = 20  # Number of points to be sampled for each
# contour interval.
contour_probs = 10 ** (-1 * np.array([1, 2, 2.5, 3, 3.5, 4, 4.5, 5, 5.5, 6]))
# Probabilities defining sampling contour bounds.
random_seed = 2  # Random seed for sample generation

# Get samples for a full sea state long term analysis
Hs_sampleFSS, T_sampleFSS, Weight_sampleFSS = pca46022.getSamples(num_contour_points,
                                                     contour_probs, random_seed)
# Get samples for a contour approach long term analysis
T_sampleCA = np.arange(12, 26, 2)
Hs_sampleCA = pca46022.getContourPoints(T_sampleCA)

# Modify contour by steepness curve if they intersect
SteepMax = 0.07  # Optional: enter estimate of breaking steepness
T_vals = np.arange(0.1, np.amax(buoy46022.T), 0.1)

SteepH = pca46022.steepness(SteepMax, T_vals)
SteepH_Return = pca46022.steepness(SteepMax, pca46022.T_ReturnContours)

Steep_correction = np.where(SteepH_Return < pca46022.Hs_ReturnContours)
Hs_Return_Steep = copy.deepcopy(pca46022.Hs_ReturnContours)
Hs_Return_Steep[Steep_correction] = SteepH_Return[Steep_correction]

pca46022.bootStrap(boot_size=10)

# Save data in h5 file
pca46022.saveData()

# Show a plot of the data
pca46022.plotData()
