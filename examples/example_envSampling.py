import numpy as np
import WDRT.NDBCdata as NDBCdata
import WDRT.ESSC as ESSC
import os
import matplotlib.pyplot as plt
import h5py
import copy


# Pull spectral data from NDBC website
swdList, freqList, dateVals = NDBCdata.fetchFromWeb(46022, savePath='data')

# # Load data from existing text files
# swdList, freqList, dateVals = NDBCdata.loadFromText(
#     os.path.join('data', 'NDBC46022'))

# Run stats and remove NaNs
Hs, T, DateNum = NDBCdata.prepData(swdList, freqList, dateVals)

# Declare required parameters
depth = 391.4  # Depth at measurement point (m)
size_bin = 250.  # Enter chosen bin size
nb_steps = 1000.  # Enter discretization of the circle in the normal space
# used for inverse FORM calculation
Time_SS = 1.  # Sea state duration (hrs)
Time_r = np.array([100])  # Return periods (yrs) of interest
SteepMax = 0.07  # Optional: enter estimate of breaking steepness

# Contour generation example
Hs_Return, T_Return, _, _, _, _, _ = ESSC.getContours(Hs, T, depth, size_bin,
                                                      nb_steps, Time_SS, Time_r)
# Sample Generation Example
num_contour_points = 20  # Number of points to be sampled for each
# contour interval.
contour_probs = 10 ** (-1 * np.array([1, 2, 2.5, 3, 3.5, 4, 4.5, 5, 5.5, 6]))
# Probabilities defining sampling contour bounds.
random_seed = 2  # Random seed for sample generation

# Get samples for a full sea state long term analysis
Hs_sampleFSS, T_sampleFSS, Weight_sampleFSS = ESSC.getSamples(Hs, T, num_contour_points,
                                                     contour_probs, random_seed,
                                                     depth, size_bin, nb_steps,
                                                     Time_SS, Time_r)

# Get samples for a contour approach long term analysis
T_sampleCA = np.arange(12, 26, 2)
Hs_sampleCA = ESSC.getContourPoints(T_Return, Hs_Return, T_sampleCA)

# Modify contour by steepness curve if they intersect
T_vals = np.arange(0.1, np.amax(T), 0.1)
SteepH = ESSC.steepness(depth, SteepMax, T_vals)
SteepH_Return = ESSC.steepness(depth, SteepMax, T_Return)
Steep_correction = np.where(SteepH_Return < Hs_Return)
Hs_Return_Steep = copy.deepcopy(Hs_Return)
Hs_Return_Steep[Steep_correction] = SteepH_Return[Steep_correction]

# Save data for future use
with h5py.File('data/envSamples_NDBC46022.h5', 'w') as f:

    # NDBC data
    f_Hs = f.create_dataset('Hs', data=Hs)
    f_Hs.attrs['units'] = 'm'
    f_Hs.attrs['description'] = 'significant wave height'
    f_T = f.create_dataset('T', data=T)
    f_T.attrs['units'] = 'm'
    f_T.attrs['description'] = 'energy period'

    # Return contours
    f_T_Return = f.create_dataset('T_Return', data=T_Return)
    f_T_Return.attrs['units'] = 's'
    f_T_Return.attrs['description'] = 'contour, energy period'
    f_Hs_Return = f.create_dataset('Hs_Return', data=Hs_Return)
    f_Hs_Return.attrs['units'] = 'm'
    f_Hs_Return.attrs['description'] = 'contours, significant wave height'
    f_Hs_Return_Steep = f.create_dataset(
        'Hs_Return_Steep', data=Hs_Return_Steep)
    f_Hs_Return_Steep.attrs['units'] = 'm'
    f_Hs_Return_Steep.attrs['description'] = 'contours (steepness limited), significant wave height'

    # Samples for full sea state long term analysis
    f_Hs_sampleFSS = f.create_dataset('Hs_sampleFSS', data=Hs_sampleFSS)
    f_Hs_sampleFSS.attrs['units'] = 'm'
    f_Hs_sampleFSS.attrs['description'] = 'full sea state significant wave height samples'
    f_T_sampleFSS = f.create_dataset('T_sampleFSS', data=T_sampleFSS)
    f_T_sampleFSS.attrs['units'] = 's'
    f_T_sampleFSS.attrs['description'] = 'full sea state energy period samples'
    f_Weight_sampleFSS = f.create_dataset('Weight_sampleFSS', data=Weight_sampleFSS)
    f_Weight_sampleFSS.attrs['description'] = 'full sea state relative weighting samples'

    # Samples for contour approach long term analysis
    f_Hs_sampleCA = f.create_dataset('Hs_sampleCA', data=Hs_sampleCA)
    f_Hs_sampleCA.attrs['units'] = 'm'
    f_Hs_sampleCA.attrs['description'] = 'contour approach significant wave height samples'
    f_T_sampleCA = f.create_dataset('T_sampleCA', data=T_sampleCA)
    f_T_sampleCA.attrs['units'] = 's'
    f_T_sampleCA.attrs['description'] = 'contour approach energy period samples'

# Plot data
plt.figure()
plt.plot(T, Hs, 'bo', alpha=0.1, label='NDBC data')
plt.plot(T_Return, Hs_Return, 'k-', label='100 year contour')
plt.plot(T_Return, Hs_Return_Steep, '-', color='0.65',
         label='100 year contour w/ breaking')
plt.plot(T_sampleFSS, Hs_sampleFSS, 'ro', label='full sea state samples')
plt.plot(T_sampleCA, Hs_sampleCA, 'y^', label='contour approach samples')
plt.legend(loc='lower right', fontsize='small')
plt.grid(True)
plt.xlabel('Energy period, $T_e$ [s]')
plt.ylabel('Sig. wave height, $H_s$ [m]')
plt.show()
