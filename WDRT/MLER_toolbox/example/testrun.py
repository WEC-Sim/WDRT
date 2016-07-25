#!/usr/bin/python

libpath     = '/Users/equon/WDRT/WDRT/MLER_toolbox'
RAOdir      = libpath+'/RAO_data/'
outputDir   = libpath+'/example/TestData/'

import sys
if (libpath not in sys.path): sys.path.append(libpath)
import mler


# Create the object
Test = mler.mler(H=9.0, T=15.1, numFreq=500)
Test.sim.setup()

# Setup the wave information
Test.waves.setup()
print Test.waves
Test.waves.plotSpectrum(show=False)

# Setup the RAO information
Test.readRAO(3,RAOdir+'RAO_heave_RM3float.dat')
Test.readRAO(5,RAOdir+'RAO_pitch_RM3float.dat')

# Check the input RAOs
Test.plotRAO(3,show=False)
Test.plotRAO(5)

#
# heave conditioned response
#

# Now that everything is set up, generate the MLER wave for heave.
Test.MLERcoeffsGen(3,1.0)   # generate the wave profile, 1 meter heave response desired

# At this point we can export the coefficients.  The coefficients will
# match a desired response height that was given as the second argument to
# MLERcoeffsGen.  But, if a specific height of the incoming wave is wanted,
# we can renormalize the wave amplitude now. (AP)
Test.MLERwaveAmpNormalize(Test.waves.H/2 * 1.9) # the desired peak height (peak to MSL)

# Now export the heave coefficients.
Test.MLERexportCoeffs(outputDir+'Test_MLER_heaveOpt_Coeffs.txt');

# Export the wave amplitude time series at x=x0 and the heave response.
Test.MLERexportWaveAmpTime(outputDir+'Test_MLER_heaveOpt_heave_WaveAmpTime.txt',3)

# Export the spectral info for WEC-Sim.
Test.MLERexportWECSim(outputDir+'Test_MLER_heaveOpt_WECSimInput.txt')

# make a movie of the sequence for heave
#Test.MLERanimate(3,export='Movie_Heave')        # export without plotting
Test.MLERanimate(3)                             # plot animation
Test.MLERexportMovie(outputDir+'Movie_Heave')   # export generated animation (after plotting)

# make a movie of the sequence for pitch
#Test.MLERanimate(5,export=outputDir+'Movie_Pitch') # export without plotting

# make a movie of the combined heave + pitch response
#Test.MLERanimate2D(export=outputDir+'Movie_2D_heaveOpt')


#
# pitch conditioned response
#

#Test.MLERcoeffsGen(5,1.0)   # generate the wave profile, 1 meter heave response desired
#Test.MLERwaveAmpNormalize(Test.waves.H/2 * 1.9) # the desired peak height (peak to MSL)
#Test.MLERanimate2D(export=outputDir+'Movie_2D_pitchOpt')

