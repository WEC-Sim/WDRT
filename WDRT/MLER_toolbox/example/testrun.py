#!/usr/bin/python
import MLER

# create the object
Test = MLER.focusedWave(H=9.0, T=15.1, numFreq=500)
Test.sim.setup()
print Test
print Test.waves
print Test.sim

# setup the wave information
Test.waves.setup()
Test.waves.plotSpectrum(show=True)

# setup the RAO information
#TODO: call Test.RAOread(3,'RAO_data/RAO_heave_RM3float.dat');
#TODO: call Test.RAOread(5,'RAO_data/RAO_pitch_RM3float.dat');
#TODO: call Test.RAOplot(3);        # check the RAO
#TODO: call Test.RAOplot(5);        # check the RAO

# now that everything is setup, generate the MLER wave for heave.
#TODO: call Test.MLERcoeffsGen(3,1);        # generate the wave profile, 1 meter response desired

# at this point we can export the coefficients.  The coefficients will
# match a desired response height that was given as the second argument to
# MLERcoeffsGen.  But, if a specific height of the incoming wave is wanted,
# we can renormalize the wave amplitude now
#TODO: call Test.MLERwaveAmpNormalize(Test.waves.H /2  * 1.9)           # the peak height (peak to MSL) desired

# now export the heave coefficients
#TODO: call Test.MLERexportCoeffs('TestData/Test_heave_MLER_heaveOpt_Coeffs.txt');


# export the wave amplitude time series at x=x0 for heave
#TODO: call Test.MLERexportWaveAmpTime('TestData/Test_heave_MLER_heaveOpt_heave_WaveAmpTime.txt',3)

# export the spectral info for WEC-Sim
#TODO: call Test.MLERexportWECSim('TestData/Test_heave_MLER_heaveOpt_WECSimInput.txt')


# TODO: animate, etc
# # make a movie of the sequence for heave
# MovieFramesHeave=Test.MLERanimate(3);
# # export movie of heave (don't add extension)
# Test.MLERexportMovie('Movie_Heave',MovieFramesHeave)
# 
# 
# # make a movie of the sequence for pitch
# MovieFramesPitch=Test.MLERanimate(5);
# # export movie of pitch (don't add extension)
# Test.MLERexportMovie('Movie_Pitch',MovieFramesPitch)
# 
# 
# # do some checks for heave vs. pitch optimized
# Test_heave.S = Test.S;
# 
# Test.MLERcoeffsGen(5,1);
# Test_pitch.S = Test.S;

