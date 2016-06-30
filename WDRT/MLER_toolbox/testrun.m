% create the object
Test=MLERclass;

% setup the simulation
Test.sim.simSetup;

% setup the waves
Test.waves.H=9;
Test.waves.T=15.1;
Test.waves.numFreq=500;
Test.waves.waveSetup;
Test.waves.plotSpectrum;

% setup the RAO information
Test.RAOread(3,'RAO_data/RAO_heave_RM3float.dat');
Test.RAOread(5,'RAO_data/RAO_pitch_RM3float.dat');
Test.RAOplot(3);        % check the RAO
Test.RAOplot(5);        % check the RAO

% now that everything is setup, generate the MLER wave for heave.
Test.MLERcoeffsGen(3,1);        % generate the wave profile, 1 meter response desired

% at this point we can export the coefficients.  The coefficients will
% match a desired response height that was given as the second argument to
% MLERcoeffsGen.  But, if a specific height of the incoming wave is wanted,
% we can renormalize the wave amplitude now
Test.MLERwaveAmpNormalize(Test.waves.H /2  * 1.9)           % the peak height (peak to MSL) desired

% now export the heave coefficients
Test.MLERexportCoeffs('TestData/Test_heave_MLER_heaveOpt_Coeffs.txt');


% export the wave amplitude time series at x=x0 for heave
Test.MLERexportWaveAmpTime('TestData/Test_heave_MLER_heaveOpt_heave_WaveAmpTime.txt',3)

% export the spectral info for WEC-Sim
Test.MLERexportWECSim('TestData/Test_heave_MLER_heaveOpt_WECSimInput.txt')



% make a movie of the sequence for heave
MovieFramesHeave=Test.MLERanimate(3);
% export movie of heave (don't add extension)
Test.MLERexportMovie('Movie_Heave',MovieFramesHeave)


% make a movie of the sequence for pitch
MovieFramesPitch=Test.MLERanimate(5);
% export movie of pitch (don't add extension)
Test.MLERexportMovie('Movie_Pitch',MovieFramesPitch)


% do some checks for heave vs. pitch optimized
Test_heave.S = Test.S;

Test.MLERcoeffsGen(5,1);
Test_pitch.S = Test.S;
