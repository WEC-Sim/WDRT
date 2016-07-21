#!/usr/bin/python
# TODO: put this in egg
import sys
import numpy as np
import scipy.interpolate

import wave
import simulation
import spectrum

class mler(object):
    """ Based on MLERclass.m
    """

    def __init__(self,H=None,T=None,numFreq=10001):
        self.sim   = simulation.simulation()
        self.waves = wave.wave(H,T,numFreq)

        self.desiredRespAmp = 0.0                           # [-]   Desired response amplitude.  M_d in documentation.

        # TODO: private set, public get
        self.RAO                = None                      # [-]   Complex RAO array  N x 6
        self.RAOdataReadIn      = np.zeros(6,dtype=bool)    # [-]   What RAO dimensions did we read in?
        self.RAOdataFileName    = ['','','','','','']       # [-]   Name of the data file read in for RAO
        self.CoeffA_Rn          = None                      # [ ]   MLER coefficients A_{R,n}
        self.S                  = None                      # [m^2] New Wave spectrum vector
        self.A                  = None                      # [m^2] 2*(new wave spectrum vector)
        self.phase              = 0.0                       # [rad] Wave phase
        self.Spect              = None                      # [-]   Spectral info
       #self.waveHeightDesired  = None                      # [m]   Height of wave desired if renormalizing amplitudes
       #self.rescaleFact        = None                      # [-]   Rescaling factor for renormalizing the amplitude

    def __repr__(self):
        s = 'MLER focused wave (desired response amplitude= {:f})'.format(self.desiredRespAmp)
        s+= '\n\tphase : {:f} deg'.format(self.phase*180./np.pi)
        #s+= '\n\trescale factor : {:f}'.format(self.rescaleFact)
        return s

    #
    # public methods
    #
    def setup(self):
        """ Shortcut for initializing simulation and waves
        """
        self.sim.setup()
        self.waves.setup()

    def readRAO(self,DOFread,RAO_File_Name):
        """ Read in the RAO from the specified file and assign it to a dimension
        DOFread : 1 - 3 (translational DOFs)
                  4 - 6 (rotational DOFs)
        """
        if self.RAO is None:
            self.RAO = np.zeros( (self.waves.numFreq,6), dtype=complex ) # set the size of the RAO matrix

        if self.RAOdataReadIn[DOFread-1] is True:
            print 'WARNING: RAO dof=',DOFread,'already read from',self.RAOdataFileName[DOFread]

        # make sure we have setup the waves info first.
        if self.waves.w is None:
            sys.exit('Call waves.waveSetup before calling ReadRAO')

        
        # Format of file to read in:
        # Column 1:    period in seconds
        # Column 2:    response amplitude (m/m for DOF 1-3; or radians/m for DOF 4-6)
        # Column 3:    response phase (radians)
        #
        print 'Reading RAO ( DOF=',DOFread,') from',RAO_File_Name
        # - get total number of lines
        with open(RAO_File_Name,'r') as f:
            for i,_ in enumerate(f.readlines()): pass
        nData = i+1
        # - get number of header lines
        dataStart = 0
        with open(RAO_File_Name,'r') as f:
            for line in f:
                try:
                    float(line.split()[0])
                    break
                except:
                    dataStart += 1
        # - preallocate and read data
        nData -= dataStart
        tmpRAO = np.zeros((nData,3))
        with open(RAO_File_Name,'r') as f:
            lines = f.readlines()
            for i,line in enumerate(lines[dataStart:]):
                tmpRAO[i,:] = [ float(val) for val in line.split() ]
        
        # convert from period in seconds to frequency in rad/s
        T = tmpRAO[:,0]
        tmpRAO[:,0] = 2*np.pi / T
        tmpRAO = tmpRAO[ np.argsort(tmpRAO[:,0]) ]  # sort by frequency
        
        # Add at w=0, amp=0, phase=0
        # TODO check:
        #### Questionable if the amplitude should be 1 or 0.  If set to 1, we count
        #### on the spectrum having nothing out there.  For dimensions
        #### other than in heave, this should be valid.  Heave should be
        #### set to 1. (AP)
        if DOFread == 3: #heave
            tmp = np.array( [[0,1,0]] )
        else:
            tmp = np.array( [[0,0,0]] )
        # TODO rewrite w/o concatenate
        tmpRAO = np.concatenate( (tmp,tmpRAO), axis=0 )
        
        # Now interpolate to find the values
        Amp   = scipy.interpolate.pchip_interpolate( tmpRAO[:,0], tmpRAO[:,1], self.waves.w )
        Phase = scipy.interpolate.pchip_interpolate( tmpRAO[:,0], tmpRAO[:,2], self.waves.w )

        # create the complex value to return
        self.RAO[:,DOFread-1] = Amp * np.exp(1j*Phase)
        
        # set flag so we know that this dimension was read in, and save filename
        self.RAOdataReadIn[DOFread-1] = True
        self.RAOdataFileName[DOFread-1] = RAO_File_Name;

    def plotRAO(self,DOFtoPlot,show=False):
        # make sure we have setup the waves info first.
        if self.RAOdataReadIn[DOFtoPlot-1] is False:
            sys.exit('Call waves.waveSetup and ReadRAO before RAOplot(DOF)');
        import matplotlib.pyplot as plt
        plt.figure()
        plt.plot( self.waves.w,      abs(self.RAO[:,DOFtoPlot-1]), 'b-' )
        plt.plot( self.waves.w, np.angle(self.RAO[:,DOFtoPlot-1]), 'r--' )
        plt.title( 'RAO for dimension {:d} from file {:s}'.format(DOFtoPlot,self.RAOdataFileName[DOFtoPlot-1]) )
        plt.xlabel('Frequency (rad/s)')
        if DOFtoPlot <=3:
            plt.ylabel('Response amplitude (m/m) / Response phase (rad)')
        else:
            plt.ylabel('Response amplitude (rad/m) / Response phase (rad)')
        if show is True: plt.show()
        
    def MLERcoeffsGen(self,DOFtoCalc,respDesired):
        """ This function calculates MLER (most likely extreme response) coefficients given a spectrum and RAO
        DOFtoCalc: 1 - 3 (translational DOFs)
                   4 - 6 (rotational DOFs)
        Sets self.S, self.A, self.CoeffA_Rn, self.phase
        """
        # check that we asked for something non-zero
        if respDesired == 0:
            sys.exit('Desired response amplitude (respDesired) should be non-zero.')
        self.desiredRespAmp = respDesired

        DOFtoCalc -= 1 # convert to zero-based indices (EWQ)
        
        S_tmp          = np.zeros(self.waves.numFreq);
        self.S         = np.zeros(self.waves.numFreq);
        self.A         = np.zeros(self.waves.numFreq);
        self.CoeffA_Rn = np.zeros(self.waves.numFreq);
        self.phase     = np.zeros(self.waves.numFreq);
        
        # calculate the RAO times sqrt of spectrum
        # TODO: add equation references
        # note that we could define:  a_n=(waves.A*waves.dw).^0.5; (AP)
        #S_tmp(:)=squeeze(abs(obj.RAO(:,DOFtoCalc))).*2 .* obj.waves.A;          % Response spectrum.
        # note: self.A == 2*self.S  (EWQ)
        #   i.e. S_tmp is 4 * RAO * calculatedWaveSpectrum
        S_tmp[:] = 2.0*np.abs(self.RAO[:,DOFtoCalc])*self.waves.A     # Response spectrum.

        # calculate spectral moments and other important spectral values.
        self.Spect = spectrum.stats( S_tmp, self.waves.w, self.waves.dw )
       
        # calculate coefficient A_{R,n}
        self.CoeffA_Rn[:] = np.abs(self.RAO[:,DOFtoCalc]) * np.sqrt(self.waves.A*self.waves.dw) \
                * ( (self.Spect.M2 - self.waves.w*self.Spect.M1) \
                    + self.Spect.wBar*(self.waves.w*self.Spect.M0 - self.Spect.M1) ) \
                / (self.Spect.M0*self.Spect.M2 - self.Spect.M1**2) # Drummen version.  Dietz has negative of this.
        
        # save the new spectral info to pass out
        self.phase[:] = -np.unwrap( np.angle(self.RAO[:,DOFtoCalc]) ) # Phase delay should be a positive number in this convention (AP)
        
        # for negative values of Amp, shift phase by pi and flip sign
        # TODO: verify this is legit
        self.phase[self.CoeffA_Rn < 0]     -= np.pi # for negative amplitudes, add a pi phase shift
        self.CoeffA_Rn[self.CoeffA_Rn < 0] *= -1    # then flip sign on negative Amplitudes
        
        self.S[:] = self.waves.S * self.CoeffA_Rn[:]**2 * self.desiredRespAmp**2;
        self.A[:] = self.waves.A * self.CoeffA_Rn[:]**2 * self.desiredRespAmp**2;
        
        # if the response amplitude we ask for is negative, we will add
        # a pi phase shift to the phase information.  This is because
        # the sign of self.desiredRespAmp is lost in the squaring above.
        # Ordinarily this would be put into the final equation, but we
        # are shaping the wave information so that it is buried in the
        # new spectral information, S. (AP)
        if self.desiredRespAmp < 0:
            self.phase += np.pi

    def MLERwaveAmpNormalize(self,peakHeightDesired):
        """ Renormalize the wave amplitude to some desired height of the incoming wave. 
        Uses the peak height (peak to MSL) desired rather than the full range height (peak to trough).
        """
        # check that we asked for a positive wave amplitude
        if peakHeightDesired <=0:
            sys.exit('Wave height desired during renormalization must be positive.')
        print 'Renormalizing wave peak height to {:f} m. May take some time depending on spatial and temporal resolution...'.format(peakHeightDesired)
        
        # TODO: HIGHER-ORDER CALCULATION
        tmpMaxAmp = self._MLERpeakvalue()

        # renormalization of wave amplitudes
        self.rescaleFact = np.abs(peakHeightDesired) / np.abs(tmpMaxAmp)
        self.S = self.S * self.rescaleFact**2 # rescale the wave spectral amplitude coefficients
        self.A = self.A * self.rescaleFact**2 # rescale the wave amplitude coefficients
        print 'Rescaled by {:f}'.format(self.rescaleFact)
        
    def MLERexportCoeffs(self,FileNameCoeff):
        import datetime
        
        Phase =  self.phase + self.waves.w*self.sim.T0 - self.waves.k*self.sim.X0  # note sign: overall exported phase is still backwards (AP)

        # Now export the coefficients to a file
        self._checkpath(FileNameCoeff)
        with open(FileNameCoeff,'w') as f:
        
            # Header info
            f.write('# MLER wave profile generated {:}\n'.format(datetime.datetime.now()))
            f.write('#\n')
            f.write('#\n')
            f.write('### Setup \n')
            f.write('# X0:         {:6.3f}   (m, peak position) \n'.format(self.sim.X0))
            f.write('# g:          {:6.3f}   (m/s^2, gravity) \n'.format(self.waves.g))
            f.write('# depth:      {:6.3f}   (m, water depth) \n'.format(self.waves.waterDepth))
            f.write('# T0:         {:6.3f}   (s, response peak time) \n'.format(self.sim.T0))
            f.write('#\n')
            f.write('### Wave info:\n')
            f.write('# NumFreq: {:6f}    (-, number of frequencies)\n'.format(self.waves.numFreq))
            f.write('# dW:      {:8.5g}  (rad/s, frequency spacing)\n'.format(self.waves.dw))
            f.write('# Hs:      {:6.3f}  (m, significant wave height) \n'.format(self.waves.H))
            f.write('# Tp:      {:6.3f}  (s, wave period) \n'.format(self.waves.T))
            f.write('#\n')
            f.write('# Note: Phase is calculated at X0 and t0.  For starting at a different point along x or in time, phase must be adjusted.\n')
            f.write('#\n')
            f.write('#Form of equation for wave elevation:\n')
            f.write('#   WaveElev = sum( sqrt(2*SpectAmp * dw) * cos( -k*(x-X0) + w*(t-T0) + Phase) ) \n')
            f.write('#\n')
            f.write('#\n')
            f.write('#   frequency      SpectAmp      Phase     wavenumber\n')
            f.write('#   (rad/s)          (m^2)       (rad)      (rad^2/m)\n')

            # Write coefficients, etc
            for i,wi in enumerate(self.waves.w):
                f.write('{:8.6f}   {:12.8g}   {:12.8f}   {:12.8f}\n'.format(
                        wi, self.S[i], Phase[i], self.waves.k[i] ) )

        print 'MLER coefficients written to',FileNameCoeff

    def MLERexportWaveAmpTime(self,FileNameWaveAmpTime,DOFexport):
        """ Export the wave amplitude timeseries at X0 to a file
        DOFexport: 1 - 3 (translational DOFs)
                   4 - 6 (rotational DOFs)
        """
        import datetime

        # calculate the series
        waveAmpTime = np.zeros((self.sim.maxIT,2))
        t = np.arange(self.sim.maxIT)*self.sim.dT + self.sim.startTime
        xi = self.sim.X0
        for i,ti in enumerate(t):
            
            # conditioned wave
            # TODO: check for factor of two in sqrt
            waveAmpTime[i,0] = np.sum( 
                    np.sqrt(self.A*self.waves.dw) *
                        np.cos( self.waves.w*(ti-self.sim.T0) + self.phase - self.waves.k*(xi-self.sim.X0) )
                    )
            
            # Response calculation
            # TODO: check for factor of two in sqrt
            # TODO: check for phase term?
            waveAmpTime[i,1] = np.sum( 
                    np.sqrt(self.A*self.waves.dw) * np.abs(self.RAO[:,DOFexport-1]) *
                        np.cos( self.waves.w*(ti-self.sim.T0) - self.waves.k*(xi-self.sim.X0) )
                    )
            
        print 'Exporting wave amplitude time series for DOF =',DOFexport,'at X0.'
        self._checkpath(FileNameWaveAmpTime)
        with open(FileNameWaveAmpTime,'w') as f:
            
            # Header info
            f.write('# MLER wave profile generated on {:}\n'.format(datetime.datetime.now()))
            f.write('#\n')
            f.write('#\n')
            f.write('### Setup \n')
            f.write('# X0:         {:6.3f}   (m, peak position) \n'.format(self.sim.X0))
            f.write('# g:          {:6.3f}   (m/s^2, gravity) \n'.format(self.waves.g))
            f.write('# depth:      {:6.3f}   (m, water depth) \n'.format(self.waves.waterDepth))
            f.write('# startTime:  {:6.3f}   (s) \n'.format(self.sim.startTime))
            f.write('# endTime:    {:6.3f}   (s) \n'.format(self.sim.endTime))
            f.write('# dt:         {:6.3f}   (s) \n'.format(self.sim.dT))
            f.write('# T0:         {:6.3f}   (s, response peak time) \n'.format(self.sim.T0))
            f.write('#\n')
            f.write('### Wave info:\n')
            f.write('# NumFreq: {:6f}    (-, number of frequencies)\n'.format(self.waves.numFreq))
            f.write('# dW:      {:8.5g}  (rad/s, frequency spacing)\n'.format(self.waves.dw))
            f.write('# Hs:      {:6.3f}  (m, significant wave height) \n'.format(self.waves.H))
            f.write('# Tp:      {:6.3f}  (s, wave period) \n'.format(self.waves.T))
            f.write('#\n')
            f.write('# Max wave elevation:  {:6.3f}  (m) \n'.format(np.max(waveAmpTime[:,0])))
            f.write('# Min wave elevation:  {:6.3f}  (m) \n'.format(np.min(waveAmpTime[:,0])))
            f.write('# Max WEC Response:    {:6.3f}  (m) \n'.format(np.max(waveAmpTime[:,1])))
            f.write('# Min WEC Response:    {:6.3f}  (m) \n'.format(np.min(waveAmpTime[:,1])))
            f.write('#\n')
            f.write('# Note: Phase is calculated at X0 and t0.  For starting at a different point along x or in time, phase must be adjusted.\n')
            f.write('#\n')
            f.write('#Form of equation for wave elevation:\n')
            f.write('#   WaveElev = sum( sqrt(2*SpectAmp * dw) * cos( -k*(x-X0) + w*(t-T0) + Phase) ) \n')
            f.write('#\n')
            f.write('#\n')
            f.write('#   time      WaveHeight      LinearResp\n')
            if DOFexport==5:
                f.write('#   (s)          (m)               (m)\n')
            else:
                f.write('#   (s)          (m)              (rad)\n')

            for i,ti in enumerate(t):
                f.write('{:12.8f}   {:12.8f}   {:12.8f}\n'.format(ti,waveAmpTime[i,0],waveAmpTime[i,1]))
        
        print 'MLER wave amplitude time series written to',FileNameWaveAmpTime

    def MLERexportWECSim(self,FileNameWEC):
        """ Export the coefficients to a file that WEC-Sim can read in
        """
        # note that:
        #   WaveElev = sum( sqrt(2*S * dw) * cos( -k*(x-X0) + w*(t-T0) + Phase) )
        Freq = self.waves.w / (2*np.pi)

        self._checkpath(FileNameWEC)
        with open(FileNameWEC,'w') as f:
            f.write(('{:8.6f}      '*self.waves.numFreq).format(*Freq))       # output in hertz
            f.write('\n')
            f.write(('{:8.6f}      '*self.waves.numFreq).format(*self.S))
            f.write('\n')
            f.write(('{:8.6f}      '*self.waves.numFreq).format(*self.phase))
            f.write('\n')
        
        print 'MLER coefficients for WEC-Sim written to',FileNameWEC

    #
    # protected methods
    #
    def _MLERpeakvalue(self):
        """ the maximum may not occur at X0 or T0... 
        So, we have to generate the entire time and space array, then find the maximum and minimum.
        """
        waveAmpTime = np.zeros((self.sim.maxIX,self.sim.maxIT))
        xarray = np.arange(self.sim.maxIX)*self.sim.dX + self.sim.startX
        tarray = np.arange(self.sim.maxIT)*self.sim.dT + self.sim.startTime
        for ix,x in enumerate(xarray):
            for it,t in enumerate(tarray):
                
                # conditioned wave
                # TODO: check factor of two in sqrt
                waveAmpTime[ix,it] = np.sum(
                        np.sqrt(self.A*self.waves.dw) * 
                            np.cos( self.waves.w*(t-self.sim.T0) - self.waves.k*(x-self.sim.X0) + self.phase )
                        )

        return np.max(np.abs(waveAmpTime))

    def _checkpath(self,path):
        """ Create directory if it does not exist
        """
        import os
        dirname = os.path.dirname(path)
        if not os.path.isdir(dirname):
            os.makedirs(dirname)

