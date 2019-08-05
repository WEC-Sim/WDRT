#!/usr/bin/env python
import numpy as np
import scipy.interpolate

import WDRT.mler.wave as wave
import WDRT.mler.simulation as simulation
import WDRT.mler.spectrum as spectrum

class mler(object):
    """ Based on MLERclass.m
    """

    def __init__(self,H=None,T=None,numFreq=10001):
        self.sim   = simulation.simulation()
        self.waves = wave.wave(H,T,numFreq)

        self.desiredRespAmp     = 0.0                       # [-]   Desired response amplitude.  M_d in documentation.
        self.waveHeightDesired  = None                      # [m]   Height of wave desired if renormalizing amplitudes

        # calculated variables
        self._RAO               = None                      # [-]   Complex RAO array  N x 6
        self._RAOdataReadIn     = np.zeros(6,dtype=bool)    # [-]   What RAO dimensions did we read in?
        self._RAOdataFileName   = ['','','','','','']       # [-]   Name of the data file read in for RAO
        self._CoeffA_Rn         = None                      # [ ]   MLER coefficients A_{R,n}
        self._S                 = None                      # [m^2-s] Conditioned wave spectrum vector
        self._A                 = None                      # [m^2-s] 2*(conditioned wave spectrum vector)
        self._phase             = 0.0                       # [rad] Wave phase
        self._Spect             = None                      # [-]   Spectral info
        self._rescaleFact       = None                      # [-]   Rescaling factor for renormalizing the amplitude

        self._animation         = None                      #       Object storing animation information
        self._respExtremes      = None                      # [m]   Array containing min/max of the simulated response

    @property
    def Spect(self): return self._Spect
    @property
    def S(self): return self._S
    @property
    def A(self): return self._A
    @property
    def phase(self): return self._phase

    def __repr__(self):
        s = 'MLER focused wave (desired response amplitude= {:f})'.format(self.desiredRespAmp)
        s+= '\n\tphase : {:f} deg'.format(self._phase*180./np.pi)
        s+= '\n\trescale factor : {:f}'.format(self._rescaleFact)
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
        Sets: self._RAO, self._RAOdataReadIn[DOFread], self._RAOdataFileName[DOFread]
        """
        if self._RAO is None:
            self._RAO = np.zeros( (self.waves.numFreq,6), dtype=complex ) # set the size of the RAO matrix

        if self._RAOdataReadIn[DOFread-1] is True:
            print('WARNING: RAO dof=',DOFread,'already read from',self._RAOdataFileName[DOFread-1])

        # make sure we have setup the waves info first.
        if self.waves._w is None:
            raise UnboundLocalError('Call waves.waveSetup before calling ReadRAO')

        
        # Format of file to read in:
        # Column 1:    period in seconds
        # Column 2:    response amplitude (m/m for DOF 1-3; or radians/m for DOF 4-6)
        # Column 3:    response phase (radians)
        #
        print('Reading RAO ( DOF=',DOFread,') from',RAO_File_Name)
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
        #### Questionable if the amplitude should be 1 or 0.  If set to 1, we count
        #### on the spectrum having nothing out there.  For dimensions
        #### other than in heave, this should be valid.  Heave should be
        #### set to 1. (ADP)
        if DOFread == 3: #heave
            tmp = np.array( [[0,1,0]] )
        else:
            tmp = np.array( [[0,0,0]] )
        tmpRAO = np.concatenate( (tmp,tmpRAO), axis=0 )
        
        # Now interpolate to find the values
        Amp   = scipy.interpolate.pchip_interpolate( tmpRAO[:,0], tmpRAO[:,1], self.waves._w )
        Phase = scipy.interpolate.pchip_interpolate( tmpRAO[:,0], tmpRAO[:,2], self.waves._w )

        # create the complex value to return
        self._RAO[:,DOFread-1] = Amp * np.exp(1j*Phase)
        
        # set flag so we know that this dimension was read in, and save filename
        self._RAOdataReadIn[DOFread-1] = True
        self._RAOdataFileName[DOFread-1] = RAO_File_Name;

    def plotRAO(self,DOFtoPlot,show=True):
        # make sure we have setup the waves info first.
        if self._RAOdataReadIn[DOFtoPlot-1] is False:
            raise UnboundLocalError('Call waves.waveSetup and ReadRAO before RAOplot(DOF)');
        import os
        import matplotlib.pyplot as plt
        plt.figure()
        plt.plot( self.waves._w,      abs(self._RAO[:,DOFtoPlot-1]), 'b-' , label='amplitude' )
        plt.plot( self.waves._w, np.angle(self._RAO[:,DOFtoPlot-1]), 'r--', label='phase' )
        basename = os.path.basename( self._RAOdataFileName[DOFtoPlot-1] )
        plt.title( 'RAO for dimension {:d} from file {:s}'.format(DOFtoPlot,basename) )
        plt.xlabel('Frequency (rad/s)')
        #if DOFtoPlot <=3:
        #    plt.ylabel('Response amplitude (m/m), Response phase (rad)')
        #else:
        #    plt.ylabel('Response amplitude (rad/m), Response phase (rad)')
        plt.ylabel('Response amplitude (*/m), Response phase (rad)')
        plt.legend()
        if show is True: plt.show()
        
    def MLERcoeffsGen(self,DOFtoCalc,response_desired=None,safety_factor=None):
        """ This function calculates MLER (most likely extreme response)
        coefficients given a spectrum and RAO

        DOFtoCalc: 1 - 3 (translational DOFs)
                   4 - 6 (rotational DOFs)
        response_desired: desired response, units should correspond to DOFtoCalc
            for a motion RAO or units of force for a force RAO
        safety_factor: alternative to specifying response_desired;
            non-dimensional scaling factor applied to half the significant wave
            height

        Sets self._S, self._A, self._CoeffA_Rn, self._phase
        Sets self._Spect containing spectral information
        """
        DOFtoCalc -= 1 # convert to zero-based indices (EWQ)
        
        # check that we specified a response
        if (response_desired is None) and (safety_factor is None):
            raise ValueError('Specify response_desired or safety_factor.')
        elif safety_factor:
           #RAO_Tp = np.interp(2.0*np.pi/self.waves.T,self.waves._w,np.abs(self._RAO[:,DOFtoCalc]))
            RAO_Tp = scipy.interpolate.pchip_interpolate(self.waves._w,
                                                         np.abs(self._RAO[:,DOFtoCalc]),
                                                         2.0*np.pi/self.waves.T)
            response_desired = np.abs(RAO_Tp) * safety_factor*self.waves.H/2
            print('Target wave elevation         :',safety_factor*self.waves.H/2)
            print('Interpolated RAO(Tp)          :',RAO_Tp)
            print('Desired response (calculated) :',response_desired)
        self.desiredRespAmp = response_desired

        S_R             = np.zeros(self.waves.numFreq)  # [(response units)^2-s/rad]
        self._S         = np.zeros(self.waves.numFreq)  # [m^2-s/rad]
        self._A         = np.zeros(self.waves.numFreq)  # [m^2-s/rad]
        self._CoeffA_Rn = np.zeros(self.waves.numFreq)  # [1/(response units)]
        self._phase     = np.zeros(self.waves.numFreq)
        
        # calculate the RAO times sqrt of spectrum
        # note that we could define:  a_n=(waves.A*waves.dw).^0.5; (AP)
        #S_tmp(:)=squeeze(abs(obj.RAO(:,DOFtoCalc))).*2 .* obj.waves.A;          % Response spectrum.
        # note: self.A == 2*self.S  (EWQ)
        #   i.e. S_tmp is 4 * RAO * calculatedeWaveSpectrum
       #S_tmp[:] = 2.0*np.abs(self._RAO[:,DOFtoCalc])*self.waves._A     # Response spectrum.

        # Note: waves.A is "S" in Quon2016; 'waves' naming convention matches WEC-Sim conventions (EWQ)
        S_R[:] = np.abs(self._RAO[:,DOFtoCalc])**2 * self.waves._A  # Response spectrum [(response units)^2-s/rad] -- Quon2016 Eqn. 3 

        # calculate spectral moments and other important spectral values.
        self._Spect = spectrum.stats( S_R, self.waves._w, self.waves._dw )
       
        # calculate coefficient A_{R,n} [(response units)^-1] -- Quon2016 Eqn. 8
        self._CoeffA_Rn[:] = np.abs(self._RAO[:,DOFtoCalc]) * np.sqrt(self.waves._A*self.waves._dw) \
                * ( (self._Spect.M2 - self.waves._w*self._Spect.M1) \
                    + self._Spect.wBar*(self.waves._w*self._Spect.M0 - self._Spect.M1) ) \
                / (self._Spect.M0*self._Spect.M2 - self._Spect.M1**2) # Drummen version.  Dietz has negative of this.

        # save the new spectral info to pass out
        self._phase[:] = -np.unwrap( np.angle(self._RAO[:,DOFtoCalc]) ) # Phase delay should be a positive number in this convention (AP)
        
        # for negative values of Amp, shift phase by pi and flip sign
        self._phase[self._CoeffA_Rn < 0]     -= np.pi # for negative amplitudes, add a pi phase shift
        self._CoeffA_Rn[self._CoeffA_Rn < 0] *= -1    # then flip sign on negative Amplitudes
        
        # calculate the conditioned spectrum [m^2-s/rad]
        self._S[:] = self.waves._S * self._CoeffA_Rn[:]**2 * self.desiredRespAmp**2
        self._A[:] = self.waves._A * self._CoeffA_Rn[:]**2 * self.desiredRespAmp**2 # self.A == 2*self.S
        
        # if the response amplitude we ask for is negative, we will add
        # a pi phase shift to the phase information.  This is because
        # the sign of self.desiredRespAmp is lost in the squaring above.
        # Ordinarily this would be put into the final equation, but we
        # are shaping the wave information so that it is buried in the
        # new spectral information, S. (AP)
        if self.desiredRespAmp < 0:
            self._phase += np.pi

    def MLERwaveAmpNormalize(self,peakHeightDesired):
        """ Renormalize the wave amplitude to some desired height of the incoming wave. 
        Uses the peak height (peak to MSL) desired rather than the full range height (peak to trough).
        Sets: self._rescaleFact (for reference)
        Sets: self._S, self._A
        """
        # check that we asked for a positive wave amplitude
        if peakHeightDesired <=0:
            raise ValueError('Wave height desired during renormalization must be positive.')
        self.waveHeightDesired = peakHeightDesired

        print('Renormalizing wave peak height to {:f} m. May take some time depending on spatial and temporal resolution...'.format(peakHeightDesired))
        
        tmpMaxAmp = self._MLERpeakvalue()

        # renormalization of wave amplitudes
        self._rescaleFact = np.abs(peakHeightDesired) / np.abs(tmpMaxAmp)
        self._S = self._S * self._rescaleFact**2 # rescale the wave spectral amplitude coefficients
        self._A = self._A * self._rescaleFact**2 # rescale the wave amplitude coefficients
        print('Rescaled by {:f}'.format(self._rescaleFact))
        
    def MLERexportCoeffs(self,FileNameCoeff):
        """ Export coefficients to use as input for other codes (e.g. Star, ...)
        """
        import datetime
        
        Phase =  self._phase + self.waves._w*self.sim.T0 - self.waves._k*self.sim.X0  # note sign: overall exported phase is still backwards (AP)

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
            f.write('# dW:      {:8.5g}  (rad/s, frequency spacing)\n'.format(self.waves._dw))
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
            f.write('#   (rad/s)        (m^2-s)       (rad)     (rad^2/m)\n')

            # Write coefficients, etc
            for i,wi in enumerate(self.waves._w):
                f.write('{:8.6f}   {:12.8g}   {:12.8f}   {:12.8f}\n'.format(
                        wi, self._S[i], Phase[i], self.waves._k[i] ) )

        print('MLER coefficients written to',FileNameCoeff)

    def MLERexportWaveAmpTime(self,FileNameWaveAmpTime,DOFexport):
        """Export the wave amplitude timeseries at X0 to a file and the
        response for a specified DOF.

        DOFexport: 1 - 3 (translational DOFs)
                   4 - 6 (rotational DOFs)
        """
        import datetime

        # calculate the series
        waveAmpTime = np.zeros( (self.sim.maxIT,2) )
        xi = self.sim.X0
        for i,ti in enumerate(self.sim.T):
            
            # conditioned wave
            waveAmpTime[i,0] = np.sum( 
                    np.sqrt(self._A*self.waves._dw) *
                        np.cos( self.waves._w*(ti-self.sim.T0) + self._phase - self.waves._k*(xi-self.sim.X0) )
                    )
            
            # Response calculation
            waveAmpTime[i,1] = np.sum( 
                    np.sqrt(self._A*self.waves._dw) * np.abs(self._RAO[:,DOFexport-1]) *
                        np.cos( self.waves._w*(ti-self.sim.T0) - self.waves._k*(xi-self.sim.X0) )
                    )
            
        print('Exporting wave amplitude time series for DOF =',DOFexport,'at X0.')
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
            f.write('# dW:      {:8.5g}  (rad/s, frequency spacing)\n'.format(self.waves._dw))
            f.write('# Hs:      {:6.3f}  (m, significant wave height) \n'.format(self.waves.H))
            f.write('# Tp:      {:6.3f}  (s, wave period) \n'.format(self.waves.T))
            f.write('#\n')
            f.write('# Max wave elevation:  {:6.3f}  (m) \n'.format(np.max(waveAmpTime[:,0])))
            f.write('# Min wave elevation:  {:6.3f}  (m) \n'.format(np.min(waveAmpTime[:,0])))
            f.write('# Max WEC Response:    {:6.3f}  {:s} \n'.format(np.max(waveAmpTime[:,1]),
                self.sim.DOFunits[DOFexport-1]))
            f.write('# Min WEC Response:    {:6.3f}  {:s} \n'.format(np.min(waveAmpTime[:,1]),
                self.sim.DOFunits[DOFexport-1]))
            f.write('#\n')
            f.write('# Note: Phase is calculated at X0 and t0.  For starting at a different point along x or in time, phase must be adjusted.\n')
            f.write('#\n')
            f.write('#Form of equation for wave elevation:\n')
            f.write('#   WaveElev = sum( sqrt(2*SpectAmp * dw) * cos( -k*(x-X0) + w*(t-T0) + Phase) ) \n')
            f.write('#\n')
            f.write('#\n')
            #f.write('#   time      WaveHeight      LinearResp ({:s})\n'.format(self.sim.DOFnames[DOFexport-1]))
            #f.write('#   (s)          (m)             {:s}\n'.format(self.sim.DOFunits[DOFexport-1]))
            f.write('#   time      WaveHeight      LinearResp\n'.format(self.sim.DOFnames[DOFexport-1]))
            f.write('#   (s)          (m)             (*)\n'.format(self.sim.DOFunits[DOFexport-1]))

            for i,ti in enumerate(self.sim.T):
                f.write('{:12.8f}   {:12.8f}   {:12.8f}\n'.format(ti,waveAmpTime[i,0],waveAmpTime[i,1]))
        
        print('MLER wave amplitude time series written to',FileNameWaveAmpTime)

    def MLERexportWECSim(self,FileNameWEC):
        """ Export the coefficients to a file that WEC-Sim can read in
        """
        # note that:
        #   WaveElev = sum( sqrt(2*S * dw) * cos( -k*(x-X0) + w*(t-T0) + Phase) )
        Freq = self.waves._w / (2*np.pi)

        self._checkpath(FileNameWEC)
        with open(FileNameWEC,'w') as f:
            f.write(('{:8.6f}      '*self.waves.numFreq).format(*Freq))       # output in hertz
            f.write('\n')
            f.write(('{:8.6f}      '*self.waves.numFreq).format(*self._S))
            f.write('\n')
            f.write(('{:8.6f}      '*self.waves.numFreq).format(*self._phase))
            f.write('\n')
        
        print('MLER coefficients for WEC-Sim written to',FileNameWEC)

    def MLERanimate(self,DOF=3,export=None,fps=25):
        """ Animate the MLER results so that I can see what is happening.
        DOF: 1 - 3 (translational DOFs)
             4 - 6 (rotational DOFs)
        export: specify video filename without prefix, or None to play on screen
        """
        import matplotlib.pyplot as plt
        import matplotlib.patches as mpatches
        import matplotlib.animation as anim
        print('Generating animation of wave profile and response for DOF =',DOF)

        # create the 2D dataset
        waveAmpTime = np.zeros( (self.sim.maxIX, self.sim.maxIT, 2) )
        for ix,x in enumerate(self.sim.X):
            for it,t in enumerate(self.sim.T):
                
                # conditioned wave
                waveAmpTime[ix,it,0] = np.sum( np.sqrt(self._A*self.waves._dw) * 
                        np.cos( self.waves._w*(t-self.sim.T0) + self._phase - self.waves._k*(x-self.sim.X0) )
                        )
                
                # Response calculation
                waveAmpTime[ix,it,1] = np.sum( np.sqrt(self._A*self.waves._dw) * np.abs(self._RAO[:,DOF-1]) *
                        np.cos( self.waves._w*(t-self.sim.T0) - self.waves._k*(x-self.sim.X0) )
                        )

        maxval =      max( np.max(waveAmpTime[:,:,0]), np.max(waveAmpTime[:,:,1]) )
        minval = abs( min( np.min(waveAmpTime[:,:,0]), np.min(waveAmpTime[:,:,1]) ) )

        fig, ax = plt.subplots()
        ax.axis( [ self.sim.startX, self.sim.endX, -1.1*max(maxval,minval), 1.1*max(maxval,minval) ] )
        ax.set_autoscale_on(False)
        
        # labels
        ax.set_xlabel('x (m)')
        ax.set_title('MLER '+self.sim.DOFnames[DOF-1])
        ax.set_ylabel('Amplitude (m), Displacement {:s}'.format(self.sim.DOFunits[DOF-1]))
        #ax.legend(ax,'MLER wave'); %,'MLER Response');
        
        # plot first timestep
        waveElev,  = ax.plot( self.sim.X, waveAmpTime[:,0,0], 'b-' )
        floatResp, = ax.plot( self.sim.X, waveAmpTime[:,0,1], 'm--' )
        ax.plot( [self.sim.X0, self.sim.X0], ax.get_ylim(), 'r-' ) # vertical line at X0
        
        # place a rectangle for the object on the plot
        dimX = 0.01 * np.diff(ax.get_xlim())[0]
        dimY = 0.01 * np.diff(ax.get_ylim())[0]
        tmpZ0val = np.interp( self.sim.X0, self.sim.X, waveAmpTime[:,0,1] )
        if DOF==3:
            xy = [ (self.sim.X0-dimX), (tmpZ0val-dimY) ]
            rect = mpatches.Rectangle( xy, 2*dimX, 2*dimY, edgecolor='r', fill=False ) # represents the float
            ax.add_patch(rect)
        
        # make horizontal line for the extreme values
        self._respExtremes = np.zeros(2)
        floatMin, = ax.plot( [-1.5*dimX,1.5*dimX], [0,0], 'k', marker='.', linestyle='-')
        floatMax, = ax.plot( [-1.5*dimX,1.5*dimX], [0,0], 'k', marker='.', linestyle='-')
        minText = plt.text( 1.5*dimX, 0, '', fontsize=10, horizontalalignment='left', verticalalignment='top')
        maxText = plt.text( 1.5*dimX, 0, '', fontsize=10, horizontalalignment='left', verticalalignment='bottom')
        
        # put timestamp on plot
        tstr = 't = {:5.2f} s'.format(0)
        tstamp = plt.text( .98, .97, tstr, fontsize=14,
                horizontalalignment='right', verticalalignment='top',
                transform=ax.transAxes ) # transform to axis coordinates

        def animate(it):
            # update wave surface
            waveElev.set_ydata( waveAmpTime[:,it,0] )

            # update device response surface
            floatResp.set_ydata( waveAmpTime[:,it,1] )

            # udpate device position
            tmpZ0val = np.interp( self.sim.X0, self.sim.X, waveAmpTime[:,it,1] )
            if DOF==3:
                rect.set_y(tmpZ0val-dimY)

            # update device extremes
            self._respExtremes[0] = min( self._respExtremes[0], tmpZ0val )
            self._respExtremes[1] = max( self._respExtremes[1], tmpZ0val )
            floatMin.set_ydata( [self._respExtremes[0],self._respExtremes[0]] )
            floatMax.set_ydata( [self._respExtremes[1],self._respExtremes[1]] )

            minText.set_text('{:f}'.format(self._respExtremes[0]))
            maxText.set_text('{:f}'.format(self._respExtremes[1]))
            minText.set_y(self._respExtremes[0]-dimY)
            maxText.set_y(self._respExtremes[1]+dimY)

            # update timestamp
            tstr = 't = {:5.2f} s'.format(self.sim.T[it])
            tstamp.set_text(tstr)

        def init_animate():
            self._respExtremes = np.zeros(2)

        # notes:
        # init_func: def init() only required for blitting to give a clean slate
        # interval: screen udpate rate
        # blit: blitting redraws only updated portions of the screen; causes problems on OSX
        self.animation = anim.FuncAnimation( fig, animate, frames=np.arange(self.sim.maxIT),
                init_func=init_animate, interval=1000/int(fps), repeat=False, blit=False )

        if export is None:
            plt.show()
        else:
            self.MLERexportMovie(export)

        print('Simulated response extremes:',self._respExtremes,self.sim.DOFunits[DOF-1])

    def MLERanimate2D(self,export=None,fps=25):
        """ Animate 2D (heave + pitch) MLER results
        export: specify video filename without prefix, or None to play on screen
        """
        import matplotlib.pyplot as plt
        import matplotlib.patches as mpatches
        import matplotlib.animation as anim
        heaveDOF = 2
        pitchDOF = 4

        # create the 2D dataset
        waveAmpTime = np.zeros( (self.sim.maxIX, self.sim.maxIT, 3) )
        for ix,x in enumerate(self.sim.X):
            for it,t in enumerate(self.sim.T):
                
                # conditioned wave
                waveAmpTime[ix,it,0] = np.sum( np.sqrt(self._A*self.waves._dw) * 
                        np.cos( self.waves._w*(t-self.sim.T0) + self._phase - self.waves._k*(x-self.sim.X0) )
                        )
                
                # Heave response calculation
                waveAmpTime[ix,it,1] = np.sum( np.sqrt(self._A*self.waves._dw) * np.abs(self._RAO[:,heaveDOF]) *
                        np.cos( self.waves._w*(t-self.sim.T0) - self.waves._k*(x-self.sim.X0) )
                        )

                # Pitch response calculation
                waveAmpTime[ix,it,2] = np.sum( np.sqrt(self._A*self.waves._dw) * np.abs(self._RAO[:,pitchDOF]) *
                        np.cos( self.waves._w*(t-self.sim.T0) - self.waves._k*(x-self.sim.X0) )
                        )

        maxval =      max( np.max(waveAmpTime[:,:,0]), np.max(waveAmpTime[:,:,1]) )
        minval = abs( min( np.min(waveAmpTime[:,:,0]), np.min(waveAmpTime[:,:,1]) ) )

        fig, ax = plt.subplots()
        ax.axis( [ self.sim.startX, self.sim.endX, -1.1*max(maxval,minval), 1.1*max(maxval,minval) ] )
        ax.set_autoscale_on(False)
        
        # labels
        ax.set_xlabel('x (m)')
        ax.set_title('MLER 2-DOF Simulation')
        ax.set_ylabel('Amplitude (m), Displacements {:s},{:s}'.format(
            self.sim.DOFunits[heaveDOF],self.sim.DOFunits[pitchDOF]))
        #ax.legend(ax,'MLER wave'); %,'MLER Response');
        
        # plot first timestep
        waveElev,  = ax.plot( self.sim.X, waveAmpTime[:,0,0], 'b-' )
        floatResp, = ax.plot( self.sim.X, waveAmpTime[:,0,1], 'm--' )
        ax.plot( [self.sim.X0, self.sim.X0], ax.get_ylim(), 'r-' ) # vertical line at X0
        
        # place a rectangle for the object on the plot
        dimX = 0.01 * np.diff(ax.get_xlim())[0]
        dimY = 0.01 * np.diff(ax.get_ylim())[0]
        tmpZ0val = np.interp( self.sim.X0, self.sim.X, waveAmpTime[:,0,1] )
        xy = [ (self.sim.X0-dimX), (tmpZ0val-dimY) ]
        rect = mpatches.Rectangle( xy, 2*dimX, 2*dimY, edgecolor='r', fill=False ) # represents the float
        ax.add_patch(rect)

        # put orientation on plot
        oriline, = plt.plot( [-2*dimX,2*dimX], [0,0], 'g-' )
        oritext = plt.text( self.sim.X0+2*dimX, 2*dimY, r'0.0$^\circ$', fontsize=10, color='g',
                horizontalalignment='left', verticalalignment='bottom')
        
        # put timestamp on plot
        tstr = 't = {:5.2f} s'.format(0)
        tstamp = plt.text( .98, .97, tstr, fontsize=14,
                horizontalalignment='right', verticalalignment='top',
                transform=ax.transAxes ) # transform to axis coordinates

        def animate(it):
            # update wave surface
            waveElev.set_ydata( waveAmpTime[:,it,0] )

            # update device response surface
            floatResp.set_ydata( waveAmpTime[:,it,1] )

            # udpate device position and orientation
            tmpZ0val = np.interp( self.sim.X0, self.sim.X, waveAmpTime[:,it,1] )
            tmpAngval = np.interp( self.sim.X0, self.sim.X, waveAmpTime[:,it,2] )
            ang = tmpAngval * 180./np.pi

            rect.set_y(tmpZ0val-dimY)

            xori =  2*dimX*np.cos(ang*np.pi/180.)
            yori = -2*dimX*np.sin(ang*np.pi/180.) * dimY/dimX # need to correct for aspect ratio
            oriline.set_xdata([self.sim.X0-xori,self.sim.X0+xori])
            oriline.set_ydata([   tmpZ0val-yori,   tmpZ0val+yori])

            oritext.set_y(tmpZ0val+2*dimY)
            oritext.set_text(r'{:f}$^\circ$'.format(ang))

            # update timestamp
            tstr = 't = {:5.2f} s'.format(self.sim.T[it])
            tstamp.set_text(tstr)

        # notes:
        # init_func: def init() only required for blitting to give a clean slate
        # interval: screen udpate rate
        # blit: blitting redraws only updated portions of the screen; causes problems on OSX
        self.animation = anim.FuncAnimation( fig, animate, frames=np.arange(self.sim.maxIT),
                init_func=None, interval=1000/int(fps), repeat=False, blit=False )

        if export is None:
            plt.show()
        else:
            self.MLERexportMovie(export)

    def MLERexportMovie(self,exportName):
        if self.animation: # already run MLERanimate
            fname = exportName+'.mp4'
            print('Exporting animation to',fname)
            self.animation.save(fname, writer='ffmpeg')
        else:
            raise UnboundLocalError('Need to run MLERanimate first')
        
    #
    # protected methods
    #
    def _MLERpeakvalue(self):
        """ the maximum may not occur at X0 or T0... 
        So, we have to generate the entire time and space array, then find the maximum and minimum.
        """
        waveAmpTime = np.zeros( (self.sim.maxIX, self.sim.maxIT) )
        for ix,x in enumerate(self.sim.X):
            for it,t in enumerate(self.sim.T):
                
                # conditioned wave
                waveAmpTime[ix,it] = np.sum(
                        np.sqrt(self._A*self.waves._dw) * # A == 2*S
                            np.cos( self.waves._w*(t-self.sim.T0) - self.waves._k*(x-self.sim.X0) + self._phase ) 
                        )

        return np.max(np.abs(waveAmpTime))

    def _checkpath(self,path):
        """ Create directory if it does not exist
        """
        import os
        dirname = os.path.dirname(path)
        if not os.path.isdir(dirname):
            os.makedirs(dirname)

