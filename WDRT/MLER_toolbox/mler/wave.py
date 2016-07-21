#!/usr/bin/python
import sys
import numpy as np

class wave(object):
    """ Based on waveClassMLER.m
    A stripped down version of the waveClass.m file from WEC-Sim.
    """

    def __init__(self, H=None, T=None, numFreq=10001):
        if numFreq <= 10:
            sys.exit('numFreq must be larger than 10')

        # sea state definition
        self.H               = H                    # [m]       Hs: sea state wave height
        self.T               = T                    # [s]       Tp: sea state time period of waves
        self.numFreq         = numFreq              # [-]       Number of frequencies

        # focused wave parameters
        self.startW          = 0.                   # [rad/s]   Starting frequency
        self.endW            = 2*np.pi              # [rad/s]   Ending frequency
        self.waveDir         = 0.                   # [rad]     Wave direction

        # constants
        self.waterDepth      = 70.                  # [m]       Water depth
        self.g               = 9.801                # [m/s^2]   Gravitational constant

        # TODO: private set, public get
        self.w               = None                 # [rad/s]   Wave frequency vector
        self.dw              = None                 # [rad/s]   Frequency step size
        self.type            = 'Bretschneider'      # [-]       Spectrum type
        self.S               = None                 # [m^2]     Wave spectrum vector
        self.A               = None                 # [m^2]     2*(wave spectrum vector)
        self.k               = None                 # [rad^2/m] Wavenumber array
        self.deepWaterWave   = False                # [-]       Deep water or not, depending on input water depth

    def __repr__(self):
        s = 'waveClass (Hs= {:f} m, Tp= {:f} s)'.format(self.H,self.T)
        s+= '\n\tnumber of frequencies  : {:d}'.format(self.numFreq)
        s+= '\n\twater depth            : {:f} m  (deep water={:s})'.format(self.waterDepth,str(self.deepWaterWave))
        s+= '\n\tgravitational constant : {:f} m/s^2'.format(self.g)
        s+= '\n\tfrequency range        : [ {:f} {:f} ] rad/s'.format(self.startW,self.endW)
        s+= '\n\twave direction         : {:f} deg'.format(self.waveDir*180./np.pi)
        return s

    #
    # public methods
    #
    def setup(self):
        if self.T is None:
            sys.exit('The wave time period must be defined when using MLER');
        if self.H is None:
            sys.exit('The wave height must be defined when using MLER');

        if self.H <= 0:
            sys.exit('Wave height (H) must be greater than zero to use MLER method')
        if self.T <= 0:
            sys.exit('Wave time period (T) must be greater than zero to use MLER method')

        self.dw = (self.endW - self.startW) / (self.numFreq-1)
        self.w  = np.linspace( self.startW, self.endW, self.numFreq )
            
        if self.type=='Bretschneider':
            self._BretschneiderSpectrum();
        else:
            sys.exit('Unknown spectrum type: {:s}'.format(self.type))

        self._waveNumber()

    def plotSpectrum(self,show=False):
        if self.w is None:
            sys.exit('Call waves.waveSetup before plotting spectrum');
        import matplotlib.pyplot as plt
        plt.figure()
        plt.plot(self.w,self.S)
        plt.title('{:s} spectrum for Hs = {:f} (m), Tp = {:f} (s)'.format(self.type,self.H,self.T))
        plt.xlabel('Frequency (rad/s)')
        plt.ylabel('Spectral amplitude (m^2)') #TODO: double check units
        if show is True: plt.show()

    #
    # protected methods
    #
    def _BretschneiderSpectrum(self):
        """ Calculate wave spectrum vector (obj.A)
        Used by wavesIrreg (wavesIrreg used by waveSetup)
        Sets self.S, self.A
        """
        freq = self.w / (2*np.pi)
        Tp = self.T
        Hs = self.H
        B = (1.057/Tp)**4
        A_irreg = B*(Hs/2)**2

        orig_settings = np.seterr(divide='ignore',invalid='ignore')
        S_f = (A_irreg*freq**(-5)*np.exp(-B*freq**(-4)))
        if np.isnan( S_f[0] ): # defined to avoid a NaN
            #print 'DEBUG: correcting NaN in spectrum'
            S_f[0] = 0.
        np.seterr(**orig_settings)
        assert( len(np.nonzero( np.isnan(S_f) )[0]) == 0 )

        Sf = S_f / (2*np.pi);
        self.S = Sf;
        self.A = 2 * Sf;

    def _waveNumber(self):
        """ Calculate wave number
        """
        #TODO: FIFTH-ORDER WAVENUMBER
        self.k = self.w**2 / self.g # deep water approximation
        if not self.deepWaterWave:
            lastk = self.k[:1]
            orig_settings = np.seterr(divide='ignore',invalid='ignore')
            for i in range(100):
                # TODO: more rigorous convergence
                self.k = self.w**2 / (self.g * np.tanh(self.k*self.waterDepth) );

                # check convergence
                #print i,np.max(np.abs(self.k[1:]-lastk))
                lastk = self.k[1:]
            np.seterr(**orig_settings)

        if np.isnan( self.k[0] ): # defined to avoid a NaN
            #print 'DEBUG: correcting NaN'
            self.k[0] = 0.
        assert( len(np.nonzero( np.isnan(self.k) )[0]) == 0 )

