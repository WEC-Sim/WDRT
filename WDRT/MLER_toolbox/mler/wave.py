#!/usr/bin/python
import numpy as np

class wave(object):
    """ Based on waveClassMLER.m
    A stripped down version of the waveClass.m file from WEC-Sim.
    """

    def __init__(self, H=None, T=None, numFreq=10001):
        if numFreq <= 10:
            raise ValueError('numFreq must be larger than 10')

        # sea state definition
        self.H               = H                    # [m]       Hs: sea state wave height
        self.T               = T                    # [s]       Tp: sea state time period of waves
        self.numFreq         = numFreq              # [-]       Number of frequencies
        self.type            = 'Bretschneider'      # [-]       Spectrum type

        # focused wave parameters
        self.startW          = 0.                   # [rad/s]   Starting frequency
        self.endW            = 2*np.pi              # [rad/s]   Ending frequency
        self.waveDir         = 0.                   # [rad]     Wave direction
        self.deepWaterWave   = False                # [-]       Deep water or not, depending on input water depth

        # constants
        self.waterDepth      = 70.                  # [m]       Water depth
        self.g               = 9.801                # [m/s^2]   Gravitational constant

        # calculated variables
        self._w              = None                 # [rad/s]   Wave frequency vector
        self._dw             = None                 # [rad/s]   Frequency step size
        self._S              = None                 # [m^2]     Wave spectrum vector
        self._A              = None                 # [m^2]     2*(wave spectrum vector)
        self._k              = None                 # [rad^2/m] Wavenumber array

    @property
    def w(self): return self._w
    @property
    def dw(self): return self._dw
    @property
    def k(self): return self._k
    @property
    def S(self): return self._S
    @property
    def A(self): return self._A

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
            raise UnboundLocalError('The wave time period must be defined when using MLER');
        if self.H is None:
            raise UnboundLocalError('The wave height must be defined when using MLER');

        if self.H <= 0:
            raise ValueError('Wave height (H) must be greater than zero to use MLER method')
        if self.T <= 0:
            raise ValueError('Wave time period (T) must be greater than zero to use MLER method')

        self._dw = (self.endW - self.startW) / (self.numFreq-1)
        self._w  = np.linspace( self.startW, self.endW, self.numFreq )
            
        if self.type=='Bretschneider':
            self._BretschneiderSpectrum();
        else:
            raise NotImplementedError('Unknown spectrum type: {:s}'.format(self.type))

        self._waveNumber()

    def plotSpectrum(self,show=True):
        if self._w is None:
            raise UnboundLocalError('Call waves.waveSetup before plotting spectrum');
        import matplotlib.pyplot as plt
        plt.figure()
        plt.plot(self._w,self._S)
        plt.title('{:s} spectrum for Hs = {:f} (m), Tp = {:f} (s)'.format(self.type,self.H,self.T))
        plt.xlabel('Frequency (rad/s)')
        plt.ylabel('Spectral amplitude (m^2-s)')
        if show is True: plt.show()

    #
    # protected methods
    #
    def _BretschneiderSpectrum(self):
        """ Calculate wave spectrum vector (obj.A)
        Used by wavesIrreg (wavesIrreg used by waveSetup)
        Sets self._S, self._A
        """
        freq = self._w / (2*np.pi)
        Tp = self.T
        Hs = self.H
        B = (1.057/Tp)**4
        A_irreg = B*(Hs/2)**2

        orig_settings = np.seterr(divide='ignore',invalid='ignore')
        S_f = (A_irreg*freq**(-5)*np.exp(-B*freq**(-4)))
        if np.isnan( S_f[0] ): # defined to avoid a NaN
            S_f[0] = 0.
        np.seterr(**orig_settings)
        assert( len(np.nonzero( np.isnan(S_f) )[0]) == 0 )

        Sf = S_f / (2*np.pi);
        self._S = Sf;
        self._A = 2 * Sf;

    def _waveNumber(self):
        """ Calculate wave number
        Sets self._k
        """
        self._k = self._w**2 / self.g # deep water approximation
        if not self.deepWaterWave:
            lastk = self._k[:1]
            orig_settings = np.seterr(divide='ignore',invalid='ignore')
            for i in range(100):
                self._k = self._w**2 / (self.g * np.tanh(self._k*self.waterDepth) );

                # check convergence
                #print i,np.max(np.abs(self.k[1:]-lastk))
                lastk = self._k[1:]
            np.seterr(**orig_settings)

        if np.isnan( self._k[0] ): # defined to avoid a NaN
            self._k[0] = 0.
        assert( len(np.nonzero( np.isnan(self._k) )[0]) == 0 )

