#!/usr/bin/env python
import numpy as np

class stats(object):
    """ Based on SpectralInfoClass.m
    """
    def __init__(self,S=None,w=None,dw=None):
        self.M0             = None      # [(response units)^2*(rad/s)^0] Zeroth spectral moment
        self.M1             = None      # [(response units)^2*(rad/s)^1] First spectral moment
        self.M2             = None      # [(response units)^2*(rad/s)^2] Second spectral moment
        self.M3             = None      # [(response units)^2*(rad/s)^3] Third spectral moment
        self.M4             = None      # [(response units)^2*(rad/s)^4] Fourth spectral moment

        self.AmpToExceed    = None      # [m]       Vector of amplitudes to exceed
        self.BandCoeff      = 0.0       # [-]       Width of spectral band (calculated): narrowband < 0.5 <= broadband
        self.ProbExceed     = None      # [-]       Probability of wave exceeding height corresponding to AmpToExceed

        self.Hs             = 0         # [m]       Calculated spectral amplitude --> normal distributions only
        self.wBar           = 0         # [rad/s]   Central frequency            

        if not (S is None or w is None or dw is None):
            self.calculate(S,w,dw)

    def calculate(self,S,w,dw,ampMax=25.,ampStep=0.01):
        """ S   ->  [(response units)^2-s/rad]   Wave spectrum vector
            w   ->  [rad/s]   Wave frequency vector
            dw  ->  [rad/s]   Frequency step size
        """
        self.M0 = np.trapz( S * w**0 * dw )
        self.M1 = np.trapz( S * w**1 * dw )
        self.M2 = np.trapz( S * w**2 * dw )
        self.M3 = np.trapz( S * w**3 * dw ) # (not used)
        self.M4 = np.trapz( S * w**4 * dw )
        
        # Band coefficient:  narrowband < 0.5 <= broadband
        # compare bandwidth parameter with Dietz2004 eqn 3.41
        self.BandCoeff = np.sqrt(1.-self.M2**2/(self.M0*self.M4)) # from http://ocw.mit.edu/courses/mechanical-engineering/2-019-design-of-ocean-systems-spring-2011/lecture-notes/MIT2_019S11_OWE.pdf (AP)
        # Probability of amplitude exceeding some value
        self.AmpToExceed = np.linspace(0,ampMax,ampMax/ampStep+1)
        self.ProbExceed = 2.0*np.sqrt(1.0-self.BandCoeff**2) / ( 1.0 + np.sqrt(1.0-self.BandCoeff**2) ) \
            * np.exp( -self.AmpToExceed**2 / (2.0*self.M0) ) # from http://ocw.mit.edu/courses/mechanical-engineering/2-019-design-of-ocean-systems-spring-2011/lecture-notes/MIT2_019S11_OWE.pdf (AP)
        
        # Hs -- calculated from M0 (also written as H_{1/3})
        self.Hs = 4*np.sqrt(self.M0)
        
        # wBar -- the mean spectral frequency.
        self.wBar = self.M1 / self.M0
        
    

