#!/usr/bin/python
# TODO: put this in egg
import sys
import wave
import simulation
import numpy as np
import scipy.interpolate

class focusedWave(object):
    """ Based on MLERclass.m
    """

    def __init__(self,H,T,numFreq):
        self.waves  = wave.wave(H,T,numFreq)
        self.sim    = simulation.simulation()

        self.desiredRespAmp = -1                            # [-]   Desired response amplitude.  M_d in documentation.

        # TODO: private set, public get
        self.RAO                = None                      # [-]   Complex RAO array  N x 6
        self.RAOdataReadIn      = np.zeros(6,dtype=bool)    # [-]   What RAO dimensions did we read in?
        self.RAOdataFileName    = ['','','','','','']       # [-]   Name of the data file read in for RAO
        self.CoeffA_Rn          = None                      # [ ]   MLER coefficients A_{R,n}
        self.S                  = None                      # [m^2] New Wave spectrum vector
        self.A                  = None                      # [m^2] 2*(new wave spectrum vector)
        self.phase              = 0.0                       # [rad] Wave phase
       #self.Spect              = spectralInfo.spectrum()   # [-]   Spectral info
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
        DOFread : 1 - 3 (translation DOFs)
                  4 - 6 (rotation DOFs)
        """
        if self.RAO is None:
            self.RAO = np.zeros( (self.waves.numFreq,6), dtype=complex ) # set the size of the RAO matrix

        if self.RAOdataReadIn[DOFread] is True:
            print 'WARNING: RAO dof=',DOFread,'already read from',self.RAOdataFileName[DOFread]
        
        # make sure we have setup the waves info first.
        if self.waves.w is None:
            sys.exit('Call waves.waveSetup before calling ReadRAO')
        
        # Format of file to read in:
        # Column 1:    period in seconds
        # Column 2:    response amplitude (m/m for DOF 1-3; or radians/m for DOF 4-6)
        # Column 3:    response phase (radians)
        #
        print 'Reading RAO dof=',DOFread,'from',RAO_File_Name
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
        tmpRAO = tmpRAO[ np.argsort(tmpRAO[:,0]) ] # sort by frequency
        
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
        

    #
    # protected methods
    #


