#!/usr/bin/python
# TODO: put this in egg
import wave
import simulation

class focusedWave(object):
    """ Based on MLERclass.m
    """

    def __init__(self,H,T,numFreq):
        self.waves  = wave.wave(H,T,numFreq)
        self.sim    = simulation.simulation()

        self.desiredRespAmp = -1                # [-]  Desired response amplitude.  M_d in documentation.

    def __repr__(self):
        s = 'MLER focused wave (desired response amplitude= {:f})'.format(self.desiredRespAmp)
        return s

    #
    # public methods
    #
    def setup(self):
        self.sim.setup()
        self.waves.setup()

    #
    # protected methods
    #


