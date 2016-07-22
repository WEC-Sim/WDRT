#!/usr/bin/python
import numpy as np

class simulation(object):
    """ Based on simulationClassMLER.m
    A stripped down version of the simulationClass.m file from WEC-Sim.
    """

    def __init__(self):
        # simulation parameters
        self.startTime  = -150.0        # [s]       Starting time
        self.endTime    = 150.0         # [s]       Ending time
        self.dT         = 1.0           # [s]       Time-step size
        self.T0         = 0.0           # [s]       Time of maximum event

        self.startX     = -300.0        # [m]       Start of simulation space
        self.endX       = 300.0         # [m]       End of simulation space
        self.dX         = 1.0           # [m]       Horiontal spacing
        self.X0         = 0.0           # [m]       Position of maximum event

        # calculated values
        self._maxIT     = []            # [-]       Index corresponding to last timestep
        self._maxIX     = []            # [-]       Index corresponding to last spatial position
        self._X         = []            # [m]       Array of spatial coordinates for simulation
        self._T         = []            # [s]       Array of time coordinates for simulation

    def __repr__(self):
        s = 'simulationClass'
        return s

    #
    # public methods
    #
    def setup(self):
        """ Set up domain
        Sets: self._maxIT, self._T
        Sets: self._maxIX, self._X
        """
        if self.startTime >= self.endTime:
            raise ValueError('The starting time of the simulation must occur before the ending time.')
        if self.T0 < self.startTime or self.T0 > self.endTime:
            raise ValueError('The time of T0 must be between the start and end times of the simulation.')
        if self.dT <= 1e-3:
            raise ValueError('The timestep is too small. Use a value larger than 1e-3.')
            
        self._maxIT = int(np.ceil( (self.endTime - self.startTime)/self.dT + 1 )) # maximum timestep index
        self._T     = np.linspace( self.startTime, self.endTime, self._maxIT )

        self._maxIX = int(np.ceil( (self.endX - self.startX)/self.dX + 1 ))
        self._X     = np.linspace( self.startX, self.endX, self._maxIX )

    #
    # protected methods
    #

