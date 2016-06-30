% This is a stripped down version of the simulationClass.m file from WEC-Sim.  It
% is designed to be compatible with that version.

classdef simulationClassMLER < handle
    properties (SetAccess = 'public', GetAccess = 'public')
        startTime       =   -150                % [s]       starting time
        endTime         =   150                 % [s]       ending time
        dT              =   1.0                 % [s]       timestep
        T0              =   0                   % [s]       Time of maximum event
        
        startX          =   -300                % [m]       start of simulation space
        endX            =   300                 % [m]       end of simulation space
        dX              =   1.0                 % [m]       timestep
        X0              =   0                   % [m]       Position of maximum event
    end

    
    properties (SetAccess = 'private', GetAccess = 'public')
        maxIT           =   []                  % [-]       index corresponding to last timestep
        maxIX           =   []                  % [-]       index corresponding to last spatial position
        X               =   []                  % [m]       array of spatial coordinates for simulation
        T               =   []                  % [s]       array of time coordinates for simulation
    end
    
    
    methods (Access = 'public')
        function obj=simulationClassMLER
        end
         
        function simSetup(obj)
            % setup the time part
            if obj.startTime >= obj.endTime
                error('The starting time of the simulation must occur before the ending time.');
            end
            if obj.T0 < obj.startTime || obj.T0 > obj.endTime
                error('The time of T0 must be between the start and end times of the simulation.');
            end
            if obj.dT <= 1e-3
                error('The timestep is too small. Use a value larger than 1e-3.')
            end
            
            obj.maxIT = ceil((obj.endTime - obj.startTime) / obj.dT + 1);       % maximum timestep index
            obj.T = linspace(obj.startTime,obj.endTime,obj.maxIT);
            
            obj.maxIX = ceil((obj.endX - obj.startX)/obj.dX + 1);
            obj.X = linspace(obj.startX,obj.endX,obj.maxIX);
        end
    end
end