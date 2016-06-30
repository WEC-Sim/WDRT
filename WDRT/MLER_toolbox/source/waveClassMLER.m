% This is a stripped down version of the waveClass.m file from WEC-Sim.  It
% is designed to be compatible with that version.

classdef waveClassMLER < handle
    properties (SetAccess = 'public', GetAccess = 'public')
        startW          =   0.0                 % [rad/s]   starting frequency
        endW            =   2*pi                % [rad/s]   ending frequency
        numFreq         =   10001               % [-]       Number of frequencies
        waveDir         =   0                   % [rad]     Wave direction
        T               =   'NOT DEFINED'       % [s]       Tp -- sea state time period of waves
        H               =   'NOT DEFINED'       % [m]       Hs -- sea state wave height
        waterDepth      =   70                  % [m]       water depth used
        g               =   9.801               % [kg/m^2]  Gravitational constant
    end

    
    properties (SetAccess = 'private', GetAccess = 'public')
        w               =   []                  % [rad/s]   Wave frequency vector
        dw              =   []                  % [rad/s]   frequency step size
        type            =   'Brett Schneider'   % [-]       Spectrum type
        S               =   []                  % [m^2]     Wave spectrum vector
        A               =   []                  % [m^2]     2*(wave spectrum vector)
        k               =   []                  % [rad^2/m] Wave Number array
        deepWaterWave   =   0                   % [-]       Deep water or not, depending on input water depth
     end
    
    
    methods
        function obj = set.numFreq(obj,val)
            if val <= 10
                error('numFreq must be larger than 10')
            end
            obj.numFreq = val;
        end
        
        function obj = set.H(obj,val)
            if val <= 0
                error('Wave height (H) must be greater than zero to use MLER method');
            end
            obj.H=val;
        end
        
        function obj = set.T(obj,val)
            if val <= 0
                error('Wave time period (T) must be greater than zero to use MLER method');
            end
            obj.T=val;
        end
    end
    
    methods (Access = 'public')
        function obj = waveClassMLER
        end

        function waveSetup(obj)
            if strcmp(obj.T,'NOT DEFINED')
                error('The wave time period must be defined when using MLER');
            end
            
            if strcmp(obj.H,'NOT DEFINED')
                error('The wave height must be defined when using MLER');
            end

            obj.dw= (obj.endW-obj.startW)/(obj.numFreq-1);
            obj.w = (linspace(obj.startW,obj.endW,obj.numFreq))';
            
            obj.BrettSchneiderSpectrum;
            obj.waveNumber;
        end
        
        function plotSpectrum(obj)
            if isempty(obj.w)
                error('Call waves.waveSetup before plotting spectrum');
            end
            figure
            plot(obj.w,obj.S)
            title ([ obj.type,' spectrum for Hs = ',num2str(obj.H),' (m), Tp = ' ,num2str(obj.T), ' (s)']);
            xlabel ('Frequency (rad/s)');
            ylabel ('Spectral amplitude (m^2)');
        end
        function plotWaveNumber(obj)
            if isempty(obj.w)
                error('Call waves.waveSetup before plotting the values of the wavenumber');
            end
            figure
            plot(obj.w,obj.k)
            title ([ obj.type,' wavenumber for Hs = ',num2str(obj.H),' (m), Tp = ' ,num2str(obj.T), ' (s)']);
            xlabel ('Frequency (rad/s)');
            ylabel ('wavenumber (rad^2/m)');
        end
    end
    
    methods (Access = 'protected')
        function waveNumber(obj)
            % Calculate wave number
            obj.k = obj.w.^2./obj.g;
            if obj.deepWaterWave == 0
                for i=1:100
                    obj.k = obj.w.^2./obj.g./tanh(obj.k.*obj.waterDepth);
                end
            end
            if isnan(obj.k(1))          % defined to avoid a NaN
                obj.k(1)=0;
            end
        end
        
        function BrettSchneiderSpectrum(obj)
            % Calculate wave spectrum vector (obj.A)
            % Used by wavesIrreg (wavesIrreg used by waveSetup)
            freq = obj.w/(2*pi);
            Tp = obj.T;
            Hs = obj.H;
            B = (1.057/Tp)^4;
            A_irreg = B*(Hs/2)^2;
            S_f = (A_irreg*freq.^(-5).*exp(-B*freq.^(-4)));
            if isnan( S_f(1) )            % defined to avoid a NaN
                S_f(1)=0;
            end
            Sf = S_f./(2*pi);
            obj.S = Sf;
            obj.A = 2 * Sf;
        end
    end
end