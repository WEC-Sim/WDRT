classdef SpectralInfoClass < handle


    properties (SetAccess = 'public', GetAccess = 'public')
    end
    
    properties (SetAccess = 'private', GetAccess = 'public')
        M0              =   []                  %   Zeroeth spectral moment
        M1              =   []                  %   First spectral moment
        M2              =   []                  %   Second spectral moment
        M3              =   []                  %   Third spectral moment
        M4              =   []                  %   Fourth spectral moment
        AmpToExceed     =   0:.01:25            %   [m]     amplitude to exceed
        BandCoeff       =   0                   %   width of spectral band: narrowband < 0.5 <= broadband
        ProbExceed      =   []                  %   [-]     probability of wave exceeding height corresponding to AmpToExceed
        Hs              =   0                   %   [m]     calculated spectral amplitude --> normal distributions only
        wBar            =   0                   %   [rad/s] central frequency            
    end
    
    
    methods (Access = 'public')
        function obj=SpectralInfoClass
        end
        
        function obj=SpectralInfo(obj,S,w,dw)
            % S   ->  [m^2]     Wave spectrum vector
            % w   ->  [rad/s]   Wave frequency vector
            % dw  ->  [rad/s]   frequency step size
            
            obj.M0=trapz( S .* w.^0 * dw );
            obj.M1=trapz( S .* w.^1 * dw );
            obj.M2=trapz( S .* w.^2 * dw );
            obj.M3=trapz( S .* w.^3 * dw );
            obj.M4=trapz( S .* w.^4 * dw );
            
            % Band coefficient:  narrowband < 0.5 <= broadband
            obj.BandCoeff = sqrt(1-obj.M2^2/(obj.M0*obj.M4));   % from http://ocw.mit.edu/courses/mechanical-engineering/2-019-design-of-ocean-systems-spring-2011/lecture-notes/MIT2_019S11_OWE.pdf
            
            % Probability of amplitude exceeding some value
            obj.AmpToExceed=0:.01:25;
            obj.ProbExceed=2*sqrt(1-obj.BandCoeff^2)/(1+sqrt(1-obj.BandCoeff^2)).*exp(-obj.AmpToExceed.^2/(2*obj.M0));  % from http://ocw.mit.edu/courses/mechanical-engineering/2-019-design-of-ocean-systems-spring-2011/lecture-notes/MIT2_019S11_OWE.pdf
            
            
            % Hs -- calculated from M0 (also written as H_{1/3})
            obj.Hs=4*sqrt(obj.M0);
            
            % wBar -- the mean spectral frequency.
            obj.wBar=obj.M1/obj.M0;
            
        end
        
        
    end
    
end