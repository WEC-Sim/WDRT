classdef MLERclass < handle
    properties (SetAccess = 'public', GetAccess = 'public')
        waves=waveClassMLER;
        sim=simulationClassMLER;
        desiredRespAmp  =  []                   %  [-]      Desired response amplitude.  M_d in documentation.
    end
    
    properties (SetAccess = 'private', GetAccess = 'public')
        RAO             =   []                  % [-]       Complex RAO array  N x 6
        RAOdataReadIn   =   zeros(6,1)          % [-]       What RAO dimensions did we read in?
        RAOdataFileName =   {'' '' '' '' '' ''} % [-]       Name of the data file read in for RAO
        CoeffA_Rn       =   []                  % [ ]       MLER coefficients A_{R,n}
        S               =   []                  % [m^2]     New Wave spectrum vector
        A               =   []                  % [m^2]     2*(new wave spectrum vector)
        phase           =   0                   % [rad]     Wave phase
        Spect           =   SpectralInfoClass   % [-]       Spectral info
        waveHeightDesired = []                  % [m]       Height of wave desired if renormalizing amplitudes
        rescaleFact     =   []                  % [-]       Rescaling factor for renormalizing the amplitude
    end
    
    
    methods (Access = 'public')
        function obj=MLERclass
        end
        
        
        
        function RAOread(obj,DOFread,RAO_File_Name)
            % read in the RAO from the specified file and assign it to a
            % dimension
            if isempty(obj.RAO)
                obj.RAO=complex(zeros(obj.waves.numFreq,6));        % set the size of the RAO matrix
            end

            
            % make sure we have setup the waves info first.
            if isempty(obj.waves.w)
                error('Call waves.waveSetup before calling ReadRAO');
            end
            
            % Format of file to read in:
            % Column 1:    period in seconds
            % Column 2:    response amplitude (m/m for DOF 1-3; or radians/m for DOF 4-6)
            % Column 3:    response phase (radians)
            file=importdata(RAO_File_Name);
            
            SizeData=size(file.data);
            tmpRAO=zeros(SizeData(1),SizeData(2));
            
            % convert from period in seconds to frequency in rad/s
            tmpRAO(:,1)= 2*pi ./ file.data(:,1);
            tmpRAO(:,2)= file.data(:,2);
            tmpRAO(:,3)= file.data(:,3);
            
            tmpRAO=sortrows(tmpRAO,1);
            
            % Add at w=0, amp=0, phase=0
            %%% questionable if the amplitude should be 1 or 0.  If set to 1, we count
            %%% on the spectrum having nothing out there.  For dimensions
            %%% other than in heave, this should be valid.  Heave should be
            %%% set to 1.
            if DOFread == 3
                tmp=[0 1 0];
            else
                tmp=[0 0 0];
            end
            tmpRAO=[tmp; tmpRAO];
            
            % Now interpolate to find the values
            Amp = interp1(tmpRAO(:,1),tmpRAO(:,2),obj.waves.w,'pchip');
            Phase = interp1(tmpRAO(:,1),tmpRAO(:,3),obj.waves.w,'pchip');

            % create the complex value to return
            obj.RAO(:,DOFread) = Amp .* exp(sqrt(-1) .* Phase);
            
            % set flag so we know that this dimension was read in, and save
            % filename
            obj.RAOdataReadIn(DOFread)=1;
            obj.RAOdataFileName{DOFread} = RAO_File_Name;
        end
        
        function RAOplot(obj,DOFtoPlot)
            % make sure we have setup the waves info first.
            if obj.RAOdataReadIn(DOFtoPlot) == 0
                error('Call waves.waveSetup and ReadRAO before RAOplot(DOF)');
            end
            figure
            plot(obj.waves.w,abs(obj.RAO(:,DOFtoPlot)),'-b',obj.waves.w,angle(obj.RAO(:,DOFtoPlot)),'--r');
            title (['RAO for dimension ', num2str(DOFtoPlot), ' from file ',obj.RAOdataFileName{DOFtoPlot}],'interpreter','none');
            xlabel ('Frequency (rad/s)');
            if DOFtoPlot <=3
                ylabel ('Response amplitude (m/m) / Response phase (rad)');
            else
                ylabel ('Response amplitude (rad/m) / Response phase (rad)');
            end
        end
        
        
        
        function MLERcoeffsGen(obj,DOFtoCalc,respDesired)
            
            % make sure we have a scalar value for the desired response
            % amplitude.
            if size(respDesired) == [1 1]
                obj.desiredRespAmp = respDesired;
            else
                error('Desired response amplitude (respDesired) must be a real scalar value.');
            end
            % check that we asked for something non-zero
            if respDesired == 0
                error('Desired response amplitude (respDesired) should be non-zero.');
            end
            
            % This function takes the spectrum and RAO and calculates MLER (most likely extreme response) coefficients
            % from a given RAO and spectrum
            
            S_tmp=zeros(obj.waves.numFreq,1);
            obj.S=zeros(obj.waves.numFreq,1);
            obj.A=zeros(obj.waves.numFreq,1);
            obj.CoeffA_Rn=zeros(obj.waves.numFreq,1);
            obj.phase=zeros(obj.waves.numFreq,1);
            
            % calculate the RAO times sqrt of spectrum
            % note that we could define:  a_n=(waves.A*waves.dw).^0.5;
            S_tmp(:)=squeeze(abs(obj.RAO(:,DOFtoCalc))).*2 .* obj.waves.A;          % Response spectrum.

            % calculate spectral moments and other important spectral values.
            obj.Spect.SpectralInfo(S_tmp(:),obj.waves.w,obj.waves.dw);
           
            % calculate coefficient A_{R,n}
            obj.CoeffA_Rn(:) = abs(obj.RAO(:,DOFtoCalc)) .* sqrt(obj.waves.A*obj.waves.dw) .* ((obj.Spect.M2 - obj.waves.w .* obj.Spect.M1) + obj.Spect.wBar .* (obj.waves.w .* obj.Spect.M0 - obj.Spect.M1 )) / (obj.Spect.M0*obj.Spect.M2 - obj.Spect.M1^2);        % Drummen version.  Dietz has negative of this.
            
            % save the new spectral info to pass out
            obj.phase(:)=-unwrap(angle(obj.RAO(:,DOFtoCalc))); % Phase delay should be a positive number in this convention
            
            % for negative values of Amp, shift phase by pi and flip sign
            obj.phase(obj.CoeffA_Rn(:) < 0) = obj.phase(obj.CoeffA_Rn(:) < 0) - pi;         % for negative Amplitudes, add a pi phase shift
            obj.CoeffA_Rn(obj.CoeffA_Rn(:) < 0) = -obj.CoeffA_Rn(obj.CoeffA_Rn(:) < 0);             % then flip sign on negative Amplitudes
            
            obj.S = obj.waves.S .* obj.CoeffA_Rn(:) .^2 * obj.desiredRespAmp ^2;
            obj.A = obj.waves.A .* obj.CoeffA_Rn(:) .^2 * obj.desiredRespAmp ^2;
            
            % if the response amplitude we ask for is negative, we will add
            % a pi phase shift to the phase information.  This is because
            % the sign of obj.desiredRespAmp is lost in the squaring above.
            % Ordinarily this would be put into the final equation, but we
            % are shaping the wave information so that it is buried in the
            % new spectral information, S.
            if obj.desiredRespAmp < 0
                obj.phase = obj.phase + pi;
            end
    
        end
        
        
        function MLERwaveAmpNormalize(obj,peakHeightDesired)
            % Renormalize the wave amplitude to some desired height of the
            % incoming wave. Uses the peak height (peak to MSL) desired 
            % rather than the full range height (peak to trough).
            
            display(['Renormalizing wave peak height to ' num2str(peakHeightDesired) 'm.  May take some time depending on spatial and temporal resolution.']);
            
            % check that we asked for a positive wave amplitude
            if peakHeightDesired <=0
                error('Wave height desired during renormalization must be positive.')
            end
            
            tmpMaxAmp=MLERpeakvalue(obj);

            % renormalization of wave amplitudes
            obj.rescaleFact = abs(peakHeightDesired) / abs(tmpMaxAmp);
            obj.S=obj.S .* obj.rescaleFact^2;                            % rescale the wave spectral amplitude coefficients
            obj.A=obj.A .* obj.rescaleFact^2;                            % rescale the wave spectral amplitude coefficients
            display(['Rescaling by ' num2str(obj.rescaleFact)]);
            
        end
        
        
        function MLERexportWECSim(obj,FileNameWEC)
            % Export the coefficients to a file that WEC-Sim can read in
            
            fileID = fopen(FileNameWEC,'w');
            
            % note that:
            %   WaveElev = sum( sqrt(2*S * dw) * cos( -k*(x-X0) + w*(t-T0) + Phase) )
            
            Freq=obj.waves.w ./(2*pi);
            fprintf(fileID,'%8.6f      ',Freq);       % output in hertz
            fprintf(fileID,'\n');
            fprintf(fileID,'%8.6f      ',obj.S);
            fprintf(fileID,'\n');
            fprintf(fileID,'%8.6f      ',obj.phase);
            fprintf(fileID,'\n');
            
            
            fclose(fileID);
            

        end
        
        
        function MLERexportCoeffs(obj,FileNameCoeff)
            
            % Now export the coefficients to a file
            fileID = fopen(FileNameCoeff,'w');
            
            % Header info
            fprintf(fileID,'# MLER wave profile generated on %s\n',char(datetime('now')));
            fprintf(fileID,'#\n');
            fprintf(fileID,'#\n');
            fprintf(fileID,'### Setup \n');
            fprintf(fileID,'# X0:         %6.3f   (m, peak position) \n',obj.sim.X0);
%             fprintf(fileID,'# x_start:    %6.3f   (m, start position) \n',obj.sim.startX);
%             fprintf(fileID,'# Rho:        %6.3f   (kg/m^3, density of water) \n',obj.sim.rho);      %not needed!
            fprintf(fileID,'# g:          %6.3f   (m/s^2, gravity) \n',obj.waves.g);
            fprintf(fileID,'# depth:      %6.3f   (m, water depth) \n',obj.waves.waterDepth);
%             fprintf(fileID,'# startTime:  %6.3f   (s) \n',obj.sim.startTime);
%             fprintf(fileID,'# endTime:    %6.3f   (s) \n',obj.sim.endTime);
%             fprintf(fileID,'# dt:         %6.3f   (s) \n',obj.sim.dT);
            fprintf(fileID,'# T0:         %6.3f   (s, response peak time) \n',obj.sim.T0);
            fprintf(fileID,'#\n');
            fprintf(fileID,'### Wave info:\n');
            fprintf(fileID,'# NumFreq: %6f    (-, number of frequencies)\n',obj.waves.numFreq);
            fprintf(fileID,'# dW:      %8.5g  (rad/s, frequency spacing)\n',obj.waves.dw);
            fprintf(fileID,'# Hs:      %6.3f  (m, significant wave height) \n',obj.waves.H);
            fprintf(fileID,'# Tp:      %6.3f  (s, wave period) \n',obj.waves.T);
            fprintf(fileID,'#\n');
%             fprintf(fileID,'# Max wave elevation:  %6.3f  (m) \n',max(max(MLER.waveAmpTime(:,:,DOFexport,1))));
%             fprintf(fileID,'# Min wave elevation:  %6.3f  (m) \n',min(min(MLER.waveAmpTime(:,:,DOFexport,1))));
%             fprintf(fileID,'# Max WEC Response:    %6.3f  (m) \n',max(max(MLER.waveAmpTime(:,:,DOFexport,2))));
%             fprintf(fileID,'# Min WEC Response:    %6.3f  (m) \n',min(min(MLER.waveAmpTime(:,:,DOFexport,2))));
%             fprintf(fileID,'#\n');
%             fprintf(fileID,'# Probability of wave exceeding Max wave elevation (normal distribution):  %g\n',ProbExceed);
%             fprintf(fileID,'#\n');
            fprintf(fileID,'# Note: Phase is calculated at X0 and t0.  For starting at a different point along x or in time, phase must be adjusted.\n');
            fprintf(fileID,'#\n');
            fprintf(fileID,'#Form of equation for wave elevation:\n');
            fprintf(fileID,'#   WaveElev = sum( sqrt(2*SpectAmp * dw) * cos( -k*(x-X0) + w*(t-T0) + Phase) ) \n');
            fprintf(fileID,'#\n');
            fprintf(fileID,'#\n');
            fprintf(fileID,'#   frequency      SpectAmp      Phase     wavenumber\n');
            fprintf(fileID,'#   (rad/s)          (m^2)       (rad)      (rad^2/m)\n');
            Phase=  obj.phase(:) + obj.waves.w .* obj.sim.T0  - obj.waves.k(:) .* obj.sim.X0;  % note sign: overall exported phase is still backwards
            for ww=1:obj.waves.numFreq
                fprintf(fileID,'%8.6f   %12.8g   %12.8f   %12.8f\n',obj.waves.w(ww),obj.S(ww),Phase(ww),obj.waves.k(ww));
            end
            
            fclose(fileID);
        end
        
        
        function MLERexportWaveAmpTime(obj,FileNameWaveAmpTime,DOFexport)
            % Export the wave amplitude timeseries at X0 to a file
            
            fileID = fopen(FileNameWaveAmpTime,'w');
            
            display(['Exporting wave amplitude time series for DOF = ' num2str(DOFexport) ' at X0.']);
            
            % calculate the series
            waveAmpTime=zeros(obj.sim.maxIT,2);
            t=zeros(obj.sim.maxIT,1);
            x = obj.sim.X0;
            for tt=1:obj.sim.maxIT;
                t(tt) = (tt-1)*obj.sim.dT+obj.sim.startTime;
                
                % conditioned wave
                waveAmpTime(tt,1) = sum(sqrt(obj.A * obj.waves.dw) .* cos( obj.waves.w .* (t(tt)-obj.sim.T0) + obj.phase - obj.waves.k .* (x - obj.sim.X0)));
                
                % Response calculation
                waveAmpTime(tt,2) = sum(sqrt(obj.A .* obj.waves.dw) .* abs(obj.RAO(:,DOFexport)) .* cos( obj.waves.w .* (t(tt)-obj.sim.T0) - obj.waves.k .* (x - obj.sim.X0)));
                
            end

            
            % Header info
            fprintf(fileID,'# MLER wave profile generated on %s\n',char(datetime('now')));
            fprintf(fileID,'#\n');
            fprintf(fileID,'#\n');
            fprintf(fileID,'### Setup \n');
            fprintf(fileID,'# X0:         %6.3f   (m, peak position) \n',obj.sim.X0);
%             fprintf(fileID,'# x_start:    %6.3f   (m, start position) \n',obj.sim.startX);
%             fprintf(fileID,'# Rho:        %6.3f   (kg/m^3, density of water) \n',obj.sim.rho);      %not needed!
            fprintf(fileID,'# g:          %6.3f   (m/s^2, gravity) \n',obj.waves.g);
            fprintf(fileID,'# depth:      %6.3f   (m, water depth) \n',obj.waves.waterDepth);
            fprintf(fileID,'# startTime:  %6.3f   (s) \n',obj.sim.startTime);
            fprintf(fileID,'# endTime:    %6.3f   (s) \n',obj.sim.endTime);
            fprintf(fileID,'# dt:         %6.3f   (s) \n',obj.sim.dT);
            fprintf(fileID,'# T0:         %6.3f   (s, response peak time) \n',obj.sim.T0);
            fprintf(fileID,'#\n');
            fprintf(fileID,'### Wave info:\n');
            fprintf(fileID,'# NumFreq: %6f    (-, number of frequencies)\n',obj.waves.numFreq);
            fprintf(fileID,'# dW:      %8.5g  (rad/s, frequency spacing)\n',obj.waves.dw);
            fprintf(fileID,'# Hs:      %6.3f  (m, significant wave height) \n',obj.waves.H);
            fprintf(fileID,'# Tp:      %6.3f  (s, wave period) \n',obj.waves.T);
            fprintf(fileID,'#\n');
             fprintf(fileID,'# Max wave elevation:  %6.3f  (m) \n',max(max(waveAmpTime(:,1))));
             fprintf(fileID,'# Min wave elevation:  %6.3f  (m) \n',min(min(waveAmpTime(:,1))));
             fprintf(fileID,'# Max WEC Response:    %6.3f  (m) \n',max(max(waveAmpTime(:,2))));
             fprintf(fileID,'# Min WEC Response:    %6.3f  (m) \n',min(min(waveAmpTime(:,2))));
             fprintf(fileID,'#\n');
%             fprintf(fileID,'# Probability of wave exceeding Max wave elevation (normal distribution):  %g\n',ProbExceed);
%             fprintf(fileID,'#\n');
            fprintf(fileID,'# Note: Phase is calculated at X0 and t0.  For starting at a different point along x or in time, phase must be adjusted.\n');
            fprintf(fileID,'#\n');
            fprintf(fileID,'#Form of equation for wave elevation:\n');
            fprintf(fileID,'#   WaveElev = sum( sqrt(2*SpectAmp * dw) * cos( -k*(x-X0) + w*(t-T0) + Phase) ) \n');
            fprintf(fileID,'#\n');
            fprintf(fileID,'#\n');            fprintf(fileID,'#   time      WaveHeight      LinearResp\n');
            if DOFexport==5
                fprintf(fileID,'#   (s)          (m)               (m)\n');
            else
                fprintf(fileID,'#   (s)          (m)              (rad)\n');
            end
            for tt=1:obj.sim.maxIT
                fprintf(fileID,'%12.8f   %12.8f   %12.8f\n',t(tt),waveAmpTime(tt,1),waveAmpTime(tt,2));
            end
            
            fclose(fileID);
            
            
            
        end
        
        
        function MovieFrames=MLERanimate(obj,DOF)
            % Animate the MLER results so that I can see what is happening.
            
            display(['Generating animation of wave profile and response for DOF = ' num2str(DOF) '.']);
            % create the 2D dataset
            waveAmpTime=zeros(obj.sim.maxIX,obj.sim.maxIT,2);
            t=zeros(obj.sim.maxIT,1);
            x=zeros(obj.sim.maxIX,1);
            for tt=1:obj.sim.maxIT;
                t(tt) = (tt-1)*obj.sim.dT+obj.sim.startTime;

                for xx=1:obj.sim.maxIX+1;
                    x(xx) = (xx-1)*obj.sim.dX + obj.sim.startX;
                    
                    % conditioned wave
                    waveAmpTime(xx,tt,1) = sum(sqrt(obj.A * obj.waves.dw) .* cos( obj.waves.w .* (t(tt)-obj.sim.T0) + obj.phase - obj.waves.k .* (x(xx) - obj.sim.X0)));
                    
                    % Response calculation
                    waveAmpTime(xx,tt,2) = sum(sqrt(obj.A .* obj.waves.dw) .* abs(obj.RAO(:,DOF)) .* cos( obj.waves.w .* (t(tt)-obj.sim.T0) - obj.waves.k .* (x(xx) - obj.sim.X0)));
                end
            end

            
            figure
            maxval=    max(max(max(waveAmpTime(:,:,1))),max(max(waveAmpTime(:,:,2))));
            minval=abs(min(min(min(waveAmpTime(:,:,1))),min(min(waveAmpTime(:,:,2)))));
            axis([x(1),x(obj.sim.maxIX),-1.1*max(maxval,minval),1.1*max(maxval,minval)])
            axis manual
            ax = gca;
            ax.NextPlot = 'replaceChildren';
            
            % plot first timestep
            plot(x,waveAmpTime(:,1,1),[obj.sim.X0 obj.sim.X0],get(ax,'YLim'),x,waveAmpTime(:,1,2))
            
            % labels
            xlabel(ax,'x (m)');
            if DOF==3
                title(ax,'MLER heave');
                ylabel(ax,'Amplitude (m)');
            elseif DOF==5
                title(ax,'MLER pitch');
                ylabel(ax,'Amplitude (m), Pitch (rad)');
            end
            %legend(ax,'MLER wave'); %,'MLER Response');
            
            aspect=daspect;
            
            % place a rectangle for the object on the plot
            dimX=0.01 * (x(obj.sim.maxIX) - x(1));
            dimY=dimX * aspect(2) / aspect(1);
            tmpX0val=interp1(x,waveAmpTime(:,1,2),obj.sim.X0);
            dim=[(obj.sim.X0-dimX) (tmpX0val-dimY) 2*dimX 2*dimY];
            float=rectangle('Position',dim,'EdgeColor','r');
            
            % make horizontal line for the extreme values
            Extremes=zeros(2,1);
            floatMin=line([-1.5*dimX 1.5*dimX],[Extremes(1) Extremes(1)],'Marker','.','LineStyle','-');
            floatMax=line([-1.5*dimX 1.5*dimX],[Extremes(2) Extremes(2)],'Marker','.','LineStyle','-');
            
            % plot envelope of the extreme wave profile
            WaveExtremes.Min=zeros(obj.sim.maxIX+1,1);
            WaveExtremes.Max=zeros(obj.sim.maxIX+1,1);
            ResponseExtremes.Min=zeros(obj.sim.maxIX+1,1);
            ResponseExtremes.Max=zeros(obj.sim.maxIX+1,1);
            
            % put timestamp on plot
            str = sprintf('t = %5.2f s',0);
            timestamp=annotation('textbox',[.72 .86 .17 .05],'String',str,'FitBoxToText','off','FontSize',14,'Margin',8);
            
            loops = obj.sim.maxIT;
            MovieFrames(loops) = struct('cdata',[],'colormap',[]);
            for tt = 1:loops
                
                % draw lines for extreme values of wave elevation
                WaveExtremes.Max(:)=max(WaveExtremes.Max,squeeze(waveAmpTime(:,tt,1)));
                WaveExtremes.Min(:)=min(WaveExtremes.Min,squeeze(waveAmpTime(:,tt,1)));
                
                % draw lines for extreme values of response
                ResponseExtremes.Max(:)=max(ResponseExtremes.Max,squeeze(waveAmpTime(:,tt,2)));
                ResponseExtremes.Min(:)=min(ResponseExtremes.Min,squeeze(waveAmpTime(:,tt,2)));
                
                % plot the wave profile along x
                plot(x,WaveExtremes.Max,'--b',x,WaveExtremes.Min,'--b',x,waveAmpTime(:,tt,1),'-b',x,ResponseExtremes.Max,':m',x,ResponseExtremes.Min,':m',x,waveAmpTime(:,tt,2),'--m');
                
                str = sprintf('t = %5.2f s',t(tt));
                delete(timestamp);
                timestamp=annotation('textbox',[.72 .86 .17 .05],'String',str,'FitBoxToText','off','FontSize',14,'Margin',8);
                
                % make a vertical line at X0
                line([obj.sim.X0 obj.sim.X0],get(ax,'YLim'),'LineStyle','-','Color','r');
                
                % make a box for representing the object
                tmpX0val=interp1(x,waveAmpTime(:,tt,2),obj.sim.X0);
                dim=[(obj.sim.X0-dimX) (tmpX0val-dimY) 2*dimX 2*dimY];
                float=rectangle('Position',dim,'EdgeColor','r');
                
                % draw lines for extreme values of float travel
                Extremes(2)=max(tmpX0val,Extremes(2));
                Extremes(1)=min(tmpX0val,Extremes(1));
                %delete(floatMin);
                %delete(floatMax);
                floatMin=line([-1.5*dimX 1.5*dimX],[Extremes(1) Extremes(1)],'Marker','.','LineStyle','-');
                floatMax=line([-1.5*dimX 1.5*dimX],[Extremes(2) Extremes(2)],'Marker','.','LineStyle','-');
                
                % make sure things get updated
                drawnow
                MovieFrames(tt) = getframe(gcf);
            end
            

        end

        
        function MLERexportMovie(obj,FileNameMovie,MovieFrames)
            % Now to save out the movie animation
            vidName=strcat(FileNameMovie,'.mp4');
            display(['Exporting video ' vidName]);
            myVideo=VideoWriter(vidName,'MPEG-4');
            myVideo.FrameRate=15;
            open(myVideo);
            writeVideo(myVideo,MovieFrames);
            close(myVideo);
        end
        
        
        function MLERplotT(obj,x,DOF)
            % create a time series plot at location X
            
            % add some error checking
%             if isempty(obj.w)
%                 error('Call waves.waveSetup before plotting spectrum');
%             end
            % calculate the series
            waveAmpTime=zeros(obj.sim.maxIT,2);
            t=zeros(obj.sim.maxIT,1);
            for tt=1:obj.sim.maxIT;
                t(tt) = (tt-1)*obj.sim.dT+obj.sim.startTime;
                % conditioned wave
                waveAmpTime(tt,1) = sum(sqrt(obj.A * obj.waves.dw) .* cos( obj.waves.w .* (t(tt)-obj.sim.T0) + obj.phase - obj.waves.k .* (x - obj.sim.X0)));
                % Response calculation
                waveAmpTime(tt,2) = sum(sqrt(obj.A .* obj.waves.dw) .* abs(obj.RAO(:,DOF)) .* cos( obj.waves.w .* (t(tt)-obj.sim.T0) - obj.waves.k .* (x - obj.sim.X0)));
            end
            figure
            plot(t,waveAmpTime(:,1),'-b',t,waveAmpTime(:,2),'--r')
            title ([ 'MLER response for DOF = ' num2str(DOF)]);
            xlabel ('Time (s)');
            if DOF <= 3
                ylabel ('Wave height (m), response (m)');
            else
                ylabel ('Wave height (m), response (rad)');
            end
        end
    end
    
    
    methods (Access = 'private')
        function peakval=MLERpeakvalue(obj)
            
            % the maximum may not occur at X0 or T0.  So, we have to
            % generate the entire time and space array, then find the
            % maximum and minimum.
            waveAmpTime=zeros(obj.sim.maxIX,obj.sim.maxIT);
            for tt=1:obj.sim.maxIT+1;
                t = (tt-1)*obj.sim.dT+obj.sim.startTime;
                
                for xx=1:obj.sim.maxIX+1;
                    x = (xx-1)*obj.sim.dX + obj.sim.startX;
                    
                    % conditioned wave
                    waveAmpTime(xx,tt) = sum(sqrt(obj.A * obj.waves.dw) .* cos( obj.waves.w .* (t-obj.sim.T0) + obj.phase - obj.waves.k .* (x - obj.sim.X0)));
                    
                    %                     % Response calculation
                    %                     waveAmpTime(xx,tt) = sum(sqrt(obj.S(:) .* obj.waves.dw) .* abs(obj.RAO(:,jj)) .* cos( obj.waves.w .* (t-obj.sim.T0) - obj.waves.k .* (x - obj.sim.X0)));
                    
                end
            end
            peakval = max(max(abs(waveAmpTime(:,:))));
        end
    end
    
end


