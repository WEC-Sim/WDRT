%%

% Simulation time and setup -- most of this is from wecSim.m, or made up as
% I saw fit.
sim.dt=.25;
sim.startTime=-150;
sim.endTime=150;
sim.T0=(sim.startTime+sim.endTime)/2;
sim.maxIt = ceil((sim.endTime - sim.startTime) / sim.dt + 1);   % maximum timestep
sim.T0idx=ceil((sim.T0 - sim.startTime) / sim.dt)+1;              % index of T0

sim.waveAmpTime = zeros(sim.maxIt+1,2);
sim.numFreq = 500;
sim.startW = 0;
sim.endW = 2*pi/1.0;
sim.dW = (sim.endW - sim.startW ) / (sim.numFreq);
sim.startW=sim.dW;
sim.w = sim.startW:sim.dW:sim.endW;

sim.X0=0;
sim.startX=-600;    % meters
sim.endX=600;       % meters
sim.dX=5;
sim.maxIx= ceil((sim.endX - sim.startX)/sim.dX);
sim.X0idx=ceil((sim.X0 - sim.startX) / sim.dX)+1;

sim.rho=1025;   % kg / m^3 -- density of sea water
sim.g=9.801;    % m/s^2
sim.waveDir=0;    % head on only case.
sim.water_depth=70;

waves=waveClass('irregular');
waves.numFreq=sim.numFreq;
waves.spectrumType='BS';
if exist('Case','var')
    waves.H=Case.H;
    waves.T=Case.T;
    %MLER.CaseNum=Case.Num;
    DOF=Case.DOF;
else
    waves.H=9.0;           % Test case of Hs=6, Tp=12.3.  Should be typical Hs Tp ratio
    waves.T=15.1;           % seconds
end

%% Wave info for MLW, MLER, and MLWR methods
sim.MLER.MultFact=1.9;      % 1.9 for H100
sim.MLER.RespAmpDesired=[1 0 1 0 -1 0] * waves.H/2 * sim.MLER.MultFact;       % desired wave amplitude (by DOF, m, rad)
sim.MLER.WaveAmpDesired=[1 0 1 0 1 0] * waves.H/2 * sim.MLER.MultFact;       % desired wave amplitude (m)



% options for WEC-Sim imported pieces.
sim.nlHydro=0;      % no non-linear hydro
sim.ssCalc=0;       % no state space, do convolution.


% Use the rm3 body and setup -- modified from wecSim.m --> not actually
% used other than in setting up for the RAO calculation (not used now)
% % body(1)=bodyClass('/Users/aplatt/software-development/github/WEC-Sim/tutorials/RM3/hydroData/rm3.h5',1);
% % body(1).momOfInertia = [20907301,0,0; 0, 21306090.66, 4305; 0, 4305, 37085481.11];      % Note this is unused here.
% % body(1).mass = 'equilibrium';
% % body(1).geometryFile = '/Users/aplatt/software-development/github/WEC-Sim/tutorials/RM3/geometry/float.stl';
% % body(1).bodyNumber=1;
% % 
% % % verify things and read the H5 file contents in (wecSim.m)
% % body(1).checkinputs;
% % body(1).readH5File;
% % body(1).checkBemio;
% % 
% % % Non-linear hydro info -- we won't use this here.
% % body(1).bodyGeo(body(1).geometryFile);
% % 
% setup the waves.
waves.checkinputs;
waves.waveSetup(sim.w, sim.water_depth, 0, sim.dt, sim.maxIt, sim.g, sim.endTime);


% % %% Calculate the RAO
% % %RAO=RAO_calc_BEM(sim,body);
% % 


%% Read RAO from file
RAO=complex(zeros(sim.numFreq,6));
if exist('Case','var')
   if DOF==3
      RAOfile=['../' Case.Dir '/' 'RAO_heave_RM3float.dat'];
      RAO(:,3)=RAO_read_file(sim,RAOfile);
   elseif DOF==5
      RAOfile=['../' Case.Dir '/' 'RAO_pitch_RM3float.dat'];
      RAO(:,5)=RAO_read_file(sim,RAOfile);
   end
else
   if DOF==3
      RAO(:,3)=RAO_read_file(sim,'data/RAO_heave_RM3float.dat');
   elseif DOF==5
      RAO(:,5)=RAO_read_file(sim,'data/RAO_pitch_RM3float.dat');
   else
      RAO(:,3)=RAO_read_file(sim,'data/RAO_heave_RM3float.dat');
      RAO(:,5)=RAO_read_file(sim,'data/RAO_pitch_RM3float.dat');
   end
end


% Generate some info on the spectrum
if exist('SpectralInfo')
    clear SpectralInfo;
end
WaveSpectInfo=SpectralInfo(waves.A,waves.w,waves.dw);


        
%% Calculate the MLER method:
% get the MLER ECM spectrum
MLER=MLER_calc(RAO,sim,waves);
if exist('Case','var')
    MLER.SeaNum=Case.SeaNum;
    MLER.Config=Case.Config;
    MLER.RAOfile=RAOfile;
end

% MLER Waves:
MLER.x=zeros(sim.maxIx,1);
MLER.t=zeros(sim.maxIt,1);
MLER.waveAmpTime=zeros(sim.maxIx,sim.maxIt,6,2);
NormWaveTest.waveAmpTime=zeros(sim.maxIx,sim.maxIt,6,2);
MLER.S_rn=zeros(sim.numFreq,6);
rephase=zeros(sim.numFreq,1);
MLER.NegResponse=zeros(6,1);
for jj=1:6

    MLER.S_rn(:,jj) = waves.A .* (MLER.Coeff(:,jj) .* sim.MLER.RespAmpDesired(jj)) .^2;
    MLER.SpectInfo(jj)=SpectralInfo(squeeze(MLER.S_rn(:,jj)),MLER.w,MLER.dw);

    
    if sign(sim.MLER.RespAmpDesired(jj)) == -1
        MLER.NegResponse(jj) = 1;
    else
        MLER.NegResponse(jj) = 0;
    end
    if MLER.NegResponse(jj) ~= 0
        rephase(:) = pi;
        MLER.phase(:,jj) = MLER.phase(:,jj) + pi;
    end
    
    for tt=1:sim.maxIt+1;
        t = (tt-1)*sim.dt+sim.startTime;
        MLER.t(tt)=t;
        
        for xx=1:sim.maxIx+1;
            x = (xx-1)*sim.dX + sim.startX;
            MLER.x(xx)=x;
            
            
            % conditioned wave
            MLER.waveAmpTime(xx,tt,jj,1) = sum(sqrt(MLER.S_rn(:,jj) .* waves.dw) .* cos( MLER.w .* (t-sim.T0) - MLER.phase(:,jj) - waves.k .* (x - sim.X0)));   % rephase included in MLER.phase

            % Response calculation
            MLER.waveAmpTime(xx,tt,jj,2) = sum(sqrt(MLER.S_rn(:,jj) .* waves.dw) .* abs(RAO(:,jj)) .* cos( MLER.w .* (t-sim.T0) - waves.k .* (x - sim.X0) + rephase));
            
            % normal wave with phasing
            NormWaveTest.waveAmpTime(xx,tt,jj,1) = sum(sqrt(waves.A*waves.dw) .* cos( waves.w .* (t-sim.T0) - MLER.phase(:,jj) - waves.k .* (x - sim.X0)));     % rephase included in MLER.phase

            % normal wave response with phasing
            NormWaveTest.waveAmpTime(xx,tt,jj,2) = sum(sqrt(waves.A*waves.dw) .* abs(RAO(:,jj)) .* cos( waves.w .* (t-sim.T0) - waves.k .* (x - sim.X0) + rephase));

        end
    end
    
    % renormalization of wave amplitudes
    if MLER.NegResponse(jj) == 0
        tmpMaxAmp = abs(max(max(MLER.waveAmpTime(:,:,jj,1))));
    else
        display(sprintf('generating most negative response for DOF = %i',DOF))
        tmpMaxAmp = abs(min(min(MLER.waveAmpTime(:,:,jj,1))));
    end
    rescale(jj) = abs(sim.MLER.WaveAmpDesired(jj)) / tmpMaxAmp;
    MLER.waveAmpTime(:,:,jj,:) = MLER.waveAmpTime(:,:,jj,:) * rescale(jj);      % rescale the wave amplitude result
    MLER.S_rn(:,jj)=MLER.S_rn(:,jj) * rescale(jj)^2;                            % rescale the wave spectral amplitude coefficients
    
    
end


% normal wave with the calculated phasing
NormWaveTest.w=MLER.w;
NormWaveTest.x=MLER.x;
NormWaveTest.t=MLER.t;

MLER.waveAmp.Max = max(max(MLER.waveAmpTime(:,:,DOF,1)));
MLER.waveAmp.Min = min(min(MLER.waveAmpTime(:,:,DOF,1)));
MLER.RespAmp.Max = max(max(MLER.waveAmpTime(:,:,DOF,2)));
MLER.RespAmp.Min = min(min(MLER.waveAmpTime(:,:,DOF,2)));


%figure out probability of wave crest at max or higher
[junk idx]=histc(MLER.waveAmp.Max,WaveSpectInfo.AmpToExceed);
ProbExceed=trapz(WaveSpectInfo.ProbExceed(idx:end))/trapz(WaveSpectInfo.ProbExceed);
clear junk idx;

% summary information
display(' ');
display(sprintf('-------- MLER Summary -- DOF = %i --------------------------',DOF));
if isfield(MLER,'SeaNum')
    display(sprintf('  Sea state:  %d, %g (m), %g (s)',MLER.SeaNum,waves.H,waves.T));
    if isfield(MLER,'Config')
        display(sprintf('  Config:     %s',MLER.Config));
    end
    if isfield(MLER,'RAOfile')
        display(sprintf('  RAO:        %s',MLER.RAOfile));
    end
end
if DOF>=4
    unitDOF = 'rad';
else
    unitDOF = 'm';
end

display('--------------------------------------------------------');
display(strcat(sprintf('  Requested wave amplitude:   %6.3f   (',sim.MLER.WaveAmpDesired(DOF)),unitDOF,')'));
display(strcat(sprintf('  wave rescaling:             %6.3f   (',rescale(DOF)),'-)'));
display(strcat(sprintf('--> Max wave elevation:       %6.3f   (',MLER.waveAmp.Max),unitDOF,')'));
display(strcat(sprintf('--> Min wave elevation:       %6.3f   (',MLER.waveAmp.Min),unitDOF,')'));
display(strcat(sprintf('  Effective wave height:      %6.3f   (',MLER.waveAmp.Max-MLER.waveAmp.Min),unitDOF,')'));
display(sprintf('  Probability of crest of %6.3f m or higher:   %6.3g',MLER.waveAmp.Max,ProbExceed));
display('--------------------------------------------------------');
display(strcat(sprintf('  Requested WEC response:     %6.3f   (',sim.MLER.RespAmpDesired(DOF)),unitDOF,')'));
display(strcat(sprintf('--> Max WEC response:         %6.3f   (',MLER.RespAmp.Max),unitDOF,')'));
display(strcat(sprintf('--> Min WEC response:         %6.3f   (',MLER.RespAmp.Min),unitDOF,')'));
display(strcat(sprintf('  Effective response height:  %6.3f   (',MLER.RespAmp.Max-MLER.RespAmp.Min),unitDOF,')'));



display('>> MovieFrames=AnimateMLER(sim,MLER,DOF)');
display('>> MLER_export(waves,sim,MLER,DOF,ProbExceed,MovieFrames)');

MovieFrames=AnimateMLER(sim,MLER,DOF);
MLER_export(waves,sim,MLER,DOF,ProbExceed,MovieFrames);

%% cleanup
clear ii jj t tmp tmp1;
 
 
