%% 2016.03.21m  -- This script was taken from /Users/aplatt/software-development/scratchspace/OMAE2016--20160107/ECM-Wave for testing purposes
%%
%% Config      Surge    Heave    Pitch          DIR
%    A           -        y        -         RegularWaves_heaveOnly
%    B           -        y        -         RegularWaves_heaveTuned
%    C           -        y        y         RegularWaves_PitchTranslation


%% Heave cases
%%clear; close all; Case.SeaNum=1;Case.DOF=3;Case.H=5;     Case.T=8.2;  Case.Config='C';     Case.Dir='RegularWaves_PitchTranslation';     TestCase; 
%%clear; close all; Case.SeaNum=2;Case.DOF=3;Case.H=7;     Case.T=11.7; Case.Config='C';     Case.Dir='RegularWaves_PitchTranslation';     TestCase; 
%clear; close all; Case.SeaNum=3;Case.DOF=3;Case.H=9;     Case.T=15.1; Case.Config='C';     Case.Dir='RegularWaves_PitchTranslation';     TestCase; 
%%clear; close all; Case.SeaNum=4;Case.DOF=3;Case.H=10.5;  Case.T=20.5; Case.Config='C';     Case.Dir='RegularWaves_PitchTranslation';     TestCase; 
%%clear; close all; Case.SeaNum=5;Case.DOF=3;Case.H=8;     Case.T=25.9; Case.Config='C';     Case.Dir='RegularWaves_PitchTranslation';     TestCase; 

%clear; close all; Case.SeaNum=1;Case.DOF=3;Case.H=5;     Case.T=8.2;  Case.Config='C';     Case.Dir='RAO_data';     TestCase; 
%clear; close all; Case.SeaNum=2;Case.DOF=3;Case.H=7;     Case.T=11.7; Case.Config='C';     Case.Dir='RAO_data';     TestCase; 
clear; close all; Case.SeaNum=3;Case.DOF=3;Case.H=9;     Case.T=15.1; Case.Config='C';     Case.Dir='RAO_data';     TestCase; 
%clear; close all; Case.SeaNum=4;Case.DOF=3;Case.H=10.5;  Case.T=20.5; Case.Config='C';     Case.Dir='RAO_data';     TestCase; 
%clear; close all; Case.SeaNum=5;Case.DOF=3;Case.H=8;     Case.T=25.9; Case.Config='C';     Case.Dir='RAO_data';     TestCase; 


%% Pitch cases
%%clear; close all; Case.SeaNum=1;Case.DOF=5;Case.H=5;     Case.T=8.2;  Case.Config='C';     Case.Dir='RegularWaves_PitchTranslation';     TestCase; 
%%clear; close all; Case.SeaNum=2;Case.DOF=5;Case.H=7;     Case.T=11.7; Case.Config='C';     Case.Dir='RegularWaves_PitchTranslation';     TestCase; 
%clear; close all; Case.SeaNum=3;Case.DOF=5;Case.H=9;     Case.T=15.1; Case.Config='C';     Case.Dir='RegularWaves_PitchTranslation';     TestCase; 
%%clear; close all; Case.SeaNum=4;Case.DOF=5;Case.H=10.5;  Case.T=20.5; Case.Config='C';     Case.Dir='RegularWaves_PitchTranslation';     TestCase; 
%%clear; close all; Case.SeaNum=5;Case.DOF=5;Case.H=8;     Case.T=25.9; Case.Config='C';     Case.Dir='RegularWaves_PitchTranslation';     TestCase; 

%clear; close all; Case.SeaNum=1;Case.DOF=5;Case.H=5;     Case.T=8.2;  Case.Config='C';     Case.Dir='RAO_data';     TestCase; 
%clear; close all; Case.SeaNum=2;Case.DOF=5;Case.H=7;     Case.T=11.7; Case.Config='C';     Case.Dir='RAO_data';     TestCase; 
clear; close all; Case.SeaNum=3;Case.DOF=5;Case.H=9;     Case.T=15.1; Case.Config='C';     Case.Dir='RAO_data';     TestCase; 
%clear; close all; Case.SeaNum=4;Case.DOF=5;Case.H=10.5;  Case.T=20.5; Case.Config='C';     Case.Dir='RAO_data';     TestCase; 
%clear; close all; Case.SeaNum=5;Case.DOF=5;Case.H=8;     Case.T=25.9; Case.Config='C';     Case.Dir='RAO_data';     TestCase; 


