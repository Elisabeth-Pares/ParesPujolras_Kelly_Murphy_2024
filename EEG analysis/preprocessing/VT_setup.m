function [exp, anPar, plt] = VT_setup(exp)

%% Path variables
exp.codePath = ['C:\Users\elisa\Desktop\VolatilityTask\Analysis\'];
exp.dataPath = ['E:\Elisabeth EEG Data\VolatilityTask\Data\']; %P01/S1/EEG ... // P01/S1/Behavioural ... //P01/S1/Eyelink ...
exp.modPath = ['E:\Elisabeth EEG Data\VolatilityTask\Modelling\']; 
exp.figPath = ['E:\Elisabeth EEG Data\VolatilityTask\Figures\']; 

%These are nested
exp.procBehPath = ['E:\Elisabeth EEG Data\VolatilityTask\Modelling\PreprocessedStim\']
exp.regPath = ['E:\Elisabeth EEG Data\VolatilityTask\Modelling\Regression results\'];

addpath(genpath(exp.dataPath));
addpath(genpath(exp.codePath));
addpath(genpath(exp.modPath));

%% =======================================================================
%% Experimental variables 
%% =======================================================================

exp.nSessions = 2; 
exp.nBlocksxSession = 9; 
exp.nEEGchans = 128; 

exp.nSamples = 10;
exp.H = 0.1;

exp.sampBounds = [1 2; 3 4; 5 6; 7 8; 9 10]; %To use during analysis, in VT_pullBehFits_optimised

%% EEG 
var = setup_EEGtriggers(exp); %Function in Task code.
exp.EEG.triggers = var.trigg;
exp.EEG.srate = 512;

exp.FrontChans = [71 72 80 81 93 94 103 70 73 79 82 92 95 102]; % at the very front edge SD can be high not necessarily because channels are bad but because person is blinking /moving eyes. We don't want to interpolate GOOD electrodes that later will HELP us pick up artifacts that we want to throw out.

exp.movCluster = [108,107,115,116,124,123;... %Left
    62,63,54,55,50,49]; %Right

exp.frontalCluster = [89,85,76,...
    88,86,75,...
    98,87,66];
exp.posteriorCluster = [4,19,20,21,22,23];

exp.frontBackCluster = [exp.frontalCluster, exp.posteriorCluster];

centralFocus = [86,75,64,53,51,4,113,114,109,88,...
    87,66,52,34,3,112,110,98,...
    65,33,2,111,97,...
    1];

parietalFocus = [4,32,31,30,21,17,18,5,...
    19,20];

exp.centroParietal = unique([centralFocus, parietalFocus]);

%% EL 
exp.EL.srate = 500;
      
exp.EL.sacc = 30; %Threshold (in pixels) to mark saccades
%% =======================================================================
%% Analysis parameters 
%% =======================================================================

%% EEG
anPar.EEG.highpassFilter = [0.1];
anPar.EEG.notchFilters = [50,100,150]; %EU power line

anPar.EEG.epochs = {'SL', 'RL'}; %Stimulus locked (trial onset); frame locked (sample onset); response-locked (action onset);
anPar.EEG.epochLims = [-0.5, 5.4; -2,0.5]; %For each of the epochs
anPar.EEG.epochTriggs = {exp.EEG.triggers.premask_on(1); {exp.EEG.triggers.resp_left(1), exp.EEG.triggers.resp_right(1)}}

%FieldTrip configuration
anPar.TFA             = [];
anPar.TFA.method      = 'mtmconvol';%multitaper method
anPar.TFA.taper       = 'hanning'; % with Hanning window
anPar.TFA.output      = 'pow';
anPar.TFA.precision   = 'single'; %saves disk space
anPar.TFA.channel     = [1:128];%[115,54,3];%[1:128];, 'Left', 'Right'alternative 'EEG', 'Parietal, where CPP is.
anPar.TFA.foi         = 1:1:35;   % frequencies of interest
anPar.TFA.toi         = 'all' %Time of interest
anPar.TFA.t_ftimwin   = repmat(0.4,1,length(anPar.TFA.foi));
anPar.TFA.keeptrials  = 'yes';
anPar.TFA.keeptapers  = 'no';
anPar.TFA.pad         = 'nextpow2';

%Frequencies to average over for regression
anPar.TFA.alpha       = [8:12];
anPar.TFA.beta        = [13:30];
anPar.TFA.gamma       = [30:60];
anPar.TFA.highGamma   = [60:90];

anPar.TFA.freqNames   = {'alpha', 'beta', 'gamma', 'highGamma'};
anPar.TFA.freqLabs    = {'Alpha (8-12Hz)', 'Beta (13-30Hz)', 'Gamma (30-60Hz)', 'High gamma (60-90Hz)'}


%% EL 
anPar.EL.baseline = 0.1*exp.EL.srate; %convert to sample space
anPar.EL.vars = {'gx', 'gy', 'pa'}; %variables that we want to extract from Eyelink data
anPar.EL.varLabs = {'X pos', 'Y pos', 'Pupil size'};
anPar.EL.eye = 1; %1 = left, 2 = right;

%% =======================================================================
%% Plotting parameters 
%% =======================================================================

plt = [];
plt.cpLabs = {'Early (<=3)', 'Mid (4-7)', 'Late (>=8)'}
plt.cpShade = {[1,1,4,4], [4,4,7,7], [7,7,10,10]}; %To indicate where the CP happened
plt.alpha = 0.15;
plt.topoTimes = [1,3,5,7,9];

plt.respLabs = {'Left', 'Right'};
plt.accLabs = {'Error', 'Correct'};

%Generated using "shades" in https://mycolor.space/
plt.newPalettes = {    {'#09325A','#395682', '#637DAC', '#8DA7D8','#B9D2FF'},... %Blue
    {'#09115A', '#413788', '#7161B8','#A28EEA','#D4BEFF'},...%Purple
    {'#5A0909','#88352D','#B85E53','#E9897B','#FFB6A7'},... %Red
    {'#5A3909','#815B2B','#A97F4D','#D3A571','#FECD97'},... %Orange
    {'#115A09','#397C2C','#5D9F4D','#81C46F','#A6EA92'},...} %Green
    {'#561647', '#BD4ABF', '#DA6CC0', '#F196F3', '#F8C9E9'}} %Pink


plt.brightPalettes = {{'#0080BF','#00ACDF', '#55D0FF', '#7CE8FF','#CCF9FF'},... %Blue
      {'#DC1C13','#EA4C46','#F07470','#F1959B','#F6BDC0'},...} %Red
      {'#00671A','#04A12B','#00CC33','#0BF446','#AAAAAA'}} %Green

end
