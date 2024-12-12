%% ========================================================================
% VOLATILITY TASK MASTER SCRIPT - EEG analysis 
%% ========================================================================
% Author: Elisabeth Parés-Pujolràs
% Date: 09/March/2022
% This script reproduces the analysis in this paper: 

% Parés-Pujolràs, E., Kelly, S. P., & Murphy, P. R. (2024). Dissociable 
% encoding of evolving beliefs and momentary belief updates in distinct neural decision signals. 
% bioRxiv, 2024-05.[https://www.biorxiv.org/content/10.1101/2024.05.15.594345v2]

%% Set up analysis parameters 
clear; close all; 
cd('C:\Users\elisa\Desktop\VolatilityTask\Analysis');
addpath(genpath('C:\Users\elisa\Desktop\VolatilityTask\Analysis\'));

exp.recordEEG = 0;
[exp, anPar, plt] = VT_setup(exp);

exp.sub_id = {'01', '02', '03', '04', '05', '06', '07', '08', '09', '10', '11', '12', '13', '14', '15', '16', '17', '18', '19', '20'};

%% ========================================================================
% 1. EEG preprocessing
%% ========================================================================
%% STEP 0: Convert data to EEGlab format
% Convert from .bdf to .set and .fdt; add metadata [trial & block numbers]
for sub = exp.sub_id([1:end])
    for sess = 1:exp.nSessions(1:end)
        VT_convertData(sub{1}, sess, exp)
    end
end

%% STEP 1: Filter data
% Applies high pass, low pass & notch filter to EEG data
for sub = exp.sub_id([1:end])
    for sess = 1:exp.nSessions(1:end)
        EEG = VT_filter(sub{1},sess, exp, anPar);
    end
end

%% STEP 2: Epoch, detect bad channels & interpolate, re-reference
%  Epoch data, interpolate bad channels & re-reference to common average 
for sub = exp.sub_id([1:end])
    for sess = 1:exp.nSessions(1:end)
        VT_epochData(sub{1},sess,exp,anPar);
    end
end
%% STEP 3A:
% Run ICA to remove EyeBlink components
for sub = exp.sub_id([1:end])
    for sess = 1:exp.nSessions(1:end)
        EEG = VT_runICA(sub{1},sess,exp, anPar)
    end
end
%% MANUALLY REMOVE EYEBLINK COMPONENTS %%
%% STEP 3B:
% 3.1. TBT interpolation of triasl where fluctuations >150uV in any channel
% Relying on TBT toolbox (https://github.com/mattansb/TBT)
% 3.2. Load EL data, append to EEG dataset & mark trials with blinks for
% rejection

for sub = exp.sub_id([1:end])
    for sess = 1:exp.nSessions(1:end)
        EEG = VT_cleanEpochs_correct(sub{1},sess,exp, anPar); % Threshold; remove epochs with > 150uV drifts in any channel
    end
end

%% STEP 4: CSD filtering
% Apply Common Source Density (CSD) to increase spatial resolution of EEG signals. 
% Relying on CSD toolbox (https://psychophysiology.cpmc.columbia.edu/software/csdtoolbox/)
par.ds = 1; %Downsample data
for sub = exp.sub_id([1:end])
    for sess = 1:exp.nSessions(1:end)
        VT_applyCSD(sub{1},sess, exp, 'SL',par)
    end
end

%% STEP 5:
% Decompose data  
par.epoch = 'SL';
par.baseline = 'preTrial'; %Options = 'preTrial'; to implement: 'preEpoch'
par.csdLab = 'csd_';
par.filter = 'f01';
par.debug = 0;

for sub = exp.sub_id(1:end)
    for sess = 1:2
        VT_tfDecomposition(sub{1},sess,exp,anPar, par); %Check paper for methods - hanning window & taper window
    end
end

%% ========================================================================
% 2. EEG regressions
%% ========================================================================
% Implements a variety of regression analysis. 
% To reproduce the anaylsis in the paper, set the following settings: 
% - Fig. 2D: par.regName: full; par.dataType = 'TFA'; par.chans = 'lateralisation'; 
% - Fig. 2E: par.regName: psi_deltaPsi; par.dataType = 'TFA'; par.chans = 'lateralisation'; 

% - Fig. 3B: par.regName: full; par.dataType = 'ERP'; par.chans = 'parietal'; 
% - Fig. 3C: par.regName: deltaPsi; par.dataType = 'ERP'; par.chans = 'parietal'; 
% Additional options for supplementary Fig. 2 include: surprise, LLR, 

par.regName = 'deltaPsi'%Which regression to run?
par.dataType = 'ERP';  %TFA (for MBL) or ERP (for CPP);
par.chans = 'parietal' %Lateralisation (for MBL) or parietal (for CPP); 
par.save_aR2 = 0;  %Do we want to save goodness of fit (adjRsq) for this regression?
par.saveData =1; 
par.ds = 1;
par.excludeSacc = 1; %Exclude samples with saccades?

par.lowFilt = 1; %Apply low-pass filter (<6Hz). Always set to 1: reduces noise in single-trial regressions. 
par.hpFilt = 0; %Apply high-pass filter (>1Hz). Set to 1 for residuals analysis to remove slow drift. 

par.epoch = 'SL';
par.baseline = 'none'; %Options = 'preTrial'; 
par.baselineMethod = 'allT_subtr'; %Options: allT_subtr.
par.baselineTime = [0.1];
par.csdLab = 'csd_';

par.db = 0; %DB transform data

%Which LLR & prior values to use in regression: 'rel' for relative, signed; 'abs' for absolute, unsigned.
if strcmp(par.dataType, 'ERP')
    par.llr = 'abs'; %For CPP analysis
elseif strcmp(par.dataType, 'TFA')
    par.llr = 'rel'; %For MBL analysis
end

par.avFrex = 1;
par.fit = 1; 


for sub = 1:length(exp.sub_id)  
    [sumData] = VT_fitRegression(sub,exp,par);
    if par.excludeSacc; percReject(sub,1) = sumData.percReject; end
end

%% Plot results
% Motor beta lateralisation:
VT_figure2(exp)
% CPP:
VT_figure3(exp)

%% Compare adjusted Rsquares
% Produces plots comparing fits of objective evidence (LLR),
% effective evidence (deltaPsi), and belief (updatedPsi).
par.dataType = 'ERP'; %ERP for CPP analysis, TFA for MBL analysis.
VT_adjRsquaredPlots(exp,par)

%% ========================================================================
% 3. EEG residuals analysis
%% ========================================================================
%% Get residuals from regression of interest 
% CPP: Get residuals from HP filtered CPP deltaPsi regression
% MBL: Get residuals from psi_deltaPsi regression 

par.reg = 'deltaPsi'
for sub = 1:length(exp.sub_id) 
    disp(['Loading P' num2str(sub) '...'])
    [residuals] = VT_getResiduals(sub, exp, par)
    allRes{sub} = residuals;
end

%% EEG TO BEHAVIOUR REGRESSIONS
% Tests effect of single-trial fluctuations (residuals) on choice behaviour. 
% For CPP: "CPPRes_choice_timeRes"
% Choice ~ LLR + LLR*surprise + LLR*uncertainty + LLR*CPPresidual 

% For MBL: "BetaRes_choice_timeRes"
% Choice ~ LLR + LLR*surprise + LLR*uncertainty + MBLresidual
par.thisModel = 'CPPRes_choice_timeRes' %Options: BetaRes_choice_timeRes; CPPRes_choice_timeRes; 
close all 
  
GA = [];
GA = VT_eeg2behRegression(exp, par)

%% Plot
par.res = {'timeResolved_frontPostCluster', 'timeResolved_frontalCluster', 'timeResolved_posteriorCluster'};
VT_figure4A(exp, par) %Plot CPP residuals regression 
VT_figure4B(exp, par) %Plot motorBeta residuals regression 

%% EEG to EEG fit 
% DeltaBeta  ~ DeltaPsi + DeltaPsi*CPPres
par.regName = 'effEv_effEvxCPPres'; %effEv_effEvxCPPres
par.thisModel = par.regName; 
par.db = 1;
par.dataType = 'TFA';
par.freqBand = 'beta';
par.baselineTime = [0.1];
par.fit = 1;
par.res = {'timeResolved_frontPostCluster', 'timeResolved_frontalCluster', 'timeResolved_posteriorCluster'};
  
for nsubj = 1:20
[dat] = VT_eeg2eegRegression(nsubj,exp,par)
end
%% Plot residuals analysis 
par.minS = 1; par.maxS = 10; %Samples to include: 1 to 10.
VT_figure4C(exp, par)
