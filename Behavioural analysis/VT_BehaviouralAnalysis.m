%% ========================================================================
% VOLATILITY TASK MASTER SCRIPT - BEHAVIOURAL ANALYSIS
%% ========================================================================
% Author: Elisabeth Parés-Pujolràs
% Date: 09/March/2022
% This script reproduces the behavioural analysis in this paper: 

% Parés-Pujolràs, E., Kelly, S. P., & Murphy, P. R. (2024). Dissociable 
% encoding of evolving beliefs and momentary belief updates in distinct neural decision signals. 
% bioRxiv, 2024-05.[https://www.biorxiv.org/content/10.1101/2024.05.15.594345v2]

%% ========================================================================
%% Setup parameters 
clear 
cd('C:\Users\elisa\Desktop\VolatilityTask\Analysis');
addpath(genpath('C:\Users\elisa\Desktop\VolatilityTask\Analysis\'));

exp.recordEEG = 0;
[exp, anPar, plt] = VT_setup(exp);

exp.sub_id = {'01', '02', '03', '04', '05', '06', '07', '08','09','10','11','12', '13', '14', '15', '16', '17', '18', '19', '20'};

%% ========================================================================
% BEHAVIOURAL: Extract & concatenate behavioural data for analysis on R
%% ========================================================================

dat = struct(); dat.allBehav = []; dat.sumBehav = [];
for sub = exp.sub_id(1:end)
    for sess = 1:exp.nSessions(1:end)
        dat = VT_concatBehav(sub{1}, sess, exp, dat);
    end
end

%Plot bar graph with accuracies
allBehaviour = array2table(dat.allBehav);
allBehaviour.Properties.VariableNames = {'ID', 'sess', 'block', 'trial', 'dist', 'resp', 'acc', 'rt', 'timestamp', 'sacc'};
writetable(allBehaviour,[exp.dataPath, 'VT_BehaviouralData.csv'])

%% ========================================================================
% BEHAVIOURAL MODELLING: 
%% ========================================================================
  
%% Glaze model 
for subj = exp.sub_id(1:end)
    [pm_fit,err] = Glaze_fit_fixed(exp,subj)
end

%% Perfect accumulation 
for subj = 1:length(exp.sub_id)
    [pm_fit,err] = Perf_fit_fixed(exp,subj)
end

%% Perfect accumulation w/non-absorbing bounds
for subj = 1:length(exp.sub_id)
    [pm_fit,err] = Perf_naBounds_fit_fixed(exp,subj)
end

%% Leaky accumulation 
for subj = 1:length(exp.sub_id)
    [pm_fit,err] = Leak_fit_fixed(exp,subj)
end

%% Compare normative model to data 
par.include = 1;
par.stype = 'pCP'; %Surprise type. Different ways of calculating how "surprising" a sample is, given previous evidence. In the paper main text, pCP was used. take that as default. 
par.models2fit = {'ideal', 'perf','fitted_glaze','fitted_perf','fitted_perf_naBounds', 'fitted_leak','last_sample', 'data'};%,,'norm_glaze', 'fitted_glaze'};
par.datFields = {'oLLR','oLLR_full','LLR', 'LPR','surprise','psi','choices','acc','fsmps'...
    'psi_final', 'LPR_final', 'psi_final_fsmps','LPR_final_fsmps', 'LPR_final_full', 'LLR_final_fsmps','LLR_full','LPR_full','final_fsmps','final_fsmps',...
    'surprise_full','psi_full','choices_full','fdist_full','acc_full','fsmps_full', ...
    'fCPpos', 'fCPpos_full',...
    'finalDist_full', ...
    'nCP', 'nCP_full'};

par.regression = 'Murphy';
par.changeTime = 'last';

GA = struct();
for m = 1:length(par.models2fit)
    par.thisModel = par.models2fit{m};
    for sub = 1:length(exp.sub_id)
        [GA] = VT_pullBehFits_optimised(sub,exp,par,GA);
    end
end

%% Plot behaviour
% Reproduces all plots and statistics related to behavioural kernels &
% model fits in Fig 1  and Fig S1
par.thisAcc = 'acc';
VT_plotBehaviour(exp, par,GA)