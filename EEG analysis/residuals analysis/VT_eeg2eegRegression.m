function [dat, allVif] = VT_eeg2eegRegression(nsubj,exp,par,allVif) %Adapted from runSRanalysis_LIhalf_ext

allsubj = exp.sub_id; % ORDER IS IMPORTANT FOR TRIAL-IDs!!!!
sub = allsubj{nsubj};

subj = ['P' allsubj{nsubj}];%find(strcmp(allsubj,subj));

% Paths dervied from processing options
smpwin = [0 1.4];  % window for sample-wise analyses
trlwin = [-0.1 5.4];  % window for full-trial-wise analyses

% ==================================================================
% LOAD DATA
% ==================================================================
% loop through sessions and load data
trl_data=[];
LLR_full=[]; LPR_full=[]; surprise_full=[]; psi_full=[];
LLR_all=[]; LPR_all=[]; surprise_all=[]; psi_all=[];

choices_full=[]; sess_full=[]; allT_sess_data = [];
sumData_all = []; sumData_full = [];
latData_all = []; latData_full =[];
l_data_full = []; r_data_full = [];

figure

thisPath = ['P' sub '\S2\EEG\'];
fullPath = ([exp.dataPath, thisPath]);
semifullPath = ([exp.dataPath, 'P' sub '\']);
cd(fullPath);

   
for sess = 1:exp.nSessions %(strcmp(allsubj,sub))
    trialNum = [];
    thisPath = ['P' sub '\S' num2str(sess) '\EEG\'];
    fullPath = ([exp.dataPath, thisPath]);
    semifullPath = ([exp.dataPath, 'P' sub '\']);
    cd(fullPath);
    
    % Load preprocessed & epoched EEG data
    EEG = pop_loadset( [fullPath '\csd_DSacICASL_f01fc_P' sub '_S' num2str(sess) '.set'] );
    
    % If TFA, load TFAdata too
    if strcmp(par.freqBand, 'alphaBeta')
        filename = ['tfa_', sub '_', num2str(sess) '_' par.freqBand '.mat']
    else
        filename = ['tfa_', sub '_', num2str(sess) '.mat']
    end
    load(filename);
    
    % and load CPP residuals
    load([semifullPath, '\VT_CPPresiduals_deltaPsi.mat']);
    
    
    plot(nanmean(EEG.data(19,:,:),3)')
    xline([0.4*EEG.srate]); xline([0.5*EEG.srate]); sgtitle(['P' sub, ', S' num2str(sess)])
    hold on;
    
    %Permute order so it fits Peter's default
    sess_data = freqData; %dimensions are trials*freqs*times;
    
    % CHECK TRIAL ORDER
    for t = 1:size([EEG.epoch],2)
        trialNum(t)= ([EEG.epoch(t).eventntrial{1}]);
    end
    
    if ~issorted(trialNum)  % if trials aren't in ascending order, sort - this happens in participant 01, session 02
        [mne_ids,sortI] = sort(trialNum);
        sess_data = sess_data(sortI,:,:);
    end
    
    % If dB transform
    if par.db == 1 & ~strcmp(par.dataType, 'ERP')
        baselineTime = round([(0.5-par.baselineTime)*200]):round([0.5*200]);
        EEG.baselines= squeeze(mean(sess_data(:,:,baselineTime),3)); %Recalculate baselines
        sess_data = 10.*log10(sess_data./nanmean(EEG.baselines));  % convert to dB
    end
    
    %Keep all trials for plotting
    all_sess_data = sess_data;
    
    times = round(EEG.times)/1000;
    
    % trim data to desired epoch length and concatenate across sessions
    sess_data = sess_data(:,:,times>=trlwin(1) & times<=trlwin(2)); %sess_data = [trials x channels x timepoints]
    trl_data = cat(1,trl_data,sess_data);
    allT_sess_data = cat(1,allT_sess_data, all_sess_data);
    
    sess_full = [sess_full; repmat(sess, size(sess_data,1), 1)];
end
    

% load and concatenate behaviour/pupil
load(['C:\Users\elisa\Desktop\VolatilityTask\Analysis\ProcStim_fitted_glaze_' subj '_S2.mat']) % This MDLVARS. comes from here_
choices = dat.fitted_glaze.choices;
choices_full = dat.fitted_glaze.choices_full;
psi_full = dat.fitted_glaze.psi_full;
llr_full = dat.fitted_glaze.oLLR_full;
surprise_full = dat.fitted_glaze.surprise_full;

%Check that length of behaviour & EEG are the same
%assert(isempty(find((mne_ids'-trialID)~=0, 1)),'ERROR: Mismatch in eeg/behaviour trial IDs')
if strcmp(sub, '01'); idcs = [225:235];  psi_full(idcs,:) = []; llr_full(idcs,:) = []; surprise_full(idcs,:) = []; choices_full(idcs,:)=[];
    wholeIdx = [321:335]; dat.fitted_glaze.LLR(wholeIdx,:) = [];
end

% ==================================================================
% ENSURE USE OF DESIRED FORM OF PRIOR BELIEF
% ==================================================================
prior_full(:,:) = psi_full(:,2:end); 

% ==================================================================
% MAKE CATEGORICAL SESSION REGRESSORS
% ==================================================================
sessions = unique(sess_full);
sess_r = zeros(length(sess_full),length(sessions)-1);
for s = 1:length(sessions)-1
    sess_r(sess_full==sessions(s+1),s) = ones(length(find(sess_full==sessions(s+1))),1);
end


if ~exist('deltaLat', 'var')
    % EPP edited: pick ONLY trials with full sequence -otherwise regression
    % includes all trials, even those where sample times are out of trial
    % bounds.
    [ridx, cidx] = find(isnan(dat.fitted_glaze.LLR)); %Returns all columns without nans
    trials2exclude = unique(ridx);
    
    %Remove non-full trials;
    trl_data(trials2exclude,:,:) = [];
    sess_full(trials2exclude,:) = [];
    
    sess_r(trials2exclude,:) = [];
    assert(length(choices_full)==size(trl_data,1),'ERROR: Trial counts in eeg/behaviour are unequal')
    
    % ==================================================================
    % SEGMENT DATA AROUND INDIVIDUAL SAMPLES
    % ==================================================================
    fprintf('Concatenating trial- & sample-wise data segments...\n')
    times = times(times>=trlwin(1) & times<=trlwin(2));
    onsets = 0.4:0.4:0.4*10;  onsets=round(onsets,1);  % vector of all sample onset times relative to pre-mask; epoched around pre-mask on
    smptimes = times(times>=smpwin(1) & times<=smpwin(2));  % getting vector of sample times relative to dot onset
    
    for s = 1:length(onsets)
        stsmp = find(times>=onsets(s),1,'first');
        smp_data(:,:,:,s) = trl_data(:,:,stsmp:stsmp+length(smptimes)-2); %Smp_data dimensions are [trials x channels x timepoints x samples]
        smp_baselines(:,:,:,s) = nanmean(trl_data(:,:,stsmp-20:stsmp),3);
    end
    
    newSmpData =  nanmean(smp_data(:,exp.movCluster(2,:),:,:),2) - nanmean(smp_data(:,exp.movCluster(1,:),:,:),2); %Before, 1:6
    smp_data = newSmpData;
    % Baseline lateralisation %
    smp_data = smp_data - nanmean(nanmean(smp_baselines(:,exp.movCluster(2,:),1,1),2)-nanmean(smp_baselines(:,exp.movCluster(1,:),1,1),2));
    
    for s = 1:10
        if s == 10
            thisStsmp = find(times>=onsets(s),1,'first') + 0.4*200;
            thisDat = trl_data(:,:,thisStsmp:end); %Smp_data dimensions are [trials x channels x timepoints x samples]
            newThisDat=  nanmean(thisDat(:,exp.movCluster(2,:),:),2) - nanmean(thisDat(:,exp.movCluster(1,:),:),2); %Before, 1:6
            thisDat = newThisDat; 
            deltaLat(:,s) = squeeze(nanmean(thisDat(:,:,1:0.1*200),3))-squeeze(nanmean(smp_data(:,1,1:0.1*200,s),3));
        else
            deltaLat(:,s) = squeeze(nanmean(smp_data(:,1,1:0.1*200,s+1),3)) - squeeze(nanmean(smp_data(:,1,1:0.1*200,s),3));
            deltaLat_50ms(:,s) = squeeze(nanmean(smp_data(:,1,1:0.05*200,s+1),3)) - squeeze(nanmean(smp_data(:,1,1:0.05*200,s),3));
            deltaLat_2samp(:,s) = squeeze(nanmean(smp_data(:,1,1:0.1*200,s+1),3)) - squeeze(nanmean(smp_data(:,1,0.8*200:0.9*200,s),3));
        end
        %     deltaLat_n2Samps(:,s) = squeeze(nanmean(smp_data(:,1,1:0.1*200,s+1),3)) - squeeze(nanmean(smp_data(:,1,1:0.1*200,s),3));
    end
    save(['tfa_deltaLat_' subj '_' par.freqBand '.mat'], 'deltaLat', 'deltaLat_50ms', 'deltaLat_2samp', 'sess_r')
end



% Effective prior
effPrior = abs(diff(psi_full,1,2));

% Effective ev_signed
effPrior_signed = (diff(psi_full,1,2));


% Do for all residuals 
for r = 1:length(par.res)
    CPPresidual = []; 
    thisRes = par.res{r}; 
    CPPresidual = res.(thisRes);


if par.fit == 1
    
    % ==================================================================
    % RUN SINGLE-TRIAL REGRESSIONS
    % ==================================================================
    % Sample-wise regressions, sample*electrode*time
    
    %% Dummy code to check what happens during zscoring -abs
    fprintf('Running regressions...\n Sample ')
%     smp_data = squeeze(smp_data);
    for s = 1:10  % looping through samples
        fprintf('%d, ',s)
        for e = 1%1:128 %1:size(smp_data,2)  % looping through electrodes OR frequencies
            for t = 1:EEG.srate% Do for 1 second, not all! size(smp_data,3)  % looping through time-points
                if strcmp(par.regName, 'effEv_effEvxCPPres')
                    regNames = {'effEv', 'effEvxCPPres'};
                    regressors = nanzscore([(effPrior_signed(:,s)) nanzscore((effPrior_signed(:,s))).*nanzscore(squeeze(CPPresidual(s,t,:))) sess_r]);
                    m = regstats(nanzscore(deltaLat(:,s)),regressors,'linear',{'beta','tstat','rsquare'});  
                end
                
                % Save beta coefficients and residuals
                for r = 1:length(regNames)
                    thisReg = regNames{r};
                    dat.(thisReg).(thisRes).beta(e,t,s) = m.beta(r+1); %Beta coefficints
                end
                
            end
        end
    end
    fprintf('Done.\n')
    
    thisPath = ['P' sub];
    fullPath = ([exp.dataPath, thisPath]);
    
end
end
save([fullPath, '\VT_cICAEEGreg_allCh_CPPRes2Beta_' par.regName '_allRes_P' sub '.mat'], 'dat');

end




