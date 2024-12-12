function [gaDat] = VT_psiQuartiles_ERP(nsubj,exp,par, gaDat) 
allsubj = exp.sub_id; % ORDER IS IMPORTANT FOR TRIAL-IDs!!!!
sub = allsubj{nsubj};

subj = ['P' allsubj{nsubj}];%find(strcmp(allsubj,subj));

% Paths dervied from processing options

smpwin = [0 1.4];  % window for sample-wise analyses
trlwin = [-0.1 5.4];  % window for full-trial-wise analysesd
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
for sess = 1:exp.nSessions %(strcmp(allsubj,sub))
    trialNum = [];
    thisPath = ['P' sub '\S' num2str(sess) '\EEG\'];
    fullPath = ([exp.dataPath, thisPath]);
    semifullPath = ([exp.dataPath, 'P' sub '\']);
    cd(fullPath);
    
    % Load preprocessed & epoched EEG data
    EEG = pop_loadset( [fullPath '\csd_DSacICASL_f01fc_P' sub '_S' num2str(sess) '.set'] );
    
   
    xline([0.4*EEG.srate]); xline([0.5*EEG.srate]); sgtitle(['P' sub, ', S' num2str(sess)])
    hold on;
    
    %Permute order so it fits Peter's default
    if strcmp(par.dataType, 'ERP')
        % EPP adapted: construct data matrix (trials*electrodes*times)
        sess_data = permute(EEG.data,[3,1,2]);
    elseif strcmp(par.dataType, 'TFA')
        sess_data = freqData; %dimensions are trials*freqs*times;
    end
    
    % CHECK TRIAL ORDER
    for t = 1:size([EEG.epoch],2)
        trialNum(t)= ([EEG.epoch(t).eventntrial{1}]);
    end
    
    if ~issorted(trialNum)  % if trials aren't in ascending order, sort - this happens in participant 01, session 02
        [mne_ids,sortI] = sort(trialNum);
        sess_data = sess_data(sortI,:,:);
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


% EPP edited: pick ONLY trials with full sequence -otherwise regression
% includes all trials, even those where sample times are out of trial
% bounds.
[ridx, cidx] = find(isnan(dat.fitted_glaze.LLR)); %Returns all columns without nans
trials2exclude = unique(ridx);

%Remove non-full trials;
trl_data(trials2exclude,:,:) = [];
sess_full(trials2exclude,:) = [];

assert(length(choices_full)==size(trl_data,1),'ERROR: Trial counts in eeg/behaviour are unequal')

% ==================================================================
% MAKE CATEGORICAL SESSION REGRESSORS
% ==================================================================
sessions = unique(sess_full);
sess_r = zeros(length(sess_full),length(sessions)-1);
for s = 1:length(sessions)-1
    sess_r(sess_full==sessions(s+1),s) = ones(length(find(sess_full==sessions(s+1))),1);
end

% ==================================================================
% ENSURE USE OF DESIRED FORM OF PRIOR BELIEF
% ==================================================================
% prior_full = LPR_full;
prior_full(:,:) = psi_full(:,2:end);  % N.B. this way, final term in LPR_full will always be final, untransformed belief (i.e. the quantity that the model uses to make a decision)

% ==================================================================
% Low-pass filter signal, if required
% ==================================================================
if par.lowFilt == 1 & strcmp(par.dataType, 'ERP');
    [bfilt,afilt] = butter(3, 6*2/200, 'low');   % lo-pass
    for t = 1:size(trl_data,1)
        for e = 1:size(trl_data,2)
            trl_data(t,e,:) = filtfilt(bfilt,afilt, trl_data(t,e,:));
        end
    end
end

% figure;
% plot(squeeze(nanmean(trl_data(:,19,:)))); hold on; plot(squeeze(nanmean(new_trl_data(:,19,:))))

% ==================================================================
% SEGMENT DATA AROUND INDIVIDUAL SAMPLES
% ==================================================================
fprintf('Concatenating trial- & sample-wise data segments...\n')
times = times(times>=trlwin(1) & times<=trlwin(2));
if nsubj < 13
    onsets = 0.4:0.4:0.4*10;  onsets=round(onsets,1);  % vector of all sample onset times relative to pre-mask; epoched around pre-mask on
else
    onsets = 0.386:0.386:0.386*10; % onsets=round(onsets,1);  % vector of all sample onset times relative to pre-mask; epoched around pre-mask on
end
smptimes = times(times>=smpwin(1) & times<=smpwin(2));  % getting vector of sample times relative to dot onset

for s = 1:length(onsets)
    stsmp = find(times>=onsets(s),1,'first');
    smp_data(:,:,:,s) = trl_data(:,:,stsmp:stsmp+length(smptimes)-2); %Smp_data dimensions are [trials x channels x timepoints x samples]
    smp_baselines(:,:,:,s) = nanmean(trl_data(:,:,stsmp-20:stsmp),3);
end

if strcmp(par.baseline, 'preSample')
    %Get pre-epoch baselines
    b_smp_data = smp_data;
    for ss = 1:10
        b_smp_data(:,:,:,ss) =  smp_data(:,:,:,ss)- smp_baselines(:,:,:,ss);
    end
    smp_data = [];
    smp_data = b_smp_data;
end

% Effective prior
effPrior = abs(diff(psi_full,1,2));

% Effective ev_signed
effPrior_signed = (diff(psi_full,1,2));


%% Find extreme |LLRs| & median split by deltaPsi
% to illustrate that given = LLR, high deltaPsi --> higher evoked ERP

rllr_full = llr_full;
rllr_full = reshape(rllr_full,[size(rllr_full,1)*size(rllr_full,2),1]); 
reffPrior = reshape(effPrior,[size(effPrior,1)*size(effPrior,2),1]); 

okIdx = abs(rllr_full) >= par.maxLLR;%== maxLLR;
nanIdx = abs(rllr_full) < par.maxLLR;%~= maxLLR;

maxLLR_effPrior = reffPrior; 
maxLLR_effPrior(nanIdx) = NaN;

[groupIdx] = quantileranks(abs(reffPrior(okIdx)), par.nGroups);

reshapedData = squeeze(nanmean(smp_data(:,exp.centroParietal,:,:),2));
reshapedData = permute(reshapedData, [1,3,2]);
reshapedData = reshape(reshapedData,[],size(reshapedData,3),1);

s_reffPrior = reffPrior(okIdx);
s_rllr_full = rllr_full(okIdx);
s_reshapedData = reshapedData(okIdx,:);
for i = 1:par.nGroups
    idx{i} = groupIdx == i; 
    groupLab = ['group' num2str(i)];
    gaDat.(groupLab).data(nsubj,:) = nanmean(s_reshapedData(idx{i},1:200));
    gaDat.(groupLab).mean_deltaPsi(nsubj,:) = nanmean((s_reffPrior(idx{i}))); 
    gaDat.(groupLab).mean_LLR(nsubj,:)  = nanmean(abs(s_rllr_full(idx{i})));
end


end





