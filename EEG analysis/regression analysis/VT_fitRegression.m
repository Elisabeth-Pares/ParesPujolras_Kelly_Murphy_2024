function [dat] = VT_fitRegression_final(nsubj,exp,par) 

allsubj = exp.sub_id; % ORDER IS IMPORTANT FOR TRIAL-IDs!!!!
sub = allsubj{nsubj};

subj = ['P' allsubj{nsubj}];

% Labels derived from processing options
if par.ds == 1, par.dsLab = 'DS'; else par.dsLab = ''; end

smpwin = [0 1.4];  % window for sample-wise analyses
trlwin = [-0.1 5.4];  % window for full-trial-wise analyses
percReject = [];
% ==================================================================
% LOAD DATA
% ==================================================================

trl_data = [];allT_sess_data=[];sess_full = [];
figure
for sess = 1:exp.nSessions 
    trialNum = [];
    thisPath = ['P' sub '\S' num2str(sess) '\EEG\'];
    fullPath = ([exp.dataPath, thisPath]);
    cd(fullPath);
    
    % Load preprocessed & epoched EEG data
    EEG = pop_loadset( [fullPath '\' par.csdLab par.dsLab 'acICA' par.epoch '_f01fc_P' sub '_S' num2str(sess) '.set'] );
    
    % If TFA, load TFAdata too
    if strcmp(par.dataType, 'TFA')
        filename = ['tfa_', sub '_', num2str(sess) '.mat']
        load(filename);
    end
        
    % Baseline as required %% No need- data were baselined pre-CSD
    if strcmp(par.dataType, 'ERP')
        if strcmp(par.baseline, 'preTrial')
            %Get pre-epoch baselines
            baselineTime = round([(0.5-par.baselineTime)*EEG.srate]):round([0.5*EEG.srate]);
            EEG.baselines= squeeze(mean(EEG.data(:,baselineTime,:),2)); %Recalculate baselines
            %Specify baselining method
            if strcmp(par.baselineMethod, 'tbt_subtr') %Used in all analysis up to 20/09/2023 Power - baseline
                
                for t = 1:size(EEG.baselines,2)
                    EEG.data(:,:,t) = EEG.data(:,:,t) - EEG.baselines(:,t);
                end
            end
        end
    end
    
    plot(nanmean(EEG.data(19,:,:),3)')
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
        [~,sortI] = sort(trialNum);
        sess_data = sess_data(sortI,:,:);
    end
  
    % If dB transform
    if par.db == 1 & ~strcmp(par.dataType, 'ERP')
        %Add constant if there are negative values in sess_data to
        %apply db transform
        if ~isempty(find(any(sess_data<0)))
            %Include only motor channels
            toReject = ~ismember(1:128, exp.movCluster); 
            sess_data(:,toReject,:) = NaN; 
            constant = min(sess_data(:));
            sess_data = sess_data + abs(constant)+1;
        end
        assert(isempty(find(any(sess_data<0))))
        baselineTime = round([(0.5-par.baselineTime)*200]):round([0.5*200]);
        EEG.baselines= squeeze(nanmean(sess_data(:,:,baselineTime),3)); %Recalculate baselines
        sess_data = 10.*log10(sess_data./nanmean(EEG.baselines));  % convert to dB
    end
   
    % Baseline TF data
    if strcmp(par.dataType, 'TFA')
        if strcmp(par.baseline, 'preTrial')
            %Get pre-epoch baselines
            baselineTime = round([(0.5-par.baselineTime)*200]):round([0.5*200]);
            %Specify baselining method
            if strcmp(par.baselineMethod, 'allT_subtr') % As in https://www.frontiersin.org/articles/10.3389/fpsyg.2011.00236/full
                baselines(:,:) = (nanmean(nanmean(sess_data(:,:, baselineTime),3),1));
                sess_data(:,:,:) = (sess_data(:,:,:) - baselines(:,:));
            end
        end
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
LPR_full = dat.fitted_glaze.LPR_full;

%Check that length of behaviour & EEG are the same
% Participant 1 has some missing EEG trials
if strcmp(sub, '01'); idcs = [225:235];  psi_full(idcs,:) = []; llr_full(idcs,:) = []; surprise_full(idcs,:) = []; choices_full(idcs,:)=[]; 
    LPR_full(idcs,:) = [];
    wholeIdx = [321:335]; dat.fitted_glaze.LLR(wholeIdx,:) = [];
end

% Pick ONLY trials with full sequence -otherwise regression
% includes all trials, even those where sample times are out of trial
% bounds.
[ridx, cidx] = find(isnan(dat.fitted_glaze.LLR)); %Returns all columns with nans
trials2exclude = unique(ridx);

%Remove non-full trials;
trl_data(trials2exclude,:,:) = [];
sess_full(trials2exclude,:) = [];

assert(length(choices_full)==size(trl_data,1),'ERROR: Trial counts in eeg/behaviour are unequal')

sumData_full = squeeze(nanmean(trl_data,1));

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
prior_full(:,:) = psi_full(:,2:end); 

% ==================================================================
% High-pass filter signal, if required
% ==================================================================
if par.hpFilt == 1 & strcmp(par.dataType, 'ERP') 
    [bfilt,afilt] = butter(3, 1*2/200, 'high');   % hi-pass
    for t = 1:size(trl_data,1)
        for e = 1:size(trl_data,2)
            trl_data(t,e,:) = filtfilt(bfilt,afilt, trl_data(t,e,:));
        end
    end
end
 


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

% ==================================================================
% SEGMENT DATA AROUND INDIVIDUAL SAMPLES
% ==================================================================
fprintf('Concatenating trial- & sample-wise data segments...\n')
times = times(times>=trlwin(1) & times<=trlwin(2));

% R1: In participants 13-20, there was a small timing issue whereby SOA was
% 0.386 rather than 0.4s. Accounting for that here. 
if nsubj < 13 
    onsets = 0.4:0.4:0.4*10;  onsets=round(onsets,1);  % vector of all sample onset times relative to pre-mask; epoched around pre-mask on
else
    onsets = 0.386:0.386:0.386*10; 
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
deltaPsi = abs(diff(psi_full,1,2));

% Effective ev_signed
deltaPsi_signed = (diff(psi_full,1,2));

if par.saveData == 1
        dat.allTrial_data = squeeze(allT_sess_data(:,[1,3,19,21,23,87],:));
        dat.trl_data= squeeze(trl_data(:,[1,3,19,21,23,87],:));
        dat.smp_data= squeeze(smp_data(:,[1,3,19,21,23,87],:,:));
end

%R1: exclude samples with saccades? For analysis in Fig. S5
if par.excludeSacc == 1
    exTag = 'saccEx';
    %Load 
    filename = [subj '_processedEyeData.mat'];
    thisDir = [exp.dataPath subj '/S2/Eyetracking/'];
    load([thisDir '\' filename]);

    %Find trials to exclude entirely:
    % Have <11 saccades, and X & Y positions have less than 3*SD
    % variability, on average.
    goodTrials = eyedat.numSacc < 11; %Less than 11 saccades detected
    goodTrials = (eyedat.stdY < 3*nanmean(eyedat.stdY(goodTrials)) & eyedat.nanstdX <3*nanmean(eyedat.nanstdX(goodTrials)) & goodTrials);
    badTrials = goodTrials == 0;
    
    isfull = ~isnan(eyedat.sacMat(:,10))
    badTrials(~isfull) = []; 
    saccMat = eyedat.sacMat(isfull,:); %Saccade yes/no matrix

    % Participant 1 has some missing EEG trials
    if strcmp(sub, '01'); idcs = [225:235];  
        badTrials(idcs) = []; goodTrials(idcs) = [];
        saccMat(idcs,:) = []; 
    end
    %Check dimensions match 
    assert(length(badTrials) == size(smp_data,1));
    assert(size(saccMat,1) == size(smp_data,1));

    % Exclude data:
    nSamp = numel(deltaPsi); 
    nNan = sum(sum(isnan(deltaPsi))); 
    smp_data(badTrials,:,:) = NaN;
    deltaPsi(badTrials,:) = NaN; 
    llr_full(badTrials,:) = NaN; 
    session(badTrials,:) = NaN; 
    
    %Additionally, exclude single samples where saccades were detected:
    deltaPsi(logical(saccMat)) = NaN; 
    llr_full(logical(saccMat)) = NaN; 
    
    percReject = (sum(sum(isnan(deltaPsi)))-nNan)/nSamp*100; 
    
else
    exTag = '';
end
    
if par.fit == 1
    
    % ==================================================================
    % RUN SINGLE-TRIAL REGRESSIONS
    % ==================================================================
    % Sample-wise regressions, sample*electrode*time
    
    %% Dummy code to check what happens during zscoring -abs
    fprintf('Running regressions...\n Sample ')
    
    if strcmp(par.chans, 'parietal')
       if strcmp(par.dataType, 'TFA'); electrodes = 1:128;%exp.movCluster(:)';
       else
           electrodes = 1:128;
       end
       
    elseif strcmp(par.chans, 'lateralisation') & size(smp_data,2)>1 
        electrodes = 1; %lateralisation signal
        newSmpData =  nanmean(smp_data(:,exp.movCluster(2,:),:,:),2) - nanmean(smp_data(:,exp.movCluster(1,:),:,:),2); 
        clear smp_data;
        smp_data = newSmpData;
        
    end
    
    %LOGIT transform surprise to reduce skewness
    surprise_full = log(surprise_full./(1-surprise_full));  %log(pCP/(1-pCP))
   
    for s =1:9 % looping through samples
        fprintf('%d, ',s)
        for e = electrodes % looping through electrodes OR frequencies
            for t = 1:EEG.srate % looping through time-points (1 s after each sample)
                
                if strcmp(par.llr, 'rel')
                    if strcmp(par.regName, 'full')
                        regNames = {'prior', 'zLLR', 'zLLRxSurprise'};
                        regressors = nanzscore([prior_full(:,s) llr_full(:,s+1) nanzscore(llr_full(:,s+1)).*nanzscore(surprise_full(:,s+1)) sess_r]);
                        m = regstats(nanzscore(smp_data(:,e,t,s+1)), regressors ,'linear',{'beta','tstat','rsquare','adjrsquare','r'});  % signed prior, LLR, LLR*surprise & LLR*UNcertainty (-|psi|)
                    elseif strcmp(par.regName, 'psi_deltaPsi')
                        regNames = {'prior', 'deltaPsi'};
                        regressors = nanzscore([prior_full(:,s) deltaPsi_signed(:,s+1) sess_r]);
                        m = regstats(nanzscore(smp_data(:,e,t,s+1)),regressors,'linear',{'beta','tstat','rsquare', 'r', 'adjrsquare'});
                    elseif strcmp(par.regName, 'absPsi_deltaPsi_signedLat') %R1: reply to reviewers
                        regNames = {'prior', 'deltaPsi'};
                        regressors = nanzscore([abs(prior_full(:,s)) abs(deltaPsi_signed(:,s+1)) sess_r]);
                        m = regstats(nanzscore(smp_data(:,e,t,s+1)),regressors,'linear',{'beta','tstat','rsquare', 'r', 'adjrsquare'});
                    elseif strcmp(par.regName, 'LLR')
                        regNames = {'LLR'};
                        regressors = nanzscore([llr_full(:,s) sess_r]);
                        m = regstats(nanzscore(smp_data(:,e,t,s)),regressors,'linear',{'beta','tstat','rsquare', 'r', 'adjrsquare'});
                    elseif strcmp(par.regName, 'deltaPsi')
                        regNames = {'deltaPsi'};
                        regressors = nanzscore([deltaPsi_signed(:,s) sess_r]);
                        m = regstats(nanzscore(smp_data(:,e,t,s)),regressors,'linear',{'beta','tstat','rsquare', 'r', 'adjrsquare'});
                     elseif strcmp(par.regName, 'updatedPsi')
                        regNames = {'updatedPsi'};
                        m = regstats(nanzscore(smp_data(:,e,t,s)),nanzscore([prior_full(:,s) sess_r]),'linear',{'beta','tstat','rsquare','adjrsquare','r'});  % signed prior, LLR, LLR*surprise & LLR*UNcertainty (-|psi|)

                    end
                    
                elseif strcmp(par.llr, 'abs') % For CPP analysis
                    %  Regression is the same, but with abs rather than signed LLR values.
                    LLR_full = abs(llr_full);
                    if ~strcmp(par.regName, 'signedPsi'), prior_full = abs(prior_full); end
                    
                    if strcmp(par.regName, 'full')
                        regNames = {'prior', 'zLLR', 'zLLRxSurprise'};
                        m = regstats(nanzscore(smp_data(:,e,t,s+1)),nanzscore([prior_full(:,s) nanzscore(LLR_full(:,s+1)) nanzscore(LLR_full(:,s+1)).*nanzscore(surprise_full(:,s+1)) sess_r]),'linear',{'beta','tstat','rsquare','r'});  % signed prior, LLR, LLR*surprise & LLR*UNcertainty (-|psi|)
                    elseif strcmp(par.regName, 'deltaPsi')
                        regNames = {'deltaPsi'};
                        m = regstats(nanzscore(smp_data(:,e,t,s)),nanzscore([deltaPsi(:,s) sess_r]),'linear',{'beta','tstat','adjrsquare','r'});  % signed prior, LLR, LLR*surprise & LLR*UNcertainty (-|psi|)
                    elseif strcmp(par.regName, 'psi_deltaPsi') %R1: Psi_deltaPsi - reply to reviewers
                        regNames = {'psi', 'deltaPsi'};
                        m = regstats(nanzscore(smp_data(:,e,t,s+1)),nanzscore([prior_full(:,s) deltaPsi(:,s+1) sess_r]),'linear',{'beta','tstat','adjrsquare','r'});  % signed prior, LLR, LLR*surprise & LLR*UNcertainty (-|psi|)
                    elseif strcmp(par.regName, 'signed_deltaPsi') %R1: signed_deltaPsi - reply to reviewers
                        regNames = {'deltaPsi'};
                        m = regstats(nanzscore(smp_data(:,e,t,s)),nanzscore([deltaPsi_signed(:,s) sess_r]),'linear',{'beta','tstat','adjrsquare','r'});  % signed prior, LLR, LLR*surprise & LLR*UNcertainty (-|psi|)
                    elseif strcmp(par.regName, 'updatedPsi')
                        regNames = {'updatedPsi'};
                        m = regstats(nanzscore(smp_data(:,e,t,s)),nanzscore([prior_full(:,s) sess_r]),'linear',{'beta','tstat','rsquare','adjrsquare','r'});  % signed prior, LLR, LLR*surprise & LLR*UNcertainty (-|psi|)
                    elseif strcmp(par.regName, 'prior')
                        regNames = {'prior'};
                        m = regstats(nanzscore(smp_data(:,e,t,s)),nanzscore([prior_full(:,s), sess_r]),'linear',{'beta','tstat','rsquare', 'adjrsquare'});  % signed prior, LLR, LLR*surprise & LLR*UNcertainty (-|psi|)
                    elseif strcmp(par.regName, 'LLR')
                        m = regstats(nanzscore(smp_data(:,e,t,s)),nanzscore([LLR_full(:,s), sess_r]),'linear',{'beta','tstat','rsquare', 'adjrsquare'});  % signed prior, LLR, LLR*surprise & LLR*UNcertainty (-|psi|)
                        regNames = {'LLR'};
                    elseif strcmp(par.regName, 'surprise')
                        m = regstats(nanzscore(smp_data(:,e,t,s)),nanzscore([surprise_full(:,s) sess_r]),'linear',{'beta','tstat','rsquare', 'adjrsquare'});  % signed prior, LLR, LLR*surprise & LLR*UNcertainty (-|psi|)
                        regNames = {'surprise'};
                    elseif strcmp(par.regName, 'intercept')
                        m = regstats(nanzscore(smp_data(:,e,t,s)),nanzscore([sess_r]),'linear',{'beta','tstat','rsquare', 'adjrsquare'});  % signed prior, LLR, LLR*surprise & LLR*UNcertainty (-|psi|)
                        regNames = {'intercept'};
                    end
                end
                
                % Save beta coefficients and residuals
                for r = 1:length(regNames)
                    thisReg = regNames{r};
                    dat.(thisReg).beta(e,t,s) = m.beta(r+1); %Beta coefficients

                    if r == 1 
                        if strcmp(par.regName, 'deltaPsi') & strcmp(par.dataType, 'ERP') & par.hpFilt == 1
                            dat.(thisReg).r(e,t,s,:) = m.r; %Residuals
                        end
                        if par.save_aR2 == 1
                            dat.(thisReg).adjrsquare(e,t,s,:) = m.adjrsquare; %Adjusted Rsquared
                        end
                    end
                end
                
            end
        end
    end

dat.percReject = percReject; 
dat.modelType = 'fitted_glaze';
fprintf('Done.\n')

thisPath = ['P' sub];
fullPath = ([exp.dataPath, thisPath]);

save([fullPath, '\VT_regression_' par.regName '_' par.chans '_' par.dataType, '_rsq_HP_' num2str(par.hpFilt) exTag '_P' sub, '_' num2str(sess)], 'dat','-v7.3')

end
end




