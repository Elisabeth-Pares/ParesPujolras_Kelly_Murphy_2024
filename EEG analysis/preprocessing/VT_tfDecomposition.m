function [] = VT_tfDecomposition(sub,sess,exp,anPar, par)    


% Load epoched data to transform 
thisPath = ['P' sub '\S' num2str(sess) '\EEG\'];
fullPath = ([exp.dataPath, thisPath]);
behPath = ['D:\Elisabeth EEG Data\VolatilityTask\Modelling\PreprocessedStim'];

cd(fullPath);

% Load preprocessed & epoched EEG data
EEG = pop_loadset( [fullPath '\' par.csdLab 'DSacICA' par.epoch '_f01fc_P' sub '_S' num2str(sess) '.set'] );

% Transform EEGlab data format to fieldtrip for analysis 
data = eeglab2fieldtrip(EEG,'preprocessing','none')  %Previously: raw

%Sample info
for e = 1:size(EEG.data,3)
    eventLat = [EEG.event.urevent] == EEG.epoch(e).eventurevent{1};
    disp(['Epoch ' num2str(e) 'eventLat: ' num2str(find(eventLat))]);
    if isempty(find(eventLat))
        disp('Oops!')
        data.sampleinfo(e,1:2) = [NaN, NaN];
    else
    trialOnset = EEG.event(eventLat).latency; 
    trialEnd = trialOnset + 5.4*EEG.srate; 
    data.sampleinfo(e,1:2) = [trialOnset - (0.5*EEG.srate), trialEnd];
    end
    clear trialOnset; clear trialEnd;
end

% Take subset of trials for debugging
if par.debug == 1
    data.trial(11:end) = [];
    data.time(11:end) = [];
    data.trialinfo(11:end,:) = [];
end

% Hanning decomposition 
% save those cfgs for later plotting
freq            = ft_freqanalysis(anPar.TFA, data);

if strcmp(par.epoch, 'SL')
    filename = ['tfa_', sub '_', num2str(sess) '.mat']
%     %Nanify bad trials
    freq.powspctrm(EEG.badList_TBT.bTrial_num,:,:,:) = NaN; 
        
    freqData = squeeze(mean(freq.powspctrm(:,:,13:30,:),3)); 
    save(filename, 'freqData', '-v7.3'); % single trial data
    
elseif strcmp(par.epoch, 'RL')
    %Save left v right motor responses summary only; beta band
    load([behPath '\' sub '-' num2str(sess) '-StimData.mat']) %Load stimuli
    load([behPath '\' sub '-' num2str(sess) '_vars.mat']) %Load stimuli
    fullTrials = ~isnan(mdlvars.distseq(:,10))
    for i = 1:size(EEG.data,3)
        idx(i) = [EEG.epoch(i).eventntrial{1}]
    end
    
%     Nanify bad trials
    freq.powspctrm(EEG.badList_TBT.bTrial_num,:,:,:) = NaN;
   
    %Plot Beta L-R and sort by response
    leftChoice = squeeze(nanmean(nanmean(freq.powspctrm(fullTrials(idx) & choices_full(idx) == 0,:,13:30,:),3)));
    rightChoice = squeeze(nanmean(nanmean(freq.powspctrm(fullTrials(idx) & choices_full(idx) == 1,:,13:30,:),3)));

    diffData = leftChoice-rightChoice;
    baseline = nanmean(diffData(:,1:50),2);
    bdiffData = diffData - baseline; 
    %Topography of -200 to -100 ms pre action
    t_diffData  = squeeze(nanmean(diffData(:,round(1.7*200): round(2*200)),2));
    figure; topoplot(t_diffData, EEG.chanlocs,'electrodes', 'on'); title('Left-Right choices')
    caxis([-5 5]); 

    filename = ['MotorLocaliser_RL_' sub '_' num2str(sess) '.mat']; 
    save(filename, 'leftChoice','rightChoice'); % single trial data
end

end