function [] = VT_epochData(sub,sess,exp, anPar)

thisPath = ['P' sub '\S' num2str(sess) '\EEG\'];
fullPath = ([exp.dataPath, thisPath]);

cd(fullPath);

clear chanlocs;

for e = [2]
    thisEpoch = anPar.EEG.epochs{e}
    
    % Load preprocessed data and epoch around trial start (SL), each patch (FL)
    % and the response (RL)
    EEG = pop_loadset( [fullPath '\f01fc_P' sub '_S' num2str(sess) '.set'] );
    
    % Add channel locations
    load chanlocsBioSemi_128_noext;
    elab = strings(exp.nEEGchans,1);
    for indChan = 1:exp.nEEGchans
        elab(indChan) = chanlocs(indChan).labels;
    end
    elab(end)='FCz'; clear ind*
    
    EEG.chanlocs = chanlocs;
    EEG.data(129:137,:) = []; %Remove extra channels

    if e ~= 3, thisTrigger = anPar.EEG.epochTriggs(e); else thisTrigger = [anPar.EEG.epochTriggs{e} anPar.EEG.epochTriggs(e+1)]; end %Take both left/right response triggers

    EEG = pop_epoch( EEG, thisTrigger{1}, anPar.EEG.epochLims(e,:), 'epochinfo', 'yes');
    EEG = eeg_checkset( EEG );
       
    %Get pre-epoch baselines 
    baselineTime = [0.1]*EEG.srate;
    EEG.baselines= squeeze(mean(EEG.data(:,1:baselineTime,:),2)); 
    
    %Interpolate bad channels
    if e == 1 %Use bad channels as detected in SL data to reject from all other epoching methods 
        [ch2interpolate, EEG] = VT_manualChanCheck(sub,sess,exp,EEG,fullPath);
        SLbaselines = EEG.baselines;
    else
        EEG1 = pop_loadset( [fullPath '\ICASL_f01fc_P' sub '_S' num2str(sess) '.set'] );
        EEG = eeg_interp(EEG,EEG1.interpolated); %Interpolate same channels as in RL 
    end
    
    EEG = eeg_checkset( EEG );
    
    %Re-reference epoched data
    EEG = pop_reref( EEG, [1:128] ,'exclude',[129:136] ,'keepref','on');
    EEG = eeg_checkset( EEG );
    
    if e == 2 %Upload ICA from SL data
        EEG.interpolated = EEG1.interpolated;
        EEG = eeg_checkset( EEG );
        EEG.icaact = EEG1.icaact;
        EEG.icaweights = EEG1.icaweights;
        EEG.icawinv = EEG1.icawinv;
        EEG.icachansind = EEG1.icachansind;
        EEG.icasphere = EEG1.icasphere;
    end
    %Save data
    EEG = pop_saveset( EEG,[thisEpoch '_' EEG.filename], fullPath);
    
   clear EEG; 
end

end

