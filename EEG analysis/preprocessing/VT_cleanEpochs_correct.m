function [EEG] = VT_cleanEpochs_correct(sub,sess,exp, anPar)

thisPath = ['P' sub '\S' num2str(sess) '\EEG\'];
fullPath = ([exp.dataPath, thisPath]);

cd(fullPath);

clear chanlocs;

for e = 1
    thisEpoch = anPar.EEG.epochs{e}
    
    % Load preprocessed data and epoch around trial start (SL), each patch (FL)
    % and the response (RL)
    EEG = pop_loadset( [fullPath '\cICA' thisEpoch '_f01fc_P' sub '_S' num2str(sess) '.set'] );
   
    %Change search epoch according to epoch times
    if e == 1
        EEG1 = pop_eegthresh(EEG, 1, 1:EEG.nbchan,-150 , 150 ,-0.1 , 5, 1, 0); [EEG1,comrej, badlist]= pop_TBT(EEG1,EEG1.reject.rejthreshE,20,0.3,1,EEG.chanlocs);
    elseif e == 2
        EEG1 = pop_eegthresh(EEG, 1, 1:EEG.nbchan,-150 , 150 ,-2 , 0.5, 1, 0); [EEG1,comrej, badlist]= pop_TBT(EEG1,EEG1.reject.rejthreshE,20,0.3,1, EEG.chanlocs);
    end
    % %For P1, for example, see correction on  Epoch 14 for ch 95,101
    % %figure; plot(EEG.data(50,:,14)); hold on ; plot(EEG1.data(50,:,14))
    if ~isempty(badlist.bTrial_num)
        EEG.data(:,:,badlist.bTrial_num) = NaN; %exclude rejected trials;
        %     EEG.epoch(badlist.bTrial_num) = [];
        %     EEG.event = EEG1.event; EEG.urevent = EEG1.urevent;
        EEG.trials = size(EEG.data,3);
    else
        disp('Debugging')
    end
    
    clear EEG.data; clear EEG.epochs; 
    %Rewrite EEG channels on original datafile; keep externals as they were.
    EEG.data = EEG1.data;
    EEG.trials = size(EEG.data,3); 
    EEG.epoch = EEG1.epoch; 
    EEG.badList_TBT = badlist;
    
    EEG = eeg_checkset( EEG );
    filename = ['a' EEG.filename];
    EEG = pop_saveset( EEG, filename, fullPath);
    
end
end


