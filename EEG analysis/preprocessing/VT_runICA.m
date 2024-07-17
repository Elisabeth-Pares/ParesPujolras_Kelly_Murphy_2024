function [EEG] = VT_runICA(sub,sess,exp, anPar)

thisPath = ['P' sub '\S' num2str(sess) '\EEG\'];
fullPath = ([exp.dataPath, thisPath]);

cd(fullPath);

clear chanlocs;

for e = 1
    thisEpoch = anPar.EEG.epochs{e}
    
    % Load preprocessed data and epoch around trial start (SL), each patch (FL)
    % and the response (RL)
    EEG = pop_loadset( [fullPath '\' thisEpoch '_f01fc_P' sub '_S' num2str(sess) '.set'] );
    
    [EEG1] = pop_resample(EEG,200); %Try downsampling just for ICA
    
    tic()
    EEG1 = pop_runica(EEG1, 'icatype', 'runica', 'extended',0,'interrupt','on', 'chanind', [1:128]) % exp.icaChans);
    toc()
    
    EEG.icawinv = EEG1.icawinv;
    EEG.icasphere = EEG1.icasphere;
    EEG.icaweights = EEG1.icaweights;
    EEG.icachansind = EEG1.icachansind;
    
    % Upload ICA components to full-sample dataset
    EEG=  pop_saveset(EEG, ['ICA' EEG.filename], fullPath); %  for corrected
    
end
