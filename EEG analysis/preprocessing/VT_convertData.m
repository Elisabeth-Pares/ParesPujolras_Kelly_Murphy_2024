function [] = VT_convertData(sub,sess, exp)

allEEG = [];

thisPath = ['P' sub '\S' num2str(sess) '\EEG\'];
fullPath = ([exp.dataPath, thisPath]);

files = dir(fullPath);
files = files(contains({files.name}, '.bdf'))
files = {files.name}

%Load and convert both halves    

trialCounts = [];
full_EEG = [];
for f = 1:length(files)
    %Convert data from .bdf (BioSemi) format to .set (EEGlab) format
    thisFile = files{f};
    cd(fullPath)
    if exist(thisFile, 'file')
    
        EEG = pop_fileio([fullPath thisFile], 'dataformat', 'auto');
        EEG = eeg_checkset(EEG);
        
        %Add block code 
        thisBlock = str2num(thisFile(end-4));
        if thisBlock == 0; thisBlock = 10;end            
        [EEG.event.block] = deal(thisBlock);
        
        %Add trial number (in block)       
        EEG = VT_annotateTrials(exp,EEG)
        
        %Add total trial number 
        [EEG.event.ntrial] = deal(0);
        nTrials = [EEG.event.trial] + ([EEG.event.block] -1) * 80;
        for i = 1:length(EEG.event)
            [EEG.event(i).ntrial] = nTrials(i);
        end
        
        trialCounts = [trialCounts; sess, thisBlock, max([EEG.event.trial])];
        if max([EEG.event.trial]) ~= 80; warning('Missing trials!'); end

        filename = ['c' thisFile(1:end-4)];
        EEG = pop_saveset(EEG, filename, fullPath);
        
        if isempty(full_EEG)
            full_EEG = EEG;
        else
            full_EEG = pop_mergeset(full_EEG, EEG);
        end
    end
end

% Deal with participants with special problems 
if (strcmp(sub, '01') || strcmp(sub, '02') || strcmp(sub, '03')) & sess == 1
    full_EEG = VT_checkData(full_EEG,sub,sess, exp)
end

if ~isempty(full_EEG)
    %Concatenate all blocks for each session; bad channel cleaning to be done
    %once per session.
    full_EEG.trialCounts = trialCounts; 
    filename = ['fc_P' sub '_S' num2str(sess)];
    EEG = pop_saveset(full_EEG, filename, fullPath);
end
end