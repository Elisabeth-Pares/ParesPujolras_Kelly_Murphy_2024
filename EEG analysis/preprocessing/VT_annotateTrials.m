function EEG = VT_annotateTrials(exp,EEG)

disp('Annotating trials...')

if isa([EEG.event(3).type], 'char')
   %Add trial numbers to EEG data structure
    trialEnds = find(strcmp({EEG.event.type}, num2str(exp.EEG.triggers.trial_end)));
    
    %For each trial
    for t = 1:length(trialEnds(1:end))
        thisEnd = trialEnds(t);
        prevEnd = find(strcmp({EEG.event(1:thisEnd-1).type}, num2str(exp.EEG.triggers.trial_end)), 1, 'last');
        prevStart = find(strcmp({EEG.event(1:thisEnd-1).type}, num2str(exp.EEG.triggers.premask_on)), 1, 'last');
        
        if prevEnd > prevStart %No start trigger for that trial
            [EEG.event(prevEnd+1:thisEnd).trial] = deal(t);
        else
            [EEG.event(prevStart:thisEnd).trial] = deal(t);
        end
    end
    
elseif isa([EEG.event(3).type], 'double')
   
     %Add trial numbers to EEG data structure
    trialEnds = find([EEG.event.type]== exp.EEG.triggers.trial_end);
    
    [EEG.event.trial] = deal(0); %Set default to 0;
    
    %For each trial
    for t = 1:length(trialEnds(1:end))
        thisEnd = trialEnds(t);
        prevEnd = find([EEG.event(1:thisEnd-1).type] == exp.EEG.triggers.trial_end, 1, 'last');
        prevStart = find([EEG.event(1:thisEnd-1).type]== exp.EEG.triggers.premask_on, 1, 'last');
        
        if prevEnd > prevStart %No start trigger for that trial
            [EEG.event(prevEnd+1:thisEnd).trial] = deal(t);
        else
            [EEG.event(prevStart:thisEnd).trial] = deal(t);
        end
    end
end