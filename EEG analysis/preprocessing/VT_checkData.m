function [EEG]= VT_checkData(EEG,sub,sess, exp)

%Data should have same number of block starts & ends
%Count block start & end triggers
if isa([EEG.event(3).type], 'char')
    blockStarts = find(strcmp({EEG.event.type}, num2str(exp.EEG.triggers.block_start)));
    blockEnds = find(strcmp({EEG.event.type}, num2str(exp.EEG.triggers.block_end)));
else
    blockStarts = find([EEG.event.type], (exp.EEG.triggers.block_start));
    blockEnds = find([EEG.event.type], (exp.EEG.triggers.block_end));
end

%Remove duplicate block ends
diffEnds = [1, diff(blockEnds)];
EEG.event(blockEnds(diffEnds == 1)) = []; %Remove from data
blockEnds(diffEnds == 1) = [];

if isa([EEG.event(3).type], 'char')
    blockStarts = find(strcmp({EEG.event.type}, num2str(exp.EEG.triggers.block_start)));
    blockEnds = find(strcmp({EEG.event.type}, num2str(exp.EEG.triggers.block_end)));
else
    blockStarts = find([EEG.event.type], (exp.EEG.triggers.block_start));
    blockEnds = find([EEG.event.type], (exp.EEG.triggers.block_end));
end

%Exceptions:
% P01, S1: two end of block triggers were sent; second half of eeg data
% started late; block 2 needs to be eliminated because pupil was lost.
% P03, S1: two extra start of block triggers & subsequent trials after
% block 4.

%% ===================================================================================
% General participant data annotation
%% ===================================================================================
% Add block#
EEG = VT_annotateBlocks(exp,EEG);

% Add trial#
EEG = VT_annotateTrials(exp,EEG);

%% ===================================================================================
% Special participant data cleaning - specific errors that need to be
% handled individually.
%% ===================================================================================

if strcmp(sub, '01') && sess == 1 %Peter's first session had two end of block triggers
    %Remove block 2
    block2remove = 2;
    for b = flip(block2remove(1:end))
        thisStart = blockStarts(b);
        nextStart = blockStarts(b+1);
        EEG.event(thisStart:nextStart-1) = [] %Delete that chunk
    end
    
    % Add block#
    EEG = VT_annotateBlocks(exp,EEG);
    
    % Add trial#
    allEEG = []; thisEEGevents = [];
    for b = unique([EEG.event.block])
        subset = [EEG.event([EEG.event.block] == b)];
        thisEEGevents.event = subset;
        fixedEEG = VT_annotateTrials(exp,thisEEGevents);
        if isempty(allEEG)
            allEEG = fixedEEG.event;
        else
            allEEG = [allEEG, fixedEEG.event];
        end
    end
    
    EEG.event = allEEG; clear allEEG;

    
    %Correct trial ID numbers adding 14 to all trials starting in block 5.
    %This is because total number of trials in the session should be 8*80 =
    %640, and there are only 626. Hence, 14 are missing.
    earlyTrials = [EEG.event([EEG.event.block] < 5).trial];
    lateTrials = [EEG.event([EEG.event.block] >5).trial];
    newTrials = [EEG.event([EEG.event.block] == 5).trial] + 14;
    correctedTrials = [earlyTrials, newTrials,lateTrials];
    for t = 1:length([EEG.event])
        [EEG.event(t).trial] = correctedTrials(t);
    end
    
    %Correct ntrial numbers
    [EEG.event.ntrial] = deal(0);
    trials = [EEG.event.trial] + (([EEG.event.block]-1)*80);
    for i = 1:length(EEG.event)
        [EEG.event(i).ntrial] = trials(i);
    end
    
elseif strcmp(sub, '02') && sess == 1 %P2 second session had two short extra blocksfirst session had two end of block triggers
  [EEG.event.block] = deal(0);
    for t = 1:length(EEG.event)
        if [EEG.event(t).trial] <=80; [EEG.event(t).block] = 1;
        elseif [EEG.event(t).trial] >80 & [EEG.event(t).trial] <=160; [EEG.event(t).block] = 2;
        elseif [EEG.event(t).trial] >80*2 & [EEG.event(t).trial] <=80*3; [EEG.event(t).block] = 3;
        elseif [EEG.event(t).trial] >80*3 & [EEG.event(t).trial] <=80*4; [EEG.event(t).block] = 4;
        elseif [EEG.event(t).trial] >80*4 & [EEG.event(t).trial] <=80*5; [EEG.event(t).block] = 5;
        elseif [EEG.event(t).trial] >80*5 & [EEG.event(t).trial] <=80*6; [EEG.event(t).block] = 6;
        elseif [EEG.event(t).trial] >80*6 & [EEG.event(t).trial] <=80*7; [EEG.event(t).block] = 7;
        elseif [EEG.event(t).trial] >80*7 & [EEG.event(t).trial] <=80*8; [EEG.event(t).block] = 8;
        elseif [EEG.event(t).trial] >80*8 & [EEG.event(t).trial] <=80*9; [EEG.event(t).block] = 9;
        end
    end

    %Correct ntrial numbers
    [EEG.event.ntrial] = deal(0);
    
    for i = 1:length(EEG.event)
        [EEG.event(i).ntrial] = [EEG.event(i).trial];
        [EEG.event(i).trial] = [EEG.event(i).trial] - (([EEG.event(i).block]-1)*80);
    end
    
elseif strcmp(sub, '03') && sess == 1 %P2 second session had two short extra blocksfirst session had two end of block triggers
    %Remove block 2
    block2remove = [4,6,7];
    for b = flip(block2remove(1:end))
        EEG.event([EEG.event.block] ==b) = [] %Delete that chunk
    end
    
    % Add trial#
    EEG = VT_annotateTrials(exp,EEG);

    
    [EEG.event.block] = deal(0);
    for t = 1:length(EEG.event)
        if [EEG.event(t).trial] <=80; [EEG.event(t).block] = 1;
        elseif [EEG.event(t).trial] >80 & [EEG.event(t).trial] <=160; [EEG.event(t).block] = 2;
        elseif [EEG.event(t).trial] >80*2 & [EEG.event(t).trial] <=80*3; [EEG.event(t).block] = 3;
        elseif [EEG.event(t).trial] >80*3 & [EEG.event(t).trial] <=80*4; [EEG.event(t).block] = 4;
        elseif [EEG.event(t).trial] >80*4 & [EEG.event(t).trial] <=80*5; [EEG.event(t).block] = 5;
        elseif [EEG.event(t).trial] >80*5 & [EEG.event(t).trial] <=80*6; [EEG.event(t).block] = 6;
        elseif [EEG.event(t).trial] >80*6 & [EEG.event(t).trial] <=80*7; [EEG.event(t).block] = 7;
        elseif [EEG.event(t).trial] >80*7 & [EEG.event(t).trial] <=80*8; [EEG.event(t).block] = 8;
        elseif [EEG.event(t).trial] >80*8 & [EEG.event(t).trial] <=80*9; [EEG.event(t).block] = 9;
        end
    end

    %Correct ntrial numbers
    [EEG.event.ntrial] = deal(0);
    trials = [EEG.event.trial] - (([EEG.event.block]-1)*80);
    for i = 1:length(EEG.event)
        [EEG.event(i).ntrial] = [EEG.event(i).trial];
        [EEG.event(i).trial] = trials(i);
    end
    
end

end
