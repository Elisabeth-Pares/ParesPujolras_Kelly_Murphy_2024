function  [EEG] = VT_annotateBlocks(exp,EEG)

disp('Annotating blocks...')

blockStarts = find(strcmp({EEG.event.type}, num2str(exp.EEG.triggers.block_start)));
blockEnds = find(strcmp({EEG.event.type}, num2str(exp.EEG.triggers.block_end)));

if ~isfield(EEG.event, 'block'), [EEG.event(:).block] = deal(0); end


if length(blockStarts) ~=length(blockEnds)
    disp('Checking mismatches...')
    
    % In some participants, blocks were started and stopped due to e.g. calibration errors, so some data will
    % have block starts but not block ends.
    % Delete all data between a block start trigger that has no end.
    if length(blockStarts) > length(blockEnds)
        
        for b = 1:length(blockStarts)
            thisStart = blockStarts(b);
            nextStart = find(strcmp({EEG.event(thisStart:end).type}, num2str(exp.EEG.triggers.block_start)));
            nextEnd = find(strcmp({EEG.event(thisStart:end).type}, num2str(exp.EEG.triggers.block_end)));
            
            if nextStart < nextEnd %No end trigger for that block
                EEG.event(thisStart:nextStart-1) = [] %Delete that chunk
                b = b-1; %Repeat block processing;
            else
                %Annotate block
                EEG.event(thisStart:nextEnd).block = b;
            end
            
        end
        % In some participants, the recording was started late
    elseif length(blockStarts)<length(blockEnds)
        for b = 1:length(blockEnds(1:end))
            
            thisEnd = blockEnds(b);
            prevEnd = find(strcmp({EEG.event(1:thisEnd-1).type}, num2str(exp.EEG.triggers.block_end)), 1, 'last');
            prevStart = find(strcmp({EEG.event(1:thisEnd-1).type}, num2str(exp.EEG.triggers.block_start)), 1, 'last');
            
            if prevEnd > prevStart %No start trigger for that block
                %Find last non-zero value
                %Annotate block
                [EEG.event(prevEnd+1:thisEnd).block] = deal(b);
            else
                %Annotate block
                [EEG.event(prevStart:thisEnd).block] = deal(b);
            end
            
        end
        
    end
    
end