function [ch2interp, EEG] = VT_manualChanCheck(sub,sess,exp,EEG,thisPath)
exp.ImportantChans = [4 19 85 54 115]; % electrodes around CPz (where CPP is), Fz (CNV) and left and right motor cortex
% Also mark front-most channels:
exp.FrontChans = [71 72 80 81 93 94 103 70 73 79 82 92 95 102]; % at the very front edge SD can be high not necessarily because channels are bad but because person is blinking /moving eyes. We don't want to interpolate GOOD electrodes that later will HELP us pick up artifacts that we want to throw out.
exp.rerefchan = [23]; % pick a channel to re-reference the data to and get a second picture of variance across channels - important just because sometimes the reference used to load the EEG might itself have been bad during the recording...

% Load in the epoched data
% EEG = pop_loadset([exp.filepath '/P1_fatmc' (exp.name) '_P1.set']);

% Concatenate all epochs and get Standard deviation (SD) per channel
clear SD SDoz
conc = reshape(EEG.data(1:exp.nEEGchans,:,:),[exp.nEEGchans,size(EEG.data,2)*size(EEG.data,3)]); % concatenate all trials

conc2 = conc - repmat(conc(exp.rerefchan,:),[exp.nEEGchans,1]);
for q=1:exp.nEEGchans % include externals here as well, because it might be a good idea to check whether those are noisy too
    SD(q,1) = std(conc(q,:)); % measure S.D. of each channel
    SD2(q,1) = std(conc2(q,:));
end

% Are there any channels that stick out in terms of standard deviation?
% To check, plot the SD per channel:
figure; 
subplot(2,1,1); hold on; plot(SD(1:exp.nEEGchans,:)); ylim([0 200]) % we only plot the channels in the cap because external electrodes are often higher variance (e.g. you might be recording EMG) and annoyingly set the scale so you always have to zoom in, and the purpose here is to identify channels for interpolation which is ALWAYS only the 128 cap channels
%title(['subject ' num2str(s) ' ' allsubj{s}])
% mark the important channels specified above - this is just a visual aid to know which you should particularly consider
for e=1:length(exp.ImportantChans)
    plot([1 1]*exp.ImportantChans(e),[0 max(SD(exp.ImportantChans(e),:))],'b'); 
end
% mark the front edge channels as well, again a visual aid to know where blinks are likely to create higher variance (through no fault of the electrodes)
for e=1:length(exp.FrontChans)
    plot([1 1]*exp.FrontChans(e),[0 max(SD(exp.FrontChans(e),:))],'k'); % front edge channels in BLACK
end

%Identify candidate bad channels
candidates = find(SD(1:exp.nEEGchans,:) > 50 | (SD(1:exp.nEEGchans,:) < 1 & (SD2(1:exp.nEEGchans,:) < 1)));
exclCandidates = setdiff(candidates,exp.FrontChans);

for e=1:length(exclCandidates)
    plot([1 1]*exclCandidates(e),[0 max(SD(exclCandidates(e),:))],'r'); % candidate noisy channels in RED
end

subplot(2,1,2); hold on; plot(SD2(1:exp.nEEGchans,:)); ylim([0 200])
% mark the important channels specified above - this is just a visual aid to know which you should particularly consider
for e=1:length(exp.ImportantChans)
    plot([1 1]*exp.ImportantChans(e),[0 max(SD(exp.ImportantChans(e),:))],'b'); 
end
for e=1:length(exp.FrontChans)
    plot([1 1]*exp.FrontChans(e),[0 max(SD(exp.FrontChans(e),:))],'k'); % front edge channels in BLACK
end

ch2interp = exclCandidates; % makes an empty cell array for filling in the bad channels for this subject

disp(['Channels to interpolate: ' num2str(exclCandidates')]);
prompt = 'Confirm interpolation selection? y/n [y]: ';
str = input(prompt,'s');

if strcmp(str,'y')
    EEG = eeg_interp(EEG,ch2interp);
    EEG.interpolated = ch2interp;
    EEG = eeg_checkset( EEG );
 else
    disp('Interpolation aborted');
end


%% Another way to see if a channel that seems bad for certain blocks are actually bad channels, or is it a 
%% bad subject doing something whacky on some trials, is to get SD per channel per trial:
SDct = squeeze(std(EEG.data(1:exp.nEEGchans,:,:),[],2)); 
return; 
