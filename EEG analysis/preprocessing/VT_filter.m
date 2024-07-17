function [EEG] = VT_filter(sub,sess,exp,anPar) 

thisPath = ['P' sub '\S' num2str(sess) '\EEG\'];
fullPath = ([exp.dataPath, thisPath]);

cd(fullPath);

EEG = pop_loadset( [fullPath '\fc_P' sub '_S' num2str(sess) '.set'] );

% High pass filter at 0.1Hz
EEG  = pop_basicfilter( EEG,  1:exp.nEEGchans , 'Boundary', 'boundary', 'Cutoff', [anPar.EEG.highpassFilter], 'Design', 'butter', 'Filter', 'highpass', 'order', 4);
EEG = eeg_checkset(EEG);

% Notch filter at 50, 100 & 150 Hz
for thisFilter = anPar.EEG.notchFilters(1:end)
    EEG  = pop_basicfilter( EEG,  1:exp.nEEGchans , 'Boundary', 'boundary', 'Cutoff',  thisFilter, 'Design', 'notch', 'Filter', 'PMnotch', 'Order',  180 ); % GUI: 09-Mar-2022 09:48:27
    EEG = eeg_checkset(EEG);
end

filename = ['f01' EEG.filename]; %01 for 0.1 HP filter; if no label, hp filter = 0.01

%Save filtered data
EEG.filename = filename;
EEG = pop_saveset( EEG, filename, fullPath);

end