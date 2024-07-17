function [] = VT_applyCSD(sub,sess,exp, epoch, par)

%Relies on CSD toolbox
%See https://psychophysiology.cpmc.columbia.edu/software/csdtoolbox/tutorial.html
thisPath = ['P' sub '\S' num2str(sess) '\EEG\'];
fullPath = ([exp.dataPath, thisPath]);

cd(fullPath);

if strcmp(epoch, 'SL')
    EEG = pop_loadset([fullPath '\acICASL_f01fc_P' sub '_S' num2str(sess) '.set']);
else
    EEG = pop_loadset([fullPath '\acICARL_f01fc_P' sub '_S' num2str(sess) '.set']);
end

if par.ds == 1 % Downsample data 
    keepEpoch = EEG.epoch;
    [EEG] = pop_resample(EEG,200); %Downsample data for speed 
    EEG.epoch = keepEpoch; 
    dsLabel = 'DS';
else
    dsLabel = '';
    
end
%Load CSD transform
load('CSD_coords_biosemi128.mat');

figure; subplot(1,2,1); plot(mean(EEG.data(19,:,:),3)'); title('No baseline')
EEG = eeg_checkset( EEG );

for t = 1:size(EEG.data,3)
    EEG.data(:,:,t) = EEG.data(:,:,t) - mean(EEG.data(:,round(0.4*EEG.srate):round(0.5*EEG.srate),t),2)%EEG.baselines(:,t)
end
EEG = eeg_checkset( EEG );
subplot(1,2,2); plot(mean(EEG.data(19,:,:),3)'); title('Pre-mask baseline')
sgtitle(['P' sub])
%% Apply CSD to epoched, BASELINED data, but only to EEG channels
X = [];
tic();
disp(['Running CSD for P' sub '...'])
for e = 1:length(EEG.epoch)
    X(:,:,e) = CSD(EEG.data(1:128,:,e),G,H);
end
toc();
%This takes about 3 minutes per participant!
EEG.data(1:128,:,:) = X;

EEG = eeg_checkset( EEG );
filename = ['csd_' dsLabel EEG.filename];
EEG = pop_saveset( EEG, filename, fullPath);

end