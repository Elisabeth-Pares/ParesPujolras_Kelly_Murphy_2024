function [] = VT_Figure2A(exp)
%% Load for all participants and save in the same file 
id = 0; clear leftChoice; clear rightChoice; 

for sub = 1:20
    id = id +1;
    subj = exp.sub_id{sub}
    %Load full trials
    load(['D:\Elisabeth EEG Data\VolatilityTask\Data\P' subj '\S2\EEG\VT_TFA_sumData_P' subj '_2.mat'])
     
    %Load and concatenate behaviour
    load([exp.procBehPath, 'PreprocData_' subj '-.mat']) 
    
    choices = allBeh.allChoices(logical(allBeh.fullTrials));
  
    %Compute lateralisation: LEFT - RIGHT hemisphere
    lhemi = trl_data(:,exp.movCluster(1,:),:);
    rhemi = trl_data(:,exp.movCluster(2,:),:);
   
    clear baselines; 
    baselines = nanmean(trl_data(:,:,0.4*200:0.5*200),3);

    dbdata = trl_data - baselines;
    allPbaselines(id,:) = nanmean(baselines);
   
    leftChoice_db(id,:,:) = squeeze(nanmean(nanmean(dbdata(choices == 0,exp.movCluster(1,:),:) - dbdata(choices == 0,exp.movCluster(2,:),:))));
    rightChoice_db(id,:,:) = squeeze(nanmean(nanmean(dbdata(choices ==1,exp.movCluster(1,:),:) - dbdata(choices== 1,exp.movCluster(2,:),:))));

    topoData_left(id,:,:,:) = nanmean(dbdata(choices == 0,:,:));
    topoData_right(id,:,:,:) = nanmean(dbdata(choices == 1,:,:));
end


%%
f = figure;

timeLim = [0.4,(4.5+0.4)]*200;
stdshade(leftChoice_db(:,timeLim(1):timeLim(2)), 0.15, 'k', '-',[],1) %hex2rgb(plt.brightPalettes{3}([2])
hold on;
stdshade(rightChoice_db(:,timeLim(1):timeLim(2)), 0.15, 'k', '--',[],1)

for i = 1:10
        sampleTimes(i) = [(0.1*200)+(0.4*i*200)]
end
ylim([-0.75,0.75])

xline(sampleTimes, '--k');
xlim([0,900])
xline([0.1*200])
xlabel('Time (s)')
ylabel({'MBL (dB) [Left-Right]'})

xticks(([0.1,0.5:0.4:4.5]*200));
xticklabels(([0,0.4,0.8,1.2,1.6,2,2.4,2.8,3.2,3.6,4,4.4]));
legend({'','Left response','', 'Right response'})
set(gca, 'fontsize', 11);

ax1 = gca();
% ADD second X-axis for sample labels
ax2 = axes('Position', get(ax1,'Position'), ...
    'XAxisLocation','top', ...
    'Color','none', ...
    'XColor','k',...
    'XLim', [ax1.XLim]);
ax2.YAxis.Visible = 'off';
ax2.XAxis.Visible = 'on';
ax2.XAxis.TickLength = [0,0];

ax2.XTick = [([0.1,0.5:0.4:4.5]*200)];
ax2.XTickLabel = {['Trial' ,'\newline', ' start'], 'S1', 'S2', 'S3', 'S4', 'S5', 'S6', 'S7', 'S8', 'S9', 'S10',''}
xticks(ax2.XTick);
xticklabels(ax2.XTickLabel)

set(gca, 'fontsize', 10);

% ADD topography 
% Inset topographies 
% Create smaller axes in top left, and plot on it

axes('Position',[.08 .65 .22 .22])
box off
load('C:\Users\elisa\Desktop\VolatilityTask\Analysis\functions\chanlocsBioSemi_128_noext.mat')

ttopoData = squeeze(nanmean(topoData_left - topoData_right));
ttopoData = nanmean(ttopoData(:,[(4.2+0.4)*200]:[(4.4+0.4)*200]),2);

VT_topoplot(ttopoData,chanlocs, ...
    'electrodes', 'off', ...
    'emarker', {'.','k',[],5},...
    'emarker2', {[exp.movCluster(1,:), exp.movCluster(2,:)], '*', 'k',3,0.5},...)
    'numcontour', 0, ...
    'plotrad',0.65, ...
    'headrad', 0.5)

caxis([-0.6 0.6])

f.Units = 'centimeters';
f.OuterPosition = [0.25 0.35 19 11];
exportgraphics(f,[exp.figPath, 'Fig2A.tiff'], 'Resolution', 600)

