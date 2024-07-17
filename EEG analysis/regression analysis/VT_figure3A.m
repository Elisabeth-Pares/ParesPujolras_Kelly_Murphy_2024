function [] = VT_figure3A(exp,par)

%% Panel A
for sub = 1:length(exp.sub_id)
    subj = exp.sub_id{sub};
    disp(['Processing P' num2str(sub) '...'])
    thisPath = ['P' subj '\S2\EEG\'];
    fullPath = ([exp.dataPath, thisPath]);
    load([fullPath, 'VT_ERP_sumData_P' subj '_2.mat']);
    
    % ==================================================================
    %% HP + LP filter
    % ==================================================================
    if par.hp == 1
        [bfilt,afilt] = butter(3, 1*2/200, 'high');   % hi-pass
        for t = 1:size(trl_data,1)
            for e = 1:size(trl_data,2)
                trl_data(t,e,:) = filtfilt(bfilt,afilt, trl_data(t,e,:));
            end
        end
    end
    
    
    % ==================================================================
    % Low-pass filter signal, if required
    % ==================================================================
    if par.lp == 1
        [bfilt,afilt] = butter(3, 6*2/200, 'low');   % lo-pass
        for t = 1:size(trl_data,1)
            for e = 1:size(trl_data,2)
                trl_data(t,e,:) = filtfilt(bfilt,afilt, trl_data(t,e,:));
            end
        end
    end
    
    allSum(sub,:,:) = squeeze(nanmean(trl_data));
end
% save(['VT_fullTrials_sumData_allP_hplpFiltered.mat'], 'allSum');
save([fullPath, 'VT_fullTrials_sumData_allP.mat'], 'allSum');

%%
allP_full.data = allSum;
f = figure;
ax1= axes();
colors = hex2rgb(plt.newPalettes{6});
ch = [87,1,3,19];

for i = 1:10
    sampleTimes(i) = [(0.1*200)+(0.4*i*200)]
    
end

hold on
f = figure
load('chanlocsBioSemi_128_noext.mat')

nan_dat = repmat(0,1,128)%size(allP.allDiff_allEv.fullTrials,3));
cc = 0
for c = ch(1:length(ch))
    cc = cc+1
    thisCol = colors(cc,:);
    
    stdshade(squeeze(allP_full.data(:,c,:)), 0.05, thisCol,'-', [],1);
    hold on;
    
end

xline(sampleTimes, '--k');
xlim([0,900])
xline([0.1*200])
xlabel('Time (s)')
ylabel('EEG amplitude (\muV/m^{2})')
xticks(([0.1,0.5:0.4:4.5]*200));
xticklabels(([0,0.4,0.8,1.2,1.6,2,2.4,2.8,3.2,3.6,4,4.4]));
set(gca, 'fontsize', 18);

% ADD second X-axis for sample labels
ax2 = axes('Position', get(ax1,'Position'), ...
    'XAxisLocation','top', ...
    'Color','none', ...
    'XColor','k');
xlim([0,900])
ax2.YAxis.Visible = 'off';
ax2.XTick = [([0.1,0.5:0.4:4.5]*200)];
ax2.XTickLabel = {'Mask', 'S1', 'S2', 'S3', 'S4', 'S5', 'S6', 'S7', 'S8', 'S9', 'S10',''}
xticks(ax2.XTick);
xticklabels(ax2.XTickLabel)

set(gca, 'fontsize', 18);


% Inset topographies create smaller axes in top left, and plot on it
ax3 = axes('Position',[.11 .69 .22 .22])

box off
nan_dat = repmat(0, 128,1)
colors = plt.newPalettes{6};
VT_topoplot(nan_dat,chanlocs, ...
    'style', 'contour',...
    'electrodes', 'on', ...
    'emarker', {'.','k',[],1},...
    'numcontour', 1)
for c = 1:4
    thisCol = hex2rgb(colors{c})
    hold on; VT_topoplot(nan_dat,chanlocs, ...
        'style', 'contour',...
        'electrodes', 'on', ...
        'emarker', {'.','k',[],5},...
        'emarker2', {ch(c), '*', thisCol,5,2},...)
        'numcontour', 1)
end

f.Units = 'centimeters';

f.OuterPosition = [0.25 0.25 31 16];

exportgraphics(f,[exp.figPath, 'Fig3A.tiff'], 'Resolution', 600)
