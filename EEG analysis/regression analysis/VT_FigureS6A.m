function [] = VT_FigureS6(exp, par)

%% Extract or load data
par.smpwin = [0 1.4];  % window for sample-wise analyses
par.trlwin = [-0.1 5.4];  % window for full-trial-wise analyses

par.lowFilt = 1;
par.baseline = 'preTrial'; %Options = 'preTrial'; 
par.baselineTime = [0.1];
gaDat = [];
par.maxLLR = 4.75;
par.dataType = 'ERP';

for nsubj = 1:20
    gaDat = VT_CPnoCP_ERP(nsubj,exp,par, gaDat)
end
save([exp.dataPath, 'VT_CPnoCP_' par.dataType 'gaData.mat'], 'gaDat');

%% Data for permutation test 
load([exp.dataPath, 'VT_CPnoCP_' par.dataType 'gaData.mat']);

clear permData;
for s = 1:20
    permData.cp{s} = squeeze(gaDat.cp.data(s,:)-nanmean(gaDat.cp.data(s,1:10),2));
    permData.nocp{s}= squeeze(gaDat.nocp.data(s,:)-nanmean(gaDat.nocp.data(s,1:10),2));
end 
    
par.clusterTime = 'all';%'all'; %Analysis period
par.chans = 1%exp.centroParietal;  %
par.averageChan = 0;
cluster = VT_clusterpermutation_CPnoCP(exp,permData, par)

%% Plot & save data
f = figure;
ff = fields(gaDat);
colors = flip(gray(6));
colors = colors(2:end,:);

subplot(1,3,[1,2])
for field = 1:length(ff);
    thisField = ff{field}
    baselines =  nanmean(gaDat.(thisField).data(:,1:10),2);
    stdshade(gaDat.(thisField).data-baselines, 0.15,colors(field,:), '-', [], 1);
    allPsi(field) = nanmean(gaDat.(thisField).mean_deltaPsi);
    allLLR(field) = nanmean(gaDat.(thisField).mean_LLR);
    allPsi_error(field) = nanstd(gaDat.(thisField).mean_deltaPsi)./sqrt(20);
    allLLR_error(field) = nanstd(gaDat.(thisField).mean_LLR)./sqrt(20);
    hold on;
end

hold on
xticks([1,100,200]); xticklabels({'0','0.5','1'}); xlabel('Time (s)');
yline([0,0], 'k')
ylabel('EEG amplitude (\muV/m^{2})')
set(gca, 'fontsize', 15)

ypos = -0.25;
thisField = 'stat2';
if isfield(cluster.(thisField), 'posclusters')
    if ~isempty(cluster.(thisField).posclusters)
        if cluster.(thisField).posclusters(1).prob < 0.05
            thisMark = '*';
            ypos = ypos - 0.001;
            set(gca, 'FontSize', 15);
            plot(find(cluster.(thisField).posclusterslabelmat(1,:) ==1),ypos,'Color', colors(1,:), 'marker', thisMark, 'markersize', 5)
            hold on;
        end
    end
end


if isfield(cluster.(thisField), 'negclusters')
    if ~isempty(cluster.(thisField).negclusters)
        
        if cluster.(thisField).negclusters(1).prob < 0.05
            
            thisMark = '*';
            ypos = ypos - 0.001;
            set(gca, 'FontSize', 15);
            plot(find(cluster.(thisField).negclusterslabelmat(1,:) ==1),ypos,'Color', colors(1,:), 'marker', thisMark, 'markersize', 5)
            hold on;
        end
    end
end
    ylim([-0.5,3])
    
axes('Position',[.10 .71 .2 .2])
load('chanlocsBioSemi_128_noext.mat')
nan_dat = repmat(0,1,128)%size(allP.allDiff_allEv.fullTrials,3));
VT_topoplot(nan_dat,chanlocs(1:128), ...
    'style', 'contour',...
    'electrodes', 'off', ...
    'emarker2', {exp.centroParietal, '*', 'k',3.5,0.5},...)
    'emarker', {'.','k',[],5},...
    'numcontour', 0, ...
    'plotrad',0.65,...
    'headrad', 0.5)
hold on; 

subplot(1,3,[3])
b = bar([allPsi], 'FaceColor', [0.9,0.9,0.9]);hold on;
b.FaceColor = 'flat';b.CData = colors(1:2,:);set(gca, 'fontsize', 15)
errorbar([1:2], allPsi, allPsi_error, '.k'); hold on;
ylabel(['|\Delta\Psi|']); xticks([1,2]); xticklabels({'CP', 'noCP'});


f.Units = 'centimeters';
f.OuterPosition = [5 5 16 12];
% 
% 
% 
% exportgraphics(f,[exp.figPath, 'Fig3D.tiff'], 'Resolution', 600)

%% Figure
%Plot for each subject 
figure
tiledlayout(5,4)
for s = 1:20
    nexttile;
    for field = 1:length(ff);
        thisField = ff{field}
        baselines =  0%nanmean(gaDat.(thisField).data(s,1:10),2);
        plot(squeeze(gaDat.(thisField).data(s,:)-baselines), 'Color',colors(field,:));
        hold on;
    end
    title(['Sub ' num2str(s)])
end