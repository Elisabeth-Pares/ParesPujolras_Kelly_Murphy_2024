function [] = VT_figure4A(exp, par)
%% Figure 4. Panel A - Cluster analysis of residual effects
load([exp.dataPath, 'GA_CPPxRes_regression.mat'], 'GA') %Output from VT_eeg2behRegression_timeRes_v2

clear permData
maxS = 10
for s = 1:20
    permData.LLRxRes{s} = squeeze(nanmean(GA.(par.thisModel).timeResolved_frontPostCluster.Res(s,1:maxS,:),2))';
    permData.LLRxRes_f{s} = squeeze(nanmean(GA.(par.thisModel).timeResolved_frontalCluster.Res(s,1:maxS,:),2))';
    permData.LLRxRes_b{s} = squeeze(nanmean(GA.(par.thisModel).timeResolved_posteriorCluster.Res(s,1:maxS,:),2))';
    permData.LLRxRes_fminb{s} = squeeze(nanmean(GA.(par.thisModel).timeResolved_frontalCluster.Res(s,1:maxS,:),2))'-squeeze(nanmean(GA.(par.thisModel).timeResolved_posteriorCluster.Res(s,1:maxS,:),2))'; 
end

par.clusterTime = [0.15 0.5];
cluster = VT_clusterpermutation_LLRxRes(exp,permData, par, 'LLRxRes')
cluster_fminb = VT_clusterpermutation_LLRxRes(exp,permData, par, 'LLRxRes_fminb')

%% Plot residual effects of front & back, then plot front & back separately 

f = figure;
colors = [{hex2rgb('#883bd5')}, {hex2rgb('#b88ae6')}, {hex2rgb('#631db0')}]

subplot(1,3,1)
patch([par.clusterTime(1) par.clusterTime(1) par.clusterTime(2) par.clusterTime(2)]*200,[-0.05,0.1 0.1 -0.05], 'k', 'FaceAlpha',0.05, 'EdgeColor', 'none')
hold on;
stdshade(squeeze(nanmean(GA.(par.thisModel).timeResolved_frontPostCluster.Res(:,1:maxS,:),2)),0.15,colors{1},'-', [], 1)

% stdshade(squeeze(nanmean(GA.(par.thisModel).timeResolved_c.Res(:,1:maxS,:),2)),0.15,colors{1},'-', [], 1)
yline([0,0], '--k');
ylabel('LLR weight modulation (a.u.)')
ylim([-0.05,0.075])

hold on
ypos= -0.04
labels = {'LLRxRes'}
for stats = [1]
    thisField = labels{stats}%['stat' num2str(stats)];
    thisCol = colors{stats}; hold on;
    
    if isfield(cluster.(thisField), 'posclusters')
        if cluster.(thisField).posclusters(1).prob < 0.05
            thisMark = '*';
        else
            continue
        end
    else
        continue;
    end
    ypos = ypos - 0.001;
    set(gca, 'fontsize', 10)
    hold on;
    subplot(1,3,1)
    plot(find(cluster.(thisField).posclusterslabelmat(1,:) ==1)+par.clusterTime(1)*200,ypos,'Color', thisCol, 'marker', thisMark, 'markersize', 5)
    hold on;
 
end
set(gca, 'fontsize', 10)


yline([0,0], '--k')
xticks([1,100,200]); xticklabels({'0', '0.5', '1'})
xlabel('Time from sample (s)')
set(gca, 'fontsize', 10)
box on 
hold on
axes('Position',[.18 .67 .23 .23])
load('chanlocsBioSemi_128_noext.mat')

for r = 1
    thisCol = colors{r}; 
    thisRes = par.res{r};
if strcmp(thisRes, 'timeResolved_posteriorCluster')
    chanSel = exp.posteriorCluster;
elseif strcmp(thisRes, 'timeResolved_frontalCluster')
    chanSel = exp.frontalCluster;
elseif strcmp(thisRes, 'timeResolved_frontPostCluster')
    chanSel = exp.frontBackCluster;
end

nan_dat = repmat(0,1,128)%size(allP.allDiff_allEv.fullTrials,3));
VT_topoplot(nan_dat,chanlocs(1:128), ...
    'style', 'contour',...
    'electrodes', 'off', ...
    'emarker2', {chanSel, '*', thisCol,3.5,0.5},...)
    'emarker', {'.','k',[],5},...
    'numcontour', 0, ...
    'plotrad',0.65,...
    'headrad', 0.5)
hold on; 
end

subplot(1,3,2)
patch([par.clusterTime(1) par.clusterTime(1) par.clusterTime(2) par.clusterTime(2)]*200,[-0.05,0.1 0.1 -0.05], 'k', 'FaceAlpha',0.05, 'EdgeColor', 'none')
hold on; 
stdshade(squeeze(nanmean(GA.(par.thisModel).timeResolved_frontalCluster.Res(:,1:maxS,:),2)),0.15,colors{2},'-', [], 1); hold on; %
stdshade(squeeze(nanmean(GA.(par.thisModel).timeResolved_posteriorCluster.Res(:,1:maxS,:),2)),0.15, colors{3},'-', [], 1); hold on
% ylabel('LLR weight modulation(a.u.)')
ylim([-0.05,0.075])
yline([0,0], '--k');
legend({'','','Anterior', '', 'Posterior'}, 'location', 'northwest')
box on; 
hold on
labels = {'LLRxRes'}
for stats = [1]
    thisField = labels{stats}%['stat' num2str(stats)];
    thisCol = colors{stats}; hold on;
    
    if isfield(cluster_fminb.(thisField), 'posclusters')
        if cluster_fminb.(thisField).posclusters(1).prob < 0.05
            thisMark = '*';
        else
            continue
        end
    else
        continue;
    end
    ypos = ypos - 0.001;
    set(gca, 'fontsize', 10)
    hold on;
    subplot(1,3,2)
    plot(find(cluster_fminb.(thisField).posclusterslabelmat(1,:) ==1)+par.clusterTime(1)*200,ypos,'Color', thisCol, 'marker', thisMark, 'markersize', 5)
    hold on;
       
end
xticks([1,100,200]); xticklabels({'0', '0.5', '1'})
xlabel('Time from sample (s)')
set(gca, 'fontsize', 10)
hold on
set(gca, 'fontsize', 10)

axes('Position',[.46 .67 .23 .23])
load('chanlocsBioSemi_128_noext.mat')

for r = 2:3
    thisCol = colors{r}; 
    thisRes = par.res{r};
if strcmp(thisRes, 'timeResolved_posteriorCluster')
    chanSel = exp.posteriorCluster;
elseif strcmp(thisRes, 'timeResolved_frontalCluster')
    chanSel = exp.frontalCluster;
elseif strcmp(thisRes, 'timeResolved_frontPostCluster')
    chanSel = exp.frontBackCluster;
end

nan_dat = repmat(0,1,128)%size(allP.allDiff_allEv.fullTrials,3));
VT_topoplot(nan_dat,chanlocs(1:128), ...
    'style', 'contour',...
    'electrodes', 'off', ...
    'emarker2', {chanSel, '*', thisCol,3.5,0.5},...)
    'emarker', {'.','k',[],5},...
    'numcontour', 0, ...
    'plotrad',0.65,...
    'headrad', 0.5)
hold on; 
end

subplot(1,3,3)
% patch([par.clusterTime(1) par.clusterTime(1) par.clusterTime(2) par.clusterTime(2)]*200,[-0.05,0.1 0.1 -0.05], 'k', 'FaceAlpha',0.05, 'EdgeColor', 'none')
stdshade(squeeze(nanmean(GA.(par.thisModel).timeResolved_frontPostCluster.Res(:,1:7,:),2)),0.15,colors{1},'--', [], 1); hold on; %
stdshade(squeeze(nanmean(GA.(par.thisModel).timeResolved_frontPostCluster.Res(:,8:10,:),2)),0.15, colors{1},'-', [], 1); hold on
patch([par.clusterTime(1) par.clusterTime(1) par.clusterTime(2) par.clusterTime(2)]*200,[-0.05,0.1 0.1 -0.05], 'k', 'FaceAlpha',0.05, 'EdgeColor', 'none')
box on; 
% ylabel('LLR weight modulation(a.u.)')
ylim([-0.05,0.075])
yline([0,0], '--k');
legend({'','', 'Early [1-7]', '', 'Late [8-10]'}, 'location', 'northwest')

axes('Position',[.75 .67 .23 .23])
for r = 1
    thisCol = colors{r}; 
    thisRes = par.res{r};
if strcmp(thisRes, 'timeResolved_posteriorCluster')
    chanSel = exp.posteriorCluster;
elseif strcmp(thisRes, 'timeResolved_frontalCluster')
    chanSel = exp.frontalCluster;
elseif strcmp(thisRes, 'timeResolved_frontPostCluster')
    chanSel = exp.frontBackCluster;
end

nan_dat = repmat(0,1,128)%size(allP.allDiff_allEv.fullTrials,3));
VT_topoplot(nan_dat,chanlocs(1:128), ...
    'style', 'contour',...
    'electrodes', 'off', ...
    'emarker2', {chanSel, '*', thisCol,3.5,0.5},...)
    'emarker', {'.','k',[],5},...
    'numcontour', 0, ...
    'plotrad',0.65,...
    'headrad', 0.5)
hold on; 
end

f.Units = 'centimeters';
f.OuterPosition = [0.25 0.35 24 9];

set(gca, 'fontsize', 10)

% exportgraphics(f,[exp.figPath, 'Fig4A_' par.thisModel '_Final.tiff'], 'Resolution', 600)