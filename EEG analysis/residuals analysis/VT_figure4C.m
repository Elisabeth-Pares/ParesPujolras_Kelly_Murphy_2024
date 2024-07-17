function [] = VT_figure4C(exp, par)

%% Plot
allP = [];
clear permData; 
maxS = 10;
sub = 0;
%Load data 
for subj = exp.sub_id(1:end)
    sub = sub+1;
    thisPath = ['P' subj{1}];
    fullPath = ([exp.dataPath, thisPath]);
    
    load([fullPath, '\VT_cICAEEGreg_allCh_CPPRes2Beta_' par.regName '_allRes_P' subj{1} '.mat']);
    
    for r = 1:length(par.res)
        thisRes = par.res{r};
        
        allP.effEv.(thisRes)(sub,:) = squeeze(nanmean(dat.effEv.(thisRes).beta(:,:,1:end),3));
        allP.effEvxCPPres_all.(thisRes)(sub,:,:) = dat.effEvxCPPres.(thisRes).beta(:,:,1:end);
    end
end
%%
% Cluster analysis 
maxS = 10
for s = 1:20
    permData.effEvxCPPres{s} = squeeze(nanmean(allP.effEvxCPPres_all.timeResolved_frontPostCluster(s,:,1:maxS),3));
    permData.effEvxCPPres_fminusb{s}  = squeeze(nanmean(allP.effEvxCPPres_all.timeResolved_frontalCluster(s,:,1:maxS),3)) - squeeze(nanmean(allP.effEvxCPPres_all.timeResolved_posteriorCluster(s,:,1:maxS),3)); 
    permData.effEvxCPPres_eml_fp{s}  = squeeze(nanmean(allP.effEvxCPPres_all.timeResolved_frontPostCluster(s,:,1:7),3)) - squeeze(nanmean(allP.effEvxCPPres_all.timeResolved_frontPostCluster(s,:,8:maxS),3));
    permData.effEvxCPPres_eml_f{s}  = squeeze(nanmean(allP.effEvxCPPres_all.timeResolved_frontalCluster(s,:,1:7),3)) - squeeze(nanmean(allP.effEvxCPPres_all.timeResolved_frontalCluster(s,:,8:maxS),3));
end

par.lat = [0.1 0.5];
[cluster] = VT_TFA_cluster_CPPxBeta(permData, par, 'effEvxCPPres')
[cluster_fmb] = VT_TFA_cluster_CPPxBeta(permData, par, 'effEvxCPPres_fminusb')
[cluster_eml_fp] = VT_TFA_cluster_CPPxBeta(permData, par, 'effEvxCPPres_eml_fp')
[cluster_eml_f] = VT_TFA_cluster_CPPxBeta(permData, par, 'effEvxCPPres_eml_f')


% Plot cluster results
labels = {'stat1'};
colors = [{hex2rgb('#883bd5')}, {hex2rgb('#b88ae6')}, {hex2rgb('#631db0')}]
%%
f = figure
f.Units = 'centimeters';
f.OuterPosition = [0.25 0.35 32 9];

hold on; 

subplot(1,4,1)
stdshade(squeeze(nanmean(allP.effEvxCPPres_all.timeResolved_frontPostCluster(:,:,par.minS:par.maxS),3)), 0.15, colors{1}, '-', [],1)
xticks([0,100,200]); xticklabels({'0','0.5','1'})
xlabel('Time from sample (s)')
set(gca, 'fontsize', 10)
ylabel('Regression coefficient (a.u.)')
% title(par.res)
yline([0,0], '--k')
ylim([-0.009,0.020]);

ypos= -0.005

for stats = [1]
    patch([par.lat(1) par.lat(1) par.lat(2) par.lat(2)]*200,[-0.009,0.020, 0.020, -0.009], 'k', 'FaceAlpha',0.05, 'EdgeColor', 'none')
    thisField = labels{stats}%['stat' num2str(stats)];
    thisCol = colors{stats}; hold on; 
    
    if isfield(cluster.(thisField), 'posclusters')
        if cluster.(thisField).posclusters(1).prob < 0.05
            thisMark = '*';
        else
            continue;
        end
    else 
        continue; 
    end
    ypos = ypos - 0.001;
    set(gca, 'fontsize', 10)
    plot(find(cluster.(thisField).posclusterslabelmat(1,:) ==1)+par.lat(1)*200,ypos,'Color', thisCol, 'marker', thisMark, 'markersize', 5)
    hold on;
 
end
set(gca, 'fontsize', 10)

subplot(1,4,2)
stdshade(squeeze(nanmean(allP.effEvxCPPres_all.timeResolved_frontalCluster(:,:,par.minS:par.maxS),3)), 0.15, colors{2}, '-', [],1); hold on
stdshade(squeeze(nanmean(allP.effEvxCPPres_all.timeResolved_posteriorCluster(:,:,par.minS:par.maxS),3)), 0.15, colors{3}, '-', [],1)
xticks([0,100,200]); xticklabels({'0','0.5','1'})
xlabel('Time from sample (s)')
ylabel('Regression coefficient (a.u.)')
yline([0,0], '--k')
ylim([-0.009,0.020]);

ypos= -0.005;

for stats = [1]
    patch([par.lat(1) par.lat(1) par.lat(2) par.lat(2)]*200,[-0.009,0.020, 0.020, -0.009], 'k', 'FaceAlpha',0.05, 'EdgeColor', 'none') 
    thisField = labels{stats}%['stat' num2str(stats)];
    thisCol = colors{stats}; hold on; 
    
    if isfield(cluster_fmb.(thisField), 'posclusters')
        if cluster_fmb.(thisField).posclusters(1).prob < 0.05
            thisMark = '*';
        else
            continue;
        end
    else 
        continue; 
    end
    ypos = ypos - 0.001;
    set(gca, 'fontsize', 10)
    plot(find(cluster_fmb.(thisField).posclusterslabelmat(1,:) ==1)+par.lat(1)*200,ypos,'Color', thisCol, 'marker', thisMark, 'markersize', 5)
    hold on;
    
%     patch([0.25 0.4 0.6 0.6 0.4],[0.020 -0.009 -0.009 0.020 0.020],'k','FaceAlpha',0.05);

end
set(gca, 'fontsize', 10)
legend({'','Anterior', '', 'Posterior'}, 'location', 'northwest')

subplot(1,4,3)
hold on;
stdshade(squeeze(nanmean(allP.effEvxCPPres_all.timeResolved_frontPostCluster(:,:,par.minS:7),3)), 0.15, colors{1}, '-', [],1); 
hold on; 
stdshade(squeeze(nanmean(allP.effEvxCPPres_all.timeResolved_frontPostCluster(:,:,8:10),3)), 0.15, colors{1}, '--', [],1)
box on 

xticks([0,100,200]); xticklabels({'0','0.5','1'})
xlabel('Time from sample (s)')
set(gca, 'fontsize', 10)
ylabel('Regression coefficient (a.u.)')
% title(par.res)
yline([0,0], '--k')
ylim([-0.009,0.020]);
box on 

for stats = [1]
    patch([par.lat(1) par.lat(1) par.lat(2) par.lat(2)]*200,[-0.009,0.020, 0.020, -0.009], 'k', 'FaceAlpha',0.05, 'EdgeColor', 'none')
    thisField = labels{stats}%['stat' num2str(stats)];
    thisCol = colors{stats}; hold on; 
    
    if isfield(cluster_eml_fp.(thisField), 'posclusters')
        if cluster_eml_fp.(thisField).posclusters(1).prob < 0.05
            thisMark = '*';
        else
            continue;
        end
    else 
        continue; 
    end
    ypos = ypos - 0.001;
    set(gca, 'fontsize', 10)
    plot(find(cluster_eml_fp.(thisField).posclusterslabelmat(1,:) ==1)+par.lat(1)*200,ypos,'Color', thisCol, 'marker', thisMark, 'markersize', 5)
    hold on;
 
end
legend({'','Early [1-7]', '', 'Late [8-10]'}, 'location', 'northwest');

set(gca, 'fontsize', 10)
subplot(1,4,4)
stdshade(squeeze(nanmean(allP.effEvxCPPres_all.timeResolved_frontalCluster(:,:,par.minS:7),3)), 0.15, colors{2}, '-', [],1); 
hold on; 
stdshade(squeeze(nanmean(allP.effEvxCPPres_all.timeResolved_frontalCluster(:,:,8:10),3)), 0.15, colors{2}, '--', [],1)

xticks([0,100,200]); xticklabels({'0','0.5','1'})
xlabel('Time from sample (s)')
set(gca, 'fontsize', 10)
ylabel('Regression coefficient (a.u.)')
% title(par.res)
yline([0,0], '--k')
ylim([-0.009,0.020]);
ypos= -0.005

for stats = [1]
    patch([par.lat(1) par.lat(1) par.lat(2) par.lat(2)]*200,[-0.009,0.020, 0.020, -0.009], 'k', 'FaceAlpha',0.05, 'EdgeColor', 'none')
    thisField = labels{stats}%['stat' num2str(stats)];
    thisCol = colors{stats}; hold on; 
    
    if isfield(cluster_eml_f.(thisField), 'posclusters')
        if cluster_eml_f.(thisField).posclusters(1).prob < 0.05
            thisMark = '*';
        else
            continue;
        end
    else 
        continue; 
    end
    ypos = ypos - 0.001;
    set(gca, 'fontsize', 10)
    plot(find(cluster_eml_f.(thisField).posclusterslabelmat(1,:) ==1)+par.lat(1)*200,ypos,'Color', thisCol, 'marker', thisMark, 'markersize', 5)
    hold on;
 
end
legend({'','Early [1-7]', '', 'Late [8-10]'}, 'location', 'northwest');

% Topoplots
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


box on 
axes('Position',[.145 .67 .23 .23])
load('chanlocsBioSemi_128_noext.mat')
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

axes('Position',[.35 .67 .23 .23 ])
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
box on 
axes('Position',[.555 .67 .23 .23])
load('chanlocsBioSemi_128_noext.mat')
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
box on 

for r = 2
    thisCol = colors{r}; 
    thisRes = par.res{r};
if strcmp(thisRes, 'timeResolved_posteriorCluster')
    chanSel = exp.posteriorCluster;
elseif strcmp(thisRes, 'timeResolved_frontalCluster')
    chanSel = exp.frontalCluster;
elseif strcmp(thisRes, 'timeResolved_frontPostCluster')
    chanSel = exp.frontBackCluster;
end
box on 
axes('Position',[.76 .67 .23 .23])
load('chanlocsBioSemi_128_noext.mat')
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
f.OuterPosition = [0.25 0.35 32 9];
exportgraphics(f,[exp.figPath, 'Fig4C_' par.thisModel '_Final.tiff'], 'Resolution', 600)

end