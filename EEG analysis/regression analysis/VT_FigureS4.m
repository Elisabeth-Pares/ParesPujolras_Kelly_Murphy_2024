%% Supplementary figure S4
function [] = VT_FigureS4(exp, par)

% Load single subject regression results
nsubj = 0;
for sub = exp.sub_id(1:end)
    nsubj = nsubj +1;
    disp(['Loading sub... ' num2str(nsubj)])
    
    ch = 1:128;%Whole scalp
    
    filename = (['VT_regression_full_ERP_rsq_HP_0_P' sub{1} '_2.mat'])
    
    load([exp.dataPath 'P' sub{1} '\' filename]);
    
    %Average over cluster of channels, for all samples
    s = [1:9];
    prior(nsubj,:,:) = mean(dat.prior.beta(ch,:,s),3);
    sprior(nsubj,:,:,:) = dat.prior.beta(ch,:,:);
    LLR(nsubj,:,:) = mean(dat.zLLR.beta(ch,:,s),3);
    LLRxSurprise(nsubj,:,:) = mean(dat.zLLRxSurprise.beta(ch,:,s),3);
end

%% Prepare data for fieldtrip analysis
clear permData;
for s = 1:20
    permData.prior{s} = squeeze(prior(s,:,:));
    permData.LLR{s}= squeeze(LLR(s,:,:));
    permData.LLRxSurprise{s} = squeeze(LLRxSurprise(s,:,:));
end

%% Plot figure
% Plot various channel subsets
preFrontalFocus = [93,81,80,...
    92,82,79,...
    91,83,78]; 

frontalFocus = [91,83,78,...
    90,84,77,...
    89,85,76,...
    88,86,75]

centralFocus = [1,...
    112,111,2,33,34,34,...
    3]

posteriorFocus = [4,32,31,30,21,17,18,5,...
    19,20];

%Panels A-D, G,J in Fig. S4
chanList = {exp.centroParietal, [98,87,66,88,86,75], centralFocus, posteriorFocus,[84,83,82], [1,2,3]};

par.clusterTime = [0 1];
par.averageChan = 1;

for ch = 1:6
    par.chans = chanList{ch};
    
    %% Set paramseters for cluster analysis
    cluster = VT_clusterpermutation_ERP(exp,permData, par)
    labels = {'prior', 'llr', 'llrxsurprise'}; %Regressor names
    
    f = figure;
    stdshade(squeeze(nanmean(prior(:,par.chans,1:end),2)), 0.25, hex2rgb('#04080f'),'-',[]); hold on;
    stdshade(squeeze(nanmean(LLR(:,par.chans,1:end),2)), 0.25, hex2rgb('#192bc2'),'-', []); hold on;
    stdshade(squeeze(nanmean(LLRxSurprise(:,par.chans,1:end),2)), 0.25, hex2rgb('#78c0e0'),'-',[]); hold on;
    
    yline([0 0]); set(gca, 'fontsize', 11)
    yline(0); 
    xticks([1,100,200]); xticklabels({'0','0.5','1'});
    xlabel('Time(s)'); ylabel('Regression coefficients (a.u.)')
    
    hold on;
    colors = {hex2rgb('#04080f'), hex2rgb('#192bc2'),hex2rgb('#78c0e0'), 'g'};
    ypos= -0.01;

    for stats = [1,2,3]
        thisField = labels{stats}
        thisCol = colors{stats};

        if isfield(cluster.(thisField), 'posclusters')
            for i = 1:length(cluster.(thisField).posclusters)
                if cluster.(thisField).posclusters(i).prob < 0.05
                    thisMark = '*';
                    set(gca, 'FontSize', 10);
                    plot(find(cluster.(thisField).posclusterslabelmat(1,:) ==i),ypos,'Color', thisCol, 'marker', thisMark, 'markersize', 5)
                    hold on;
                    ypos = ypos - 0.00125;

                end
            end
        end
   
        if isfield(cluster.(thisField), 'negclusters')
            for i = 1:length(cluster.(thisField).negclusters)
                if cluster.(thisField).negclusters(i).prob < 0.05
                    thisMark = '*';
                    set(gca, 'FontSize', 10);
                    plot(find(cluster.(thisField).negclusterslabelmat(1,:) ==i),ypos,'Color', thisCol, 'marker', thisMark, 'markersize', 5)
                    ypos = ypos - 0.00125;

                    hold on;
                end
            end
        end
        set(gca, 'FontSize', 10);
        yticks([-0.015:0.005:0.02])
    end
 
    ylims = [[-0.0125,0.02]];
    ylim(ylims)


    %Plot channel selection & cluster topo
    %Inset topographies create smaller axes in top left, and plot on it
    ax3 = axes('Position',[.3 .69 .22 .22]) 
    box off
    load('chanlocsBioSemi_128_noext.mat')
    
    nan_dat = repmat(0, 128,1)
    VT_topoplot(nan_dat,chanlocs, ...
        'style', 'contour',...
        'electrodes', 'on', ...
        'emarker', {'.','k',[],5},...
        'emarker2', {par.chans, '*', 'k',3,2},...)
        'numcontour', 1)
    
    f.Units = 'centimeters';
    f.OuterPosition = [0.25 0.35 8 10];
    
end
%%
%For panels E-H of Fig. S4
%Load processed data 
% load([fullPath, 'VT_fullTrials_sumData_allP_realigned_LP.mat']);

chanList = {[84,83,82], [1,2,3]};

f = figure;
ax1= axes();
colors = hex2rgb(plt.newPalettes{6});
clear sampleTimes;
for i = 1:11
    sampleTimes(i) = [[(0.1*200)+(0.4*(i-1)*200)]];
end

for ch = 1:length(chanList)
    par.chans = chanList{ch};
    
    f = figure
    load('chanlocsBioSemi_128_noext.mat')
    
    nan_dat = repmat(0,1,128)%size(allP.allDiff_allEv.fullTrials,3));
    cc = 0
    for c = par.chans(1:length(par.chans))
        cc = cc+1
        thisCol = colors(cc,:);
        
        stdshade(squeeze(allSum(:,c,:)), 0.05, thisCol,'-', [],1); hold on;
        amean=nanmean(squeeze(allSum(:,c,:)));
        astd=nanstd(squeeze(allSum(:,c,:)))/sqrt(12); % to get sem shading
        
        errorup = [amean+astd];
        errordown = [(amean-astd)];
        for s = 1:11
            x = round(sampleTimes(s):1:sampleTimes(s)+76);
            thisError = [errorup(x), fliplr(errordown(x))];
            fill([x fliplr(x)],thisError, thisCol, 'FaceAlpha', 0.05,'linestyle','none'); hold on;
        end
        hold on;        hold on;
    end
    
    xline(sampleTimes, '--k');
    xlim([1,900])
    xline([0.1*200])
    xlabel('Time (s)')
    ylabel('EEG amplitude (\muV/m^{2})')
    xticks(([0.1,0.5:0.4:4.5]*200));
    xticklabels({'Mask','S1', 'S2', 'S3', 'S4', 'S5', 'S6', 'S7', 'S8', 'S9', 'S10',''});
    set(gca, 'fontsize', 18);
    
 
    % Inset topographies create smaller axes in top left, and plot on it
    ax3 = axes('Position',[.11 .69 .22 .22])
    
    box off
    nan_dat = repmat(0, 128,1)
    % colors = jet(128)%plt.newPalettes{6};
    VT_topoplot(nan_dat,chanlocs, ...
        'style', 'contour',...
        'electrodes', 'on', ...
        'emarker', {'.','k',[],1},...
        'numcontour', 1)
    
    cc = 0;
    for c = par.chans(1:end)
        cc = cc+1;
        hold on; VT_topoplot(nan_dat,chanlocs, ...
            'style', 'contour',...
            'electrodes', 'on', ...
            'emarker', {'.','k',[],5},...
            'emarker2', {c, '*', colors(cc,:),5,2},...)
            'numcontour', 1)
    end
    
    f.Units = 'centimeters';
    f.OuterPosition = [0.25 0.25 31 16];
    
end
%% Topofigure

for ch = 1:length(chanList)
    par.chans = chanList{ch};
    figure

lastSamp = sampleTimes(10); 

lastAmp = squeeze(nanmean(nanmean(allSum(:,:,lastSamp:lastSamp+0.4*200),3),1)); 
firstAmp = squeeze(nanmean(nanmean(allSum(:,:,0.1*200:(0.1*200)+0.4*200),3),1)); 

topoDiff = lastAmp - firstAmp; 

VT_topoplot(topoDiff,chanlocs, ...
        'electrodes', 'on', ...
        'emarker', {'.','k',[],5},...
        'numcontour', 1)


for c = par.chans(1:end)
    hold on; VT_topoplot(topoDiff,chanlocs, ...
        'style', 'contour',...
        'electrodes', 'on', ...
        'emarker', {'.','k',[],5},...
        'emarker2', {c, '*', 'k',5,2},...)
        'numcontour', 1)
end
end
colorbar
