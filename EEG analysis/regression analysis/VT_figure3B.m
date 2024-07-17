%% Figure 3B
function [] = VT_figure3B(exp, par)

% Load single subject regression results
nsubj = 0;
for sub = exp.sub_id(1:end)
    nsubj = nsubj +1;
    disp(['Loading sub... ' num2str(nsubj)])
    
    ch = 1:128;%Whole scalp
    
    filename = (['VT_regression_prior_zllr_zllrxsurprise_ERP_allT_subtr_parietalrsq0_P' sub{1}, '_2_LOG.mat'])
    
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

%% Set paramseters for cluster analysis
cluster = VT_clusterpermutation_ERP(exp,permData, par)
labels = {'prior', 'llr', 'llrxsurprise'}; %Regressor names

%% Plot figure
f = figure;
ylim([-0.008,0.020])

patch([0.25,0.25,0.3,0.3]*200,[-0.008,0.020, 0.020, -0.008], 'k', 'FaceAlpha',0.05, 'EdgeColor', 'none'); hold on;
stdshade(squeeze(nanmean(prior(:,par.chans,1:end),2)), 0.25, hex2rgb('#04080f'),'-',[]); hold on;
stdshade(squeeze(nanmean(LLR(:,par.chans,1:end),2)), 0.25, hex2rgb('#192bc2'),'-', []); hold on;
stdshade(squeeze(nanmean(LLRxSurprise(:,par.chans,1:end),2)), 0.25, hex2rgb('#78c0e0'),'-',[]); hold on;

yline([0 0]); set(gca, 'fontsize', 15)
yline(0); %title('LLR*uncertainty'); %ylim([-0.01,0.065]);
xticks([1,100,200]); xticklabels({'0','0.5','1'});
xlabel('Time(s)'); ylabel('Regression coefficients (a.u.)')

hold on;
colors = {hex2rgb('#04080f'), hex2rgb('#192bc2'),hex2rgb('#78c0e0'), 'g'};
ypos= -0.005

for stats = [1,2,3]
    thisField = labels{stats}
    thisCol = colors{stats};
    if isfield(cluster.(thisField), 'posclusters')
        for i = 1:length(cluster.(thisField).posclusters)
            if cluster.(thisField).posclusters(i).prob < 0.05
                thisMark = '*';
                        ypos = ypos - 0.001;
                        plot(find((cluster.(thisField).posclusterslabelmat(1,:) ==i)),ypos,'Color', thisCol, 'marker', thisMark, 'markersize', 5)

            else
                continue;
            end
        end
        hold on;
    end
    
    if isfield(cluster.(thisField), 'negclusters')
        for i = 1:length(cluster.(thisField).negclusters)
            if cluster.(thisField).negclusters(i).prob < 0.05
                thisMark = '*';
                ypos = ypos - 0.001;
                plot(find((cluster.(thisField).negclusterslabelmat(1,:) ==i)),ypos,'Color', thisCol, 'marker', thisMark, 'markersize', 5)

            else
                continue;
            end
        end
        hold on;
    end
end

set(gca, 'fontsize', 10)
box on

f.Units = 'centimeters';
f.OuterPosition = [0.25 0.35 8 10];
set(gca, 'fontsize', 10)

exportgraphics(f,[exp.figPath, 'Fig4B.tiff'], 'Resolution', 600)
