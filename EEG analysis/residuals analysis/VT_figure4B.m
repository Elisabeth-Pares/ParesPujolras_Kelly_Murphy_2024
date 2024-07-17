function [] = VT_figure4B(exp, par)
%% Figure 4. Panel B - Cluster analysis of residual effects
load([exp.dataPath, 'GA_BetaRes_choice_timeRes_regression.mat'], 'GA') %Output from VT_eeg2behRegression

clear permData
maxS = 9
for s = 1:20
    permData.tfaRes{s} = squeeze(nanmean(GA.BetaRes_choice_timeRes.tfa_res.Res(s,1:maxS,:),2))';
    permData.tfaRes_lme{s} = squeeze(nanmean(GA.BetaRes_choice_timeRes.tfa_res.Res(s,7:9,:),2))'-squeeze(nanmean(GA.BetaRes_choice_timeRes.tfa_res.Res(s,1:6,:),2))';
end

par.clusterTime = [0,1];
cluster = VT_clusterpermutation_LLRxRes(exp,permData, par, 'tfaRes')
cluster_eminl = VT_clusterpermutation_LLRxRes(exp,permData, par, 'tfaRes_lme')

%% Plot residual effects of front & back, then plot front & back separately 

f = figure;
colors = [{hex2rgb('#883bd5')}, {hex2rgb('#b88ae6')}, {hex2rgb('#631db0')}]

subplot(1,2,1)
hold on;
stdshade(squeeze(nanmean(GA.BetaRes_choice_timeRes.tfa_res.Res(:,1:maxS,:),2)),0.15,colors{1},'-', [], 1)
ylim([-0.045,0.12])
yline([0,0], '--k');
ylabel([{'P(choice) modulation', 'by MBL residuals (a.u.)'}])

yline([0,0], '--k')
hold on
ypos= -0.035;
labels = {'LLRxRes'}
for stats = [1]
    thisField = labels{stats}
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
    subplot(1,2,1)
    plot(find(cluster.(thisField).posclusterslabelmat(1,:) ==1)+par.clusterTime(1)*200,ypos,'Color', thisCol, 'marker', thisMark, 'markersize', 5)
    hold on;
    
    
end
set(gca, 'fontsize', 10)
yline([0,0], '--k')
xticks([1,100,200]); xticklabels({'0', '0.5', '1'})
xlabel('Time from sample (s)')
set(gca, 'fontsize', 10)
box on 

subplot(1,2,2)
stdshade(squeeze(nanmean(GA.BetaRes_choice_timeRes.tfa_res.Res(:,1:6,:),2)),0.15,colors{1},'--', [], 1); hold on
stdshade(squeeze(nanmean(GA.BetaRes_choice_timeRes.tfa_res.Res(:,7:9,:),2)),0.15,colors{1},'-', [], 1)
ylim([-0.045,0.12])
yline([0,0], '--k')
hold on
ypos= -0.035;
yline([0,0], '--k');

hold on
labels = {'LLRxRes'};
for stats = [1]
    thisField = labels{stats}%['stat' num2str(stats)];
    thisCol = colors{stats}; hold on;
    
    if isfield(cluster_eminl.(thisField), 'posclusters')
        for i = 1:length([cluster_eminl.(thisField).posclusters.prob])
            if cluster_eminl.(thisField).posclusters(i).prob < 0.05
                thisMark = '*';
            else
                continue
            end
            
            ypos = ypos - 0.001;
            set(gca, 'fontsize', 10)
            hold on;
            plot(find(cluster_eminl.(thisField).posclusterslabelmat(1,:) ==i)+par.clusterTime(1)*200,ypos,'Color', thisCol, 'marker', thisMark, 'markersize', 5)
            hold on;
        end
    end
       
end
xticks([1,100,200]); xticklabels({'0', '0.5', '1'})
xlabel('Time from sample (s)')
set(gca, 'fontsize', 10)
hold on
legend({'','Early [2-7]', '', 'Late [8-10]'}, 'location', 'northwest')

f.Units = 'centimeters';
f.OuterPosition = [0.25 0.35 16 9];

set(gca, 'fontsize', 10)

exportgraphics(f,[exp.figPath, 'Fig4B_' par.thisModel '_Final.tiff'], 'Resolution', 600)