%% Plot
f = figure
stdshade(squeeze(nanmean(nanmean(effPrior_s(:,exp.centroParietal,:,2:end),2),4)), 0.25, 'k','-',[]); hold on;
stdshade(squeeze(nanmean(nanmean(effPrior_s(:,exp.frontBackCluster,:,2:end),2),4)), 0.25, 'k','--',[]); hold on;

ylims = [[-0.0325,0.06]];
ylim(ylims)

yline([0 0]); set(gca, 'fontsize', 20);
yline(0); 
xticks([0,100,200]); xticklabels({'0','0.5','1'});
xlabel('Time(s)'); ylabel('Regression coefficients (a.u.)')

hold on;
colors = {'k', 'b', 'r', 'g'};
ypos= -0.0275;
labels = {'stat1'};
for stats = [1]
    patch([0.25,0.25,0.3,0.3]*200,[ylims(1),ylims(2), ylims(2), ylims(1)], 'k', 'FaceAlpha',0.05, 'EdgeColor', 'none')
    thisField = labels{stats}%['stat' num2str(stats)];
    thisCol = colors{stats}; hold on;
    
    if isfield(cluster.(thisField), 'posclusters') & ~isempty(cluster.(thisField).posclusters)
        for i = 1:length(cluster.(thisField).posclusters)
            if cluster.(thisField).posclusters(i).prob < 0.05
                thisMark = '*';
                set(gca, 'FontSize', 12);
                plot(find(cluster.(thisField).posclusterslabelmat(1,:) ==i),ypos,'Color', thisCol, 'marker', thisMark, 'markersize', 5)
                hold on;
            end
        end
    end
    
    
    if isfield(cluster.(thisField), 'negclusters') & ~isempty(cluster.(thisField).negclusters)
        for i = 1:length(cluster.(thisField).posclusters)
            
            if cluster.(thisField).negclusters(i).prob < 0.05
                thisMark = '*';
                set(gca, 'FontSize', 12);
                plot(find(cluster.(thisField).negclusterslabelmat(1,:) ==i),ypos,'Color', thisCol, 'marker', thisMark, 'markersize', 5)
                hold on;
            end
        end
    end
    set(gca, 'FontSize', 12);
            ypos = ypos - 0.001;

    
end


legend({'|\Delta\Psi|',''});
% set(gca, 'FontSize', 12);
% f.Units = 'centimeters';
% f.OuterPosition = [0.25 0.25 13 16];
% 
% end
f.Units = 'centimeters';
f.OuterPosition = [0.25 0.35 8 10];
set(gca, 'fontsize', 10)
exportgraphics(f,[exp.figPath, 'Fig3C.tiff'], 'Resolution', 600)