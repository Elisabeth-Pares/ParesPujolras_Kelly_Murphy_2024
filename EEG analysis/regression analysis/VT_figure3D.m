function [] = VT_figure3D(exp, par)

%% Extract or load data
par.smpwin = [0 1.4];  % window for sample-wise analyses
par.trlwin = [-0.1 5.4];  % window for full-trial-wise analyses

par.lowFilt = 1;
par.baseline = 'preTrial'; %Options = 'preTrial'; 
par.baselineTime = [0.1];
gaDat = [];
par.maxLLR = 4.75;


if strcmp(par.dataType, 'TFA')
    for nsubj = 1:20
        gaDat = VT_psiQuartiles_TFA(nsubj,exp,par, gaDat)
    end
else
    for nsubj = 1:20
        gaDat = VT_psiQuartiles_ERP(nsubj,exp,par, gaDat)
    end
end
save([exp.dataPath, 'VT_medianSplit_' par.dataType 'gaData_groupN' num2str(par.nGroups) '.mat'], 'gaDat');

%% Plot & save data
f = figure;
ff = fields(gaDat);
colors = flip(gray(6));
colors = colors(2:end,:);

if strcmp(par.dataType, 'TFA')
    subplot(2,4,[1:2,5:6]);

    for field = 1:length(ff);
        thisField = ff{field}
        stdshade(gaDat.(thisField).data, 0.05,colors(field,:), '-', [], 1);
        allPsi(field) = nanmean(gaDat.(thisField).mean_Psi);
        allPsi_error(field) = nanstd(gaDat.(thisField).mean_Psi)./sqrt(20);
        hold on;
    end
    set(gca, 'fontsize', 11)
    ylabel({['MBL (dB)'], '\leftarrow Right             Left \rightarrow'})

    ylim([-0.4,0.4]);
    xticks([0,100,200]); xticklabels({'0','0.5','1'}); xlabel('Time (s)');

    subplot(2,4,[3:4,7:8]);

    set(gca, 'fontsize', 11)
    b = bar([allPsi], 'FaceColor', [0.9,0.9,0.9]);hold on; 
    b.FaceColor = 'flat';b.CData = colors(1:4,:);set(gca, 'fontsize', 11, 'YDir', 'reverse')
    errorbar([1:4], allPsi, allPsi_error, '.k'); hold on;
    xlabel(['\leftarrow Left                  Right \rightarrow']);
    xticks([])
       title(['\Psi'])
    set(gca, 'fontsize', 11);
    hold on;
    
    yline([0,0], 'k')
    set(gca, 'fontsize', 11);
 hold on
 
    f.Units = 'centimeters';
    f.OuterPosition = [10 10 19 9];

else
    subplot(2,3,[1,2,4,5])
    for field = 1:length(ff);
        thisField = ff{field}
        baselines =  nanmean(gaDat.(thisField).data(:,1:10),2);
        stdshade(gaDat.(thisField).data-baselines, 0.05,colors(field,:), '-', [], 1);
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
    set(gca, 'fontsize', 10)
    
    
    subplot(2,3,[3])
    b = bar([allPsi], 'FaceColor', [0.9,0.9,0.9]);hold on;
    b.FaceColor = 'flat';b.CData = colors(1:4,:);set(gca, 'fontsize', 10)
    errorbar([1:4], allPsi, allPsi_error, '.k'); hold on;
    ylabel(['|\Delta\Psi|']); xticks([]); xlabel('Quartile');
    
    subplot(2,3,[6])
    b = bar([allLLR], 'FaceColor', [0.9,0.9,0.9]);hold on;
    b.FaceColor = 'flat';b.CData = colors(1:4,:);set(gca, 'fontsize', 10)
    errorbar([1:4], allLLR, allLLR_error, '.k'); hold on;
    ylabel(['|LLR|']); xticks([]); xlabel('Quartile');
    
    f.Units = 'centimeters';
    f.OuterPosition = [0.25 0.35 12 10];  
    

end

exportgraphics(f,[exp.figPath, 'Fig3D.tiff'], 'Resolution', 600)
end