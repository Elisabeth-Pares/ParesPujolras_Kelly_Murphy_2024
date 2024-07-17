function [adjRsq] = VT_figure2C(exp)
%% Plot Figure 2 panel C
%% Plot GA

id = 0;
for sub = exp.sub_id([1:20])
    id = id+1;
    thisPath = ['P' sub{1}];
    fullPath = ([exp.dataPath, thisPath]);
    cd(fullPath)
        filename = ['VT_regression_standard_noUnc_TFA_allT_subtr_lateralisationrsq1_P' sub{1} '_2_LOG.mat'];

    load(filename);
    disp(['Processing P' sub{1} '...'])
    
    for s = 1:9
        allP.leftMinRight(id,s,:) =  squeeze(dat.prior.beta(:,:,s));
        allP.leftMinRight_llr(id,s,:) =  squeeze(dat.zLLR.beta(:,:,s));
        allP.leftMinRight_llrxsurp(id,s,:) =  squeeze(dat.zLLRxSurprise.beta(:,:,s));
    end
    
    adjRsq(id) = squeeze(nanmean(nanmean(dat.prior.adjrsquare)));
end

%% Do on raw beta band power
clear permData
for s = 1:20
    permData.prior{s} = squeeze(nanmean(allP.leftMinRight(s,:,:),2))';
    permData.LLR{s}= squeeze(nanmean(allP.leftMinRight_llr(s,:,:),2))';
    permData.LLRxSurprise{s} = squeeze(nanmean(allP.leftMinRight_llrxsurp(s,:,:),2))';
end
cluster = VT_TFA_cluster(permData)

%% Plot panel C
f = figure;
colors = {hex2rgb('#04080f'), hex2rgb('#192bc2'), hex2rgb('#0096c7'), hex2rgb('#ade8f4')};

stdshade(squeeze(nanmean(allP.leftMinRight(:,:,:),2)),0.25, colors{1}, '-', [], 1); hold on;
stdshade(squeeze(nanmean(allP.leftMinRight_llr(:,:,:),2)),0.25, colors{2}, '-', [], 1);hold on;
stdshade(squeeze(nanmean(allP.leftMinRight_llrxsurp(:,:,:),2)),0.25, colors{3}, '-', [], 1); hold on;

yline(0); 
xticks([1,100,200]); xticklabels({'0','0.5','1'});
xlabel('Time(s)'); ylabel('Regression coefficients (a.u.)')
ylim([-0.02,0.1])
set(gca, 'fontsize', 10) 

hold on;
ypos= -0.004
for stats = [1:3]
    ypos = ypos - 0.0025;
    thisField = ['stat' num2str(stats)];
    thisCol = colors{stats};
    for i = 1:length(cluster.(thisField).posclusters)
        if cluster.(thisField).posclusters(i).prob < 0.05
            thisMark = '*';
            plot(find((cluster.(thisField).posclusterslabelmat ==i)),ypos,'Color', colors{stats}, 'marker', thisMark, 'markersize', 3);
            hold on;
        end
    end
    
    hold on;
    set(gca, 'fontsize', 7) 
end

f.Units = 'centimeters';
f.OuterPosition = [0.25 0.35 8 10];

exportgraphics(f,[exp.figPath, 'Fig2C.tiff'], 'Resolution', 600)
