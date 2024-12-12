function [adjRsq] = VT_Figure2D(exp)
%% Plot Figure 2 panel D
clear allP
id = 0;
for sub = exp.sub_id([1:20])
    id = id+1;
    thisPath = ['P' sub{1}];
    fullPath = ([exp.dataPath, thisPath]);
    cd(fullPath)
    filename = ['VT_regression_psi_deltaPsi_TFA_rsq_HP_1_P' sub{1} '_2.mat']; 
    load(filename);
    disp(['Processing P' sub{1} '...'])

    for s = 1:9
        allP.leftMinRight_psi(id,s,:) =  squeeze(dat.prior.beta(:,:,s));
        allP.leftMinRight_effEv(id,s,:) =  squeeze(dat.deltaPsi.beta(:,:,s));
    end
    
%         adjRsq(id) = nanmean(nanmean(dat.prior.adjrsquare));

end

%% Run cluster permutation analysis
clear permData
for s = 1:20
    permData.psi{s} = squeeze(nanmean(allP.leftMinRight_psi(s,:,:),2))';
    permData.effEv{s}= squeeze(nanmean(allP.leftMinRight_effEv(s,:,:),2))';
end
cluster = VT_TFA_cluster_effEv(permData)

%%
f = figure;
hold on;
colors = {hex2rgb('#04080f'), hex2rgb('#192bc2'), hex2rgb('#78c0e0'), hex2rgb('#ade8f4')};

ypos= -0.005; 
regs = {'leftMinRight_psi', 'leftMinRight_effEv'};
labels = {'\Psi','\Delta\Psi','LLR*surprise','LLR*uncert.'};

for stats = [1:2]
    thisReg = regs{stats};
    stdshade(squeeze(nanmean(allP.(thisReg)(:,:,:),2)), 0.25, colors{stats}, '-', [], 1); hold on;
end

for stats = [1:2]
    thisReg = regs{stats};
    
    ypos = ypos - 0.0025;
    thisField = ['stat' num2str(stats)];

    for i = 1:length(cluster.(thisField).posclusters)
        if cluster.(thisField).posclusters(i).prob < 0.05
            thisMark = '*';
            plot(find((cluster.(thisField).posclusterslabelmat ==i)),ypos,'Color', colors{stats}, 'marker', thisMark, 'markersize', 3);
            hold on;
        end
    end
    plot(find((cluster.(thisField).posclusterslabelmat ==1)),ypos,'Color', colors{stats}, 'marker', thisMark, 'markersize', 3)
    hold on;
    set(gca, 'fontsize', 7)
        
end
yline([0,0], '-k')
box on 
xticks([1,100,200]); xticklabels({'0', '0.5', '1'})
xlabel('Time (s)')
ylabel([{'Regression coefficients (a.u.)'}]); 
f.Units = 'centimeters';
ylim([-0.02, 0.1])
set(gca, 'fontsize', 10)

f.OuterPosition = [0.25 0.35 8 10];
exportgraphics(f,[exp.figPath, 'Fig2D.tiff'], 'Resolution', 600)
