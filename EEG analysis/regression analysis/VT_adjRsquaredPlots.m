function [] = VT_adjRsquaredPlots(exp,par)

% For TFA
switch par.dataType
    case 'TFA'
        allP = struct();
        for s = 1:length(exp.sub_id)
            disp(['Processing sub ' num2str(s) '...'])
            thisPath = ['P' exp.sub_id{s}];
            fullPath = ([exp.dataPath, thisPath]);
            cd(fullPath)
            
            load(['VT_regression_updatedPsi_TFA_rsq_HP_1_P' exp.sub_id{s} '_2.mat'])
            allP.updatedPsi(s,:) = squeeze(nanmean(dat.updatedPsi.adjrsquare(:,:,1:10),3));
            
            load(['VT_regression_deltaPsi_TFA_rsq_HP_1_P' exp.sub_id{s} '_2.mat'])
            allP.deltaPsi(s,:) = squeeze(nanmean(dat.deltaPsi.adjrsquare(:,:,1:10),3));

            load(['VT_regression_LLR_TFA_rsq_HP_1_P' exp.sub_id{s} '_2.mat'])
            allP.LLR(s,:) = squeeze(nanmean(dat.LLR.adjrsquare(:,:,1:10),3));

        end
        
    case 'ERP'
        par.hp = 1;
        allP = struct();
        
        for s = 1:length(exp.sub_id)
            disp(['Processing sub ' num2str(s) '...'])
            
            thisPath = ['P' exp.sub_id{s}];
            fullPath = ([exp.dataPath, thisPath]);
            cd(fullPath)
   
            load(['VT_regression_updatedPsi_ERP_rsq_HP_1_P' exp.sub_id{s} '_2.mat'])
            allP.updatedPsi(s,:,:) = squeeze(nanmean(dat.updatedPsi.adjrsquare(:,:,1:10),3));
            
            load(['VT_regression_LLR_ERP_rsq_HP_1_P' exp.sub_id{s} '_2.mat'])
            allP.LLR(s,:,:) = squeeze(nanmean(dat.LLR.adjrsquare(:,:,1:10),3));
            
            load(['VT_regression_deltaPsi_ERP_rsq_HP_1_P' exp.sub_id{s} '_2.mat'])
            allP.deltaPsi(s,:,:) = squeeze((nanmean(dat.deltaPsi.adjrsquare(:,:,1:10),3)));
            allP.deltaPsi_beta(s,:,:) = squeeze((nanmean(dat.deltaPsi.beta(:,:,1:10),3)));

        end
end

%% Plot average for all channels, comparing models

f = figure;

switch par.dataType
    case 'ERP'
        testModels = {'LLR', 'deltaPsi', 'updatedPsi'}%{'zLLR', 'deltaPsi'};%{'zLLR', 'deltaPsi', 'surprise'}%, 'surprise'}%, 'zLLR', 'deltaPsi', 'surprise'}; %signed_prior', 'signed_zLLR'}%, 'signed_zLLR'}%,'signed_zLLR'};
        testModelLabs = {['deltaPsi']}%|\Delta\Psi|']}% LLR', ''}%{'LLR', ['|\Delta\Psi|'], 'surprise'};
        topo = 'deltaPsi';
        refModel = 'LLR';
        refModelLab = 'LLR';
        
        exp.occLeft = [8,9,10,11,17,16,15,14];
        exp. occRight = [30,29,28,27,37,38,39,40];
        postBroad = [17,18,19,20,21,22,30,31]
        chanSels = {postBroad}
        chanLabs = {'postBroad'}

        barLims = [0, 0.006]; barTicks = [];
        barTimes = [0.2,0.4]*200;
        ylims = [-0.002, 0.005];
        
        ytickvals = [-2,0,2,4];
        yOnset = -0.001

    case 'TFA'
        testModels = {'LLR', 'deltaPsi', 'updatedPsi'}%{'zLLR', 'deltaPsi'};%{'zLLR', 'deltaPsi', 'surprise'}%, 'surprise'}%, 'zLLR', 'deltaPsi', 'surprise'}; %signed_prior', 'signed_zLLR'}%, 'signed_zLLR'}%,'signed_zLLR'};
        topo = 'updatedPsi';
        refModel = 'LLR';
        refModelLab = 'LLR';
        
        chanLabs = {'motor'}; 
        chanSels = {1};
        
        barTimes = [1,200]; %times for bar graph
        barLims = [0.003, 0.0125]; barTicks = [0.06,0.065]
        
        ylims = [-0.007, 0.007];
        ytickvals = [-6,-3,0,3,6];
        yOnset = -0.0050;


end
colors = hot(6);

f = figure;

for cs = 1:length(chanSels)
    chans = chanSels{cs};
    handle = subplot(1, length(chanSels), cs);
    
    for model = 1:length(testModels)
        
        model1 = testModels{model};
        thisCol = colors(model,:)% colors{model}
        
        if strcmp(topo, 'deltaPsi')
            times = {{0.2*200:0.4*200}};
        elseif contains(topo, 'prior')
            times = {{1:0.2*200}};
            times = {{0.2*200:0.4*200}};
        end
        
        if cs == 1
            times = {{0.05*200:0.25*200}}; 
        else
            times = {{0.2*200:0.4*200}};
        end
       
        %Models 2 compare:
        if strcmp(par.dataType, 'ERP')
            thisTest = squeeze(nanmean(allP.(model1)(:,chans,:),2))
            thisRef = squeeze(nanmean(allP.(refModel)(:,chans,:),2))
        else
            thisTest = allP.(model1)(:,:);
            thisRef = allP.(refModel)(:,:);
        end
        
        stdshade((thisTest - thisRef), 0.05,thisCol,'-', [], 1); hold on;
        
        yticks(ytickvals/1000);
        yline([0], '-k')
        ylim(ylims)
        y = [];
        load('chanlocsBioSemi_128_noext.mat')
       ylabel(['Adj R^2 [Test - OE]']);
        
        xlabel('Time (s)')
        xticks([1,100,200]); xticklabels({'0','0.5','1'});
        set(gca, 'Fontsize', 10)
        hold on;

        if model == 1 & strcmp(par.dataType, 'ERP')
            rectangle('Position', [barTimes(1), ylims(1), barTimes(2)-barTimes(1), ylims(2)-ylims(1)], 'EdgeColor', [0 0 0 0], 'linestyle', '-', 'FaceColor',[0 0 0 0.05]); hold on
        end
        
        for s= 1:20
            permData.deltaPsi_minus_LLR.(testModels{model}).(chanLabs{cs}){s} = squeeze(thisTest(s,:) - thisRef(s,:));
            means(model,:) = squeeze(nanmean(nanmean(thisTest(:,barTimes(1):barTimes(2)))));
            stds(model,:) = squeeze(nanstd(nanmean(thisTest(:,barTimes(1):barTimes(2)))));
        end
       
    end
end

% Cluster analysis of diff in R2 differences
for cs = 1%:length(chanSels)
            ypos = yOnset;
    hold on;
    for model = 1:length(testModels)
        ypos = ypos - 0.00045;

        model1 = testModels{model};
        thisCol = colors(model,:)
        
        chanLab = chanLabs{cs};
        [cluster] = VT_clusterPerm_rSquared(exp,permData,chanLab, model1)
        hold on;
        subplot(1, length(chanSels), cs)
        
        if isfield(cluster, 'posclusters')
            for i = 1:length(cluster.posclusters)
                if cluster.posclusters(i).prob < 0.05
                    thisMark = '*';
                    
                    plot(find((cluster.posclusterslabelmat(1,:) ==i)),ypos,'Color', thisCol, 'marker', thisMark, 'markersize', 5)
                    hold on;
                end
            end
        end
        
        if isfield(cluster, 'negclusters')
            for i = 1:length(cluster.negclusters)
                if cluster.negclusters(i).prob < 0.05
                    thisMark = '*';
                    plot(find((cluster.negclusterslabelmat(1,:) ==i)),ypos,'Color', thisCol, 'marker', thisMark, 'markersize', 5)
                    hold on;
                end
            end
        end
        
    end
   
end

l=legend({'','Objective Ev (OE)','','','Effective Ev (EE)','','','Normative DV'}, 'location', 'nw')

l.ItemTokenSize(1) = 10;
l.FontSize(1) = 8;
% Inset bar graph for model comparison
% create smaller axes in top right, and plot on it
axes('Position',[.635 .46 .25 .2])

x = double(1:length(testModels))
y = double(means);
error = stds./sqrt(20);
b = bar(x,y'); hold on;
b.FaceColor = 'flat'
cols = [0 0 0; 1 0 0; 0 0 1];
b.CData= colors(1:3,:)%cols;
errorbar(x,y',error', '.k');

if strcmp(par.dataType, 'ERP'); ylabel(['Adj R^2']);
else title(['Adj R^2']); end

ylim(barLims); yticks(barTicks);
xticks([1:3]); xticklabels({'OE', 'EE', 'DV'})
xtickangle(45)
box off

if strcmp(par.dataType, 'ERP') %Plot topo of CPP effect
for cs = length(chanSels):-1:1
    hold on;

      if cs == 1
            times = {{0.05*200:0.25*200}}; 
        else
            times = {{0.2*200:0.4*200}};
      end
        
    for topoTimes = 1:length(times)
        chans = chanSels{cs};
        thisTime = barTimes(1):barTimes(2);
        axes('Position',[.62 .65 .25 .25])

        model1Topo = squeeze(nanmean(nanmean(allP.(topo)(:,:,thisTime),3)));
        model2Topo = squeeze(nanmean(nanmean(allP.(refModel)(:,:,thisTime),3)));
        nan_dat = model1Topo - model2Topo;
        
        colormap jet
        VT_topoplot(nan_dat,chanlocs(1:128), ...
            'electrodes', 'off', ...
            'emarker', {'o','k',[],5},...
            'emarker2', {chans, 'o', 'k',1,0.05},...)
            'numcontour', 0, ...
            'plotrad',0.65,...
            'headrad', 0.5)
        hold on;
        
%         caxis([-ylims(2)+0.05 ylims(2)])%300/100000]);
    end
    
end
end

f.Units = 'centimeters';
f.OuterPosition = [0.25 0.35 8 10];


