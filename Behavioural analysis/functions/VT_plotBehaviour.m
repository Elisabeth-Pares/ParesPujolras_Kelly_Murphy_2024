function [] = VT_plotBehaviour(exp, par,GA)

GA.data.nCP = GA.data.nCP*100;
%% Plot estimated pars for all participants 

theseFields = {'H' 'gainB', 'noise'}%fields(GA.fitted_glaze);
parLab = {'H', 'Gain', 'Noise'};
f = figure; 
cmap = jet(20)
for ff = 1:length(theseFields)
    thisField = theseFields{ff};
subplot(1,3,ff)

if ff == 1
[~,idx]= sort(GA.fitted_glaze.(thisField)(:)); 
    scatter([2:21],sort(GA.fitted_glaze.(thisField)(:)),10, cmap, 'filled')
else
    scatter([2:21],(GA.fitted_glaze.(thisField)(idx)),10, cmap, 'filled')
end

if ff == 2
    title('Parameter estimates')
    xlabel('Participant')
end

box on 

    hold on; 
     b = bar(nanmean(GA.fitted_glaze.(thisField)), 'FaceColor', [0.9,0.9,0.9]);%[0.9,0.9,0.9]); ga_fittedPerf;
    hold on;
    errorbar(1,nanmean(GA.fitted_glaze.(thisField)(:)), std(GA.fitted_glaze.(thisField)(:)/sqrt(20)), 'Color', 'k') %Last sample

    if ff == 1; yline([0.1], '-k'); end
    if ff == 1;    ylabel(['$$\hat{' parLab{ff} '}$$'], 'Interpreter','Latex')
    elseif ff == 2; ylabel(['Evidence ' parLab{ff}]); 
    else ylabel(parLab{ff}); end
    
    ylims = get(gca, 'ylim')
    ylim(ylims)
    set(gca, 'xTickLabels', {}, 'fontsize', 10)
    xlim([0,21]);
    set(gca, 'fontsize', 10)
    f.Units = 'centimeters';
    f.OuterPosition = [0.25 0.35 20 8];
end

% exportgraphics(f,[exp.figPath, 'FigHestimates.tiff'], 'Resolution', 600)

% figure; 
% hist(GA.fitted_glaze.H(:,1))
%% Bar graph with %cp trials 
f = figure;
b = bar(nanmean(GA.data.nCP), 'FaceColor', [0.9,0.9,0.9]);%[0.9,0.9,0.9]); ga_fittedPerf;
% colors =  {'#5A189A', '#177E89','#FFC857', '#000000',[0.5 0.5 0.5]};

% ylim([0.7,0.9]);
b.FaceColor = 'flat';
% b.CData(:,:) = colors;
hold on;
%Add error bars to averaged models
% errorbar(1,ga_fittedPerf, std(GA.fitted_perf.(par.thisAcc)), 'Color', 'k') %Perfect acc
errorbar(1,nanmean(GA.data.nCP(:,1)), std(GA.data.nCP(:,1)), 'Color', 'k') %Last sample
errorbar(2,nanmean(GA.data.nCP(:,2)), std(GA.data.nCP(:,2)), 'Color', 'k') %Perfect acc
errorbar(3,nanmean(GA.data.nCP(:,3)), std(GA.data.nCP(:,3)), 'Color', 'k') %Perfect acc

set(gca, 'xTickLabels', {'0', '1', '1<'}, 'fontsize', 9)
xlabel('# CP')
ylabel('% Trials')

f.Units = 'centimeters';
f.OuterPosition = [0.25 0.35 4 6];
exportgraphics(f,[exp.figPath, 'FigPSwitch.tiff'], 'Resolution', 600)

%% Bar graph reproducing Murphy et al. Fig2a
glazeColor = hex2rgb('#0466C8'); 
idealColor = [0.5 0.5 0.5]; 
naboundsColor = hex2rgb('#007F5F');
lastColor = hex2rgb('#DB3A34'); 
leakColor = hex2rgb('#FFC857'); 

ga_Ideal = mean(GA.ideal.(par.thisAcc));
ga_fittedGlaze= mean(GA.fitted_glaze.(par.thisAcc));
% ga_fittedPerf = mean(GA.fitted_perf.(par.thisAcc));
ga_fittedPerf_naBounds= mean(GA.fitted_perf_naBounds.(par.thisAcc));
ga_lastSample = mean(GA.last_sample.(par.thisAcc));
ga_fittedLeak = mean(GA.fitted_leak.(par.thisAcc));

%Plot data & overlay normative fit model
[orderedData, order] = sort(GA.data.acc);

f = figure;
b = bar([ga_lastSample;ga_fittedLeak;ga_fittedGlaze; ga_fittedPerf_naBounds; ga_Ideal;  orderedData; ], 'FaceColor', [0.9,0.9,0.9]);%[0.9,0.9,0.9]); ga_fittedPerf;
% colors = [[1 0 0]; [1 1 0]; [0 0 1]; [0 1 0]; [0 0 0]; repmat([1,1,1],20,1)];
% colors = [ 1,0.8,1; 1,0.8,0.8; 1,1,0.8;0.8,1,0.8;0.8,0.8,1;0.8,0.8,0.8;;repmat([1,1,1],20,1);];
% colors = [1,0.8,0.8; 1,1,0.8;0.8,0.8,1;0.8,1,0.8;0.8,0.8,0.8;;repmat([1,1,1],20,1);];
colors =  [lastColor; leakColor; glazeColor; naboundsColor; idealColor; repmat([1,1,1],20,1);];
% par.colors = {'#0466C8', '#DB3A34','#FFC857', '#000000',[0.5 0.5 0.5]}; %Blau, 'vermell', 'groc', negre %Glaze, 'na bounds', 'leak', 'data'

ylim([0.7,0.9]);
b.FaceColor = 'flat';
b.CData(:,:) = colors;
hold on;

%Add error bars to averaged models
% errorbar(1,ga_fittedPerf, std(GA.fitted_perf.(par.thisAcc)), 'Color', 'k') %Perfect acc
errorbar(1,ga_lastSample, std(GA.last_sample.(par.thisAcc))/sqrt(size(GA.last_sample.(par.thisAcc),1)), 'Color', 'k') %Last sample
errorbar(2,ga_fittedLeak, std(GA.fitted_leak.(par.thisAcc))/sqrt(size(GA.last_sample.(par.thisAcc),1)), 'Color', 'k') %Perfect acc
errorbar(3,ga_fittedGlaze, std(GA.fitted_glaze.(par.thisAcc))/sqrt(size(GA.last_sample.(par.thisAcc),1)), 'Color', 'k') %Perfect acc
errorbar(4,ga_fittedPerf_naBounds, std(GA.fitted_perf_naBounds.(par.thisAcc))/sqrt(size(GA.last_sample.(par.thisAcc),1)), 'Color', 'k') %Perfect acc
errorbar(5,ga_Ideal,  std(GA.ideal.acc)/sqrt(size(GA.last_sample.(par.thisAcc),1)), 'Color', 'k') %Ideal acc

%Add normative fits for each subject
plot([repmat(NaN, 5,1); GA.fitted_glaze.(par.thisAcc)(order,:)],'.','MarkerEdgeColor',glazeColor, 'MarkerSize', 15)

%Add perf acc , NA bounds for each subject
plot([repmat(NaN, 5,1);GA.fitted_perf_naBounds.(par.thisAcc)(order,:)], '.','MarkerEdgeColor', naboundsColor, 'MarkerSize', 15)

% %Add perf acc , NA bounds for each subject
% plot([NaN; NaN; GA.fitted_perf.acc(order,:)], '.m', 'MarkerSize', 25)

%Add perf acc , NA bounds for each subject
plot([repmat(NaN, 5,1);GA.fitted_leak.(par.thisAcc)(order,:)], '.','MarkerEdgeColor', leakColor, 'MarkerSize', 15)

set(gca, 'xTickLabels', {}, 'fontsize', 10)
xlabel('Participant data or model')
ylabel('Proportion correct')
title('Model fits')
f.Units = 'centimeters';
f.OuterPosition = [0.25 0.35 12 8];
exportgraphics(f,[exp.figPath, 'FigS1A.tiff'], 'Resolution', 600)

%% Bar graph reproducing Murphy et al. Fig2a
glazeColor = hex2rgb('#0466C8'); 
idealColor = [0.5 0.5 0.5]; 
perfColor = hex2rgb('#FF00FF');
naboundsColor = hex2rgb('#007F5F');
lastColor = hex2rgb('#DB3A34'); 
leakColor = hex2rgb('#FFC857'); 

ga_Ideal = mean(GA.ideal.(par.thisAcc));
ga_fittedGlaze= mean(GA.fitted_glaze.(par.thisAcc));
ga_perf = mean(GA.perf.(par.thisAcc));
ga_lastSample = mean(GA.last_sample.(par.thisAcc));

%Plot data & overlay normative fit model
[orderedData, order] = sort(GA.data.acc);

f = figure;
b = bar([ga_lastSample;ga_perf;ga_fittedGlaze; ga_Ideal;  orderedData; ], 'FaceColor', [0.9,0.9,0.9]);%[0.9,0.9,0.9]); ga_fittedPerf;
colors =  [lastColor; perfColor; glazeColor; idealColor; repmat([1,1,1],20,1);];
% par.colors = {'#0466C8', '#DB3A34','#FFC857', '#000000',[0.5 0.5 0.5]}; %Blau, 'vermell', 'groc', negre %Glaze, 'na bounds', 'leak', 'data'

ylim([0.7,0.9]);
b.FaceColor = 'flat';
b.CData(:,:) = colors;
hold on;

%Add error bars to averaged models
% errorbar(1,ga_fittedPerf, std(GA.fitted_perf.(par.thisAcc)), 'Color', 'k') %Perfect acc
errorbar(1,ga_lastSample, std(GA.last_sample.(par.thisAcc)/sqrt(size(GA.last_sample.(par.thisAcc),1))), 'Color', 'k') %Last sample
% errorbar(2,ga_fittedLeak, std(GA.fitted_leak.(par.thisAcc)), 'Color', 'k') %Perfect acc
errorbar(3,ga_fittedGlaze, std(GA.fitted_glaze.(par.thisAcc)/sqrt(size(GA.last_sample.(par.thisAcc),1))), 'Color', 'k') %Perfect acc
errorbar(2,ga_perf, std(GA.perf.(par.thisAcc)/sqrt(size(GA.perf.(par.thisAcc),1))), 'Color', 'k') %Perfect acc
errorbar(4,ga_Ideal, std(GA.ideal.acc)/sqrt(size(GA.last_sample.(par.thisAcc),1)), 'Color', 'k') %Ideal acc

%Add normative fits for each subject
plot([repmat(NaN,4 ,1); GA.fitted_glaze.(par.thisAcc)(order,:)],'.','MarkerEdgeColor',glazeColor, 'MarkerSize', 15)

% %Add perf acc , NA bounds for each subject
% plot([NaN; NaN; GA.fitted_perf.acc(order,:)], '.','MarkerEdgeColor','m', 'MarkerSize', 15)

set(gca, 'xTickLabels', {}, 'fontsize', 11)
xlabel('Participant data or model')
ylabel('Proportion correct')
title('Model fits')
f.Units = 'centimeters';
f.OuterPosition = [0.25 0.35 10 8];
exportgraphics(f,[exp.figPath, 'Fig1A.tiff'], 'Resolution', 600)


%% Plot GA subject beta scores (EPP)
% Line graph reproducing Murphy et al. Fig2b

par.models2plot= {'fitted_glaze', 'fitted_perf_naBounds','fitted_leak','data' };% {'fitted_perf', 'data', 'fitted_glaze', 'fitted_perf'}%,'norm_glaze', 'fitted_glaze'};
par.colors = {glazeColor, naboundsColor, leakColor, 'k'}

if strcmp(par.regression, 'Murphy')
    par.params = {'B_LLR', 'B_surprise', 'B_uncertainty'};%, 'B_CP', 'B_finalDist'};
    par.labels = {'LLR','LLR*surprise','LLR*uncertainty'};
elseif strcmp(par.regression, 'evDiscount') || strcmp(par.regression, 'priorDiscount')
    par.params = {'B_LLR', 'B_LLRxdiscount'};%, 'B_CP', 'B_finalDist'};
    par.labels = {'LLR','LLR*Scaling'};
end

fig = figure;
for f = 1:length(par.params)
    
    subplot(1,length(par.params),f)
    
    for model = 1:length(par.models2plot)
        par.thisModel = par.models2plot{model}
        
        %hold on;
        thisField = par.params{f};
        
        %Replace zeros by NaNs
        GA.(par.thisModel).(thisField)(GA.(par.thisModel).(thisField) == 0) = NaN;
        
        stdshade(GA.(par.thisModel).(thisField),0.15,par.colors{model}, '-', [], 1);
        if f > 1 
            xlim([-1,9]); xticks([0,4,9]); xticklabels({'0','5','10'});
        end
        hold on
        if strcmp(par.thisModel, 'data')
            nperm=10000; clustalpha=0.05; alpha=0.05;
            [sig_ts1,sig_ts_uncorr1,cP1,~] = cluster_permWS_fast(cat(3,GA.data.(thisField), zeros(size(GA.data.(thisField)))),nperm,clustalpha,alpha);
            if f == 1
            [sig_ts2,~,cP2,~] = cluster_permWS_fast(cat(3,nanmean(GA.data.(thisField)(:,1:5),2), nanmean(GA.data.(thisField)(:,6:10),2)),nperm,clustalpha,alpha);
            end
            if f == 1; plot(sig_ts1*-0.35, '*k', 'linewidth', 1)
            else
                plot(sig_ts1*-0.9, '*k', 'linewidth', 1); %-0.15
            end
        end
    end
    yline([0 0],'--k'); xlim([0,11])
    if f == 1; ylabel('Weight on choice (a.u.)');ylim([-0.6,4.5]); end
    if f == 2; xlabel('Sample position');  ylim([-1,1]); end;
    if f == 3; ylim([-1,1]); end;
  
    title(par.labels{f});
    set(gca, 'fontsize', 10)
    
end


%Descriptives
%Observerd
mean(GA.data.acc)*100
nanstd(GA.data.acc/sqrt(20))*100
%Last sample
mean(GA.last_sample.acc)*100
nanstd(GA.last_sample.acc/sqrt(20))*100
%Ideal
mean(GA.ideal.acc)*100
nanstd(GA.ideal.acc/sqrt(20))*100

[H, P] = ttest(GA.data.acc, GA.last_sample.acc);
[H, P] = ttest(GA.data.acc, GA.ideal.acc);

%Compare data betas
early = nanmean(GA.data.B_LLR(:,1:5),2);
late = nanmean(GA.data.B_LLR(:,6:end),2);
[H, P] = ttest(early, late);

fig.Units = 'centimeters';
fig.OuterPosition = [0.25 0.35 20 8];
exportgraphics(fig,[exp.figPath, 'FigS1B_' par.regression '.tiff'], 'Resolution', 600)

%% Plot GA subject beta scores (EPP)
% Line graph reproducing Murphy et al. Fig2b

par.models2plot= {'fitted_glaze', 'data' };% {'fitted_perf', 'data', 'fitted_glaze', 'fitted_perf'}%,'norm_glaze', 'fitted_glaze'};
par.colors = {glazeColor, 'k'}

if strcmp(par.regression, 'Murphy')
    par.params = {'B_LLR', 'B_surprise', 'B_uncertainty'};%, 'B_CP', 'B_finalDist'};
    par.labels = {'LLR','LLR*surprise','LLR*uncertainty'};
end

fig = figure;
for f = 1:length(par.params)
    
    subplot(1,length(par.params),f)
    
    for model = 1:length(par.models2plot)
        par.thisModel = par.models2plot{model}
        
        %hold on;
        thisField = par.params{f};
        
        %Replace zeros by NaNs
        GA.(par.thisModel).(thisField)(GA.(par.thisModel).(thisField) == 0) = NaN;
        
        stdshade(GA.(par.thisModel).(thisField),0.15,par.colors{model}, '-', [], 1);
        if f > 1 
            xlim([-1,9]); xticks([0,4,9]); xticklabels({'0','5','10'});
        end
        hold on
        if strcmp(par.thisModel, 'data')
            nperm=10000; clustalpha=0.05; alpha=0.05;
            [sig_ts1,sig_ts_uncorr1,cP1,~] = cluster_permWS_fast(cat(3,GA.data.(thisField), zeros(size(GA.data.(thisField)))),nperm,clustalpha,alpha);
            if f == 1
            [sig_ts2,sig_ts_uncorr2,cP2,~] = cluster_permWS_fast(cat(3,nanmean(GA.data.(thisField)(:,1:5),2), nanmean(GA.data.(thisField)(:,6:10),2)),nperm,clustalpha,alpha);
            end
            if f == 1; plot(sig_ts1*-0.35, '*k', 'linewidth', 1)
            else
                plot(sig_ts1*-0.9, '*k', 'linewidth', 1); %-0.15
            end
        end
    end
    yline([0 0],'--k');
    if f == 1; ylabel('Weight on choice (a.u.)');ylim([-0.6,4.5]); end
    if f == 2; xlabel('Sample position');  ylim([-1,1]); end;
    if f == 3; ylim([-1,1]); end;
  
    title(par.labels{f});
    set(gca, 'fontsize', 11)
%     if f == 1; legend({'Normative', 'Perfect', 'Perf. NA Bounds',  'Leaky', 'Data'}); end
    %     if f == 3; legend(par.models2plot); end
    
end


%Descriptives
%Observerd
mean(GA.data.acc)*100
nanstd(GA.data.acc/sqrt(20))*100
%Last sample
mean(GA.last_sample.acc)*100
nanstd(GA.last_sample.acc/sqrt(20))*100
%Ideal
mean(GA.ideal.acc)*100
nanstd(GA.ideal.acc/sqrt(20))*100

[H, P] = ttest(GA.data.acc, GA.last_sample.acc);
[H, P] = ttest(GA.data.acc, GA.ideal.acc);

%Compare data betas
early = nanmean(GA.data.B_LLR(:,1:5),2);
late = nanmean(GA.data.B_LLR(:,6:end),2);
[H, P] = ttest(early, late);

fig.Units = 'centimeters';
if strcmp(par.regression, 'evDiscount')
fig.OuterPosition = [0.25 0.35 13 8];
else
fig.OuterPosition = [0.25 0.35 18 8];
end
exportgraphics(fig,[exp.figPath, 'Fig1B_' par.regression '.tiff'], 'Resolution', 600)


%% Plot P(correct) as a function of time of last CP
fig = figure
par.models2plot= flip({'fitted_glaze', 'fitted_perf_naBounds','fitted_leak', 'last_sample', 'data'});% {, 'fitted_glaze', 'fitted_perf'}%,'norm_glaze', 'fitted_glaze'};
% par.colors = {'k','b', 'm',[0.2,1,0.2],'y',[1,0.2,0.2]};
par.colors = flip({glazeColor, naboundsColor, leakColor, lastColor, [0 0 0]});
%BLUE, RED, YELLOW, BLACK, GREY 

for m = 1:length(par.models2plot)
    thisModel = par.models2plot{m}
%     x = [1:5];
%     y = nanmean(GA.(thisModel).acc_fsmp);
    %     plot(x,y), 'Color', par.colors{m})
%     error = nanstd(GA.(par.thisModel).acc_fsmp)/sqrt(20);
%     errorbar(x,y,error, 'Color', par.colors{m}, 'LineWidth', 1.5); hold on;
    stdshade(GA.(thisModel).acc_fsmp, 0.15, par.colors{m}, '-', [], 1); hold on
    hold on;
end
xlabel({'Final state duration (samples)'}) %Samples from';'final change point'})
ylabel('Proportion correct')
set(gca,'fontsize', 10,'xTick', [1:5],'xTickLabels', {'1-2', '3-4', '5-6', '7-8', '9-10'})
yline([0.5,0.5], '--k')
% legend({'Data', 'Fitted Glaze', 'Fitted Perfect NA Bounds', 'Fitted leaky','Last sample'})


fig.Units = 'centimeters';
fig.OuterPosition = [0.25 0.35 7 8];
title('Temporal integration')
exportgraphics(fig,[exp.figPath, 'Fig1D.tiff'], 'Resolution', 600)
%% Plot model BICs -scatter plot

figure
subplot(1,3,1)
par.models2plot= {'fitted_perf', 'fitted_leak', 'fitted_perf_naBounds','fitted_glaze'};% {'data', 'fitted_glaze', 'fitted_perf'}%,'norm_glaze', 'fitted_glaze'};
% par.colors = {'m','y',[0.2,1,0.2],'b' };
par.colors = {[1 0 0], leakColor, naboundsColor, glazeColor}

for m = 1:length(par.models2plot)
    thisModel = par.models2plot{m}
    scatter(m, GA.(thisModel).BIC,10, [0.5,0.5,0.5], 'jitter', 'on', 'jitterAmount', 0.05);
    line([m-0.25, m+0.25], [mean(GA.(thisModel).BIC), mean(GA.(thisModel).BIC)], 'Color', par.colors{m}, 'LineWidth', 3);
    hold on
    allBics(m) =  mean(GA.(thisModel).BIC);
end
xlim([0,5]);
ylabel('BIC');
xticks([1,2,3,4])
xticklabels({'Perf', 'Leaky', 'Perf. NA bounds','Normative'})


bestBic = min(allBics);

%%
f = figure
subplot(1,3,2)

for m = 1:length(par.models2plot);
    thisModel = par.models2plot{m}
    scatter(m, GA.(thisModel).BIC-GA.fitted_glaze.BIC,10, [0.5,0.5,0.5], 'jitter', 'on', 'jitterAmount', 0.05);
    line([m-0.25, m+0.25], [mean(GA.(thisModel).BIC), mean(GA.(thisModel).BIC)]-bestBic, 'Color', par.colors{m}, 'LineWidth', 3);
    hold on
    yline(0, '-k')
    %     ylim([-40,1200])
    
end
ylim([-40,1200])
xlim([0,5]);
xticks([1,2])
xticklabels({'Perf.'})
ylabel('\DeltaBIC')
xticks([2,3,4])
xticklabels({'Perf.', 'Leaky', 'Perf. NA Bounds','Normative'})

subplot(1,3,3)
for m = 1:length(par.models2plot)
    thisModel = par.models2plot{m}
    scatter(m, GA.(thisModel).BIC-GA.fitted_glaze.BIC,10, [0.5,0.5,0.5], 'jitter', 'on', 'jitterAmount', 0.05);
    line([m-0.25, m+0.25], [mean(GA.(thisModel).BIC), mean(GA.(thisModel).BIC)]-bestBic, 'Color', par.colors{m}, 'LineWidth', 3);
    hold on
    yline(0, '-k')
    ylim([-40,160])
end
xlim([0,5]);
ylabel('\DeltaBIC')
xticks([2,3,4])
xticklabels({'Leaky', 'Perf. NA Bounds','Normative'})
ylim([-50,90]);
% Stats
[H1, P1] = ttest(GA.fitted_glaze.BIC, GA.fitted_perf.BIC);
[H2, P2] =ttest(GA.fitted_glaze.BIC, GA.fitted_leak.BIC)
[H3, P3] =ttest(GA.fitted_glaze.BIC, GA.fitted_perf_naBounds.BIC)

allP = [P1,P2,P3]*3;

%%
par.models2plot= {'fitted_leak', 'fitted_perf_naBounds','fitted_glaze'};% {'data', 'fitted_glaze', 'fitted_perf'}%,'norm_glaze', 'fitted_glaze'};
% par.colors = {'y',[0.2,1,0.2],'b' };
par.colors = {leakColor, naboundsColor, glazeColor}


f = figure; 
for m = 1:length(par.models2plot)
    thisModel = par.models2plot{m}
    scatter(m, GA.(thisModel).BIC-GA.fitted_glaze.BIC,10, [0.5,0.5,0.5], 'jitter', 'on', 'jitterAmount', 0.05);
    line([m-0.25, m+0.25], [mean(GA.(thisModel).BIC), mean(GA.(thisModel).BIC)]-bestBic, 'Color', par.colors{m}, 'LineWidth', 3);
    hold on
    yline(0, '-k')
    ylim([-40,160])
end
xlim([0,4]);
ylabel('\DeltaBIC')
xticks([1,2,3])
title('Model comparison')
xlabel('Best fitting models')
% xticklabels({'Leaky', 'NA bounds','Normative'})
xticklabels({'', '',''})

ylim([-50,90]);
% Stats
[H1, P1] = ttest(GA.fitted_glaze.BIC, GA.fitted_perf.BIC);
[H2, P2] =ttest(GA.fitted_glaze.BIC, GA.fitted_leak.BIC)
[H3, P3] =ttest(GA.fitted_glaze.BIC, GA.fitted_perf_naBounds.BIC)

set(gca, 'fontsize', 10)
allP = [P1,P2,P3]*3;
f.Units = 'centimeters';
f.OuterPosition = [0.25 0.35 6 8];
exportgraphics(f,[exp.figPath, 'Fig1C.tiff'], 'Resolution', 600)

