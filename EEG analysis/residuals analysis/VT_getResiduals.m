function [res] = VT_getResiduals(nsubj, exp, par)

%% Load regression residuals
allsubj = exp.sub_id; % ORDER IS IMPORTANT FOR TRIAL-IDs!!!!
sub = allsubj{nsubj};

subj = ['P' allsubj{nsubj}];

sess = 2;
thisPath = ['P' sub '\'];
fullPath = ([exp.dataPath, thisPath]);

cd(fullPath);

load([fullPath, 'VT_regression_deltaPsi_ERP_rsq_HP_1_' subj '_2.mat'])

thisPar = 'deltaPsi'; maxS = 10;

%Divide whole scalp in various regions 
frontalCluster = [89,85,76,...
    88,86,75,...
    98,87,66]

posteriorCluster = [4,19,20,21,22,23]; 

frontBackCluster = [frontalCluster, posteriorCluster];

load('C:\Users\elisa\Desktop\VolatilityTask\Analysis\functions\chanlocsBioSemi_128_noext.mat')
%%
figure
nan_dat = repmat(0,1,128)

for c = 1:2
    if c == 1; channels = frontalCluster; color = 'k'; 
    elseif c == 2; channels = posteriorCluster; color = 'b'; 
    end
    
    hold on; VT_topoplot(nan_dat,chanlocs(1:128), ...
        'style', 'contour',...
        'electrodes', 'on', ...
        'emarker', {'.','k',[],5},...
        'emarker2', {channels, '*', color,5,2},...)
        'numcontour', 0, ...
        'plotrad',0.65)
end

%% Get timeResolved residuals
for s = 1:maxS
    % Residuals by region, for the whole epoch time
    res.timeResolved_frontalCluster(s,:,:) = squeeze(nanmean(dat.(thisPar).r(frontalCluster,:,s,:)));
    res.timeResolved_posteriorCluster(s,:,:) = squeeze(nanmean(dat.(thisPar).r(posteriorCluster,:,s,:)));
    res.timeResolved_frontPostCluster(s,:,:) = squeeze(nanmean(dat.(thisPar).r(frontBackCluster,:,s,:)));
end
save([fullPath, '\VT_CPPresiduals_' par.reg '.mat'], 'res');

end