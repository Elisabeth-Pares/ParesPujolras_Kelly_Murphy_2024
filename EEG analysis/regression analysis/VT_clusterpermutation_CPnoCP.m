function [cluster] = VT_clusterpermutation_CPnoCP(exp,permData, par, figName)

% Run cluster-based permutation tests
% Permutation parameters 

cfg_dat = struct;
close all; 
load('C:\Users\elisa\Desktop\VolatilityTask\Analysis\functions\chanlocsBioSemi_128_noext.mat')

chans =[1:128]; 
figure; topoplot([], chanlocs(chans))

chans = par.chans; 
cfg_dat.label       = {chanlocs.labels};  %Channel locations; (exp.CPPCluster_small)
cfg_dat.label       =  cfg_dat.label(chans);  %{'MotorLat'};Channel locations; (exp.CPPCluster_small)

cfg = struct();
cfg.feedback = 'yes';

%Triangulation method
cfg.layout           = 'biosemi128.lay';
cfg_neighb           = [];
cfg_neighb.layout    = ft_prepare_layout(cfg);%ft_prepare_layout(cfg_layout);
cfg_neighb.method    = 'triangulation'; %triangulation
cfg_neighb.feedback  = 'yes';
cfg.neighbours       = ft_prepare_neighbours(cfg_neighb); %ALL channels

cfg_dat.fsample = 200;  % sampling rate
cfg_dat.time = [0.005:1/cfg_dat.fsample:1]; 
cfg_dat.dimord = 'chan_time';

cfg.minnbchan   = [0];
cfg.latency     = par.clusterTime;%'all'; %Analysis period
cfg.channel     = 1:length(chans); %[1:128];%exp.centroParietal;  %
cfg.avgoverchan = par.averageChan;
cfg.avgovertime = 'no'; %no
cfg.parameter   = 'values';
cfg.method      = 'montecarlo';
cfg.clusterthreshold    = 'parametric'%'nonparametric_common'; %EPP added

cfg.statistic   = 'depsamplesT';  
cfg.correctm    = 'cluster';
cfg.clusteralpha = 0.05;  % intial threhsold for forming clusters
cfg.clusterstatistic = 'maxsum';  % method for quantifying combined cluster test statistics

%Note on two-tailed corrections:
%Setting cfg.alpha = 0.025 & cfg.correcttail = 'no' gives the same results
%as cfg.alpha = 0.05 & cfg.correcttail = 'alpha'. This tests each tail
%against alpha = 0.025, and p-values need to be interpreted against *that*
%critical value. 

%The alternative cfg.correcttail = 'prob' tests each tail against 0.05, and
%can be interpreted against the overall p value = 0.05.
cfg.alpha       = 0.05;
cfg.tail        = 0; % two-sided test

cfg.correcttail = 'prob'; 
%this makes the p-values of ft_freqstatistics reflect the tail of the analysis
% Without the 'prob' correction, the p-values in stats.negclusters and
% stats.posclusters would otherwise be one-tailed p-values
 
cfg.numrandomization = 10000;

%Design matrix 
Nsub = 20
cfg.design(1,1:2*Nsub)  = [ones(1,Nsub) 2*ones(1,Nsub)];
cfg.design(2,1:2*Nsub)  = [1:Nsub 1:Nsub];
cfg.ivar                = 1; % the 1st row in cfg.design contains the independent variable
cfg.uvar                = 2; % the 2nd row in cfg.design contains the subject number

%Null condition 
for s = 1:20
    all_trl_cp{s} = cfg_dat; all_trl_cp{s}.values = permData.cp{s}(chans,:);%(exp.CPPCluster_small,:);
    all_trl_nocp{s} = cfg_dat; all_trl_nocp{s}.values = permData.nocp{s}(chans,:);%(exp.CPPCluster_small,:);
    all_trl_diff{s} = cfg_dat; all_trl_diff{s}.values = all_trl_cp{s}.values-all_trl_nocp{s}.values;

    all_trl_null{s} = cfg_dat; all_trl_null{s}.values = zeros(size(all_trl_cp{1}.values,1),size(all_trl_cp{1}.values,2))

end

cluster.channels = exp.centroParietal;  
stat1 = ft_timelockstatistics(cfg,all_trl_cp{:},all_trl_nocp{:});
stat2 = ft_timelockstatistics(cfg,all_trl_diff{:},all_trl_null{:});

cluster.stat1 = stat1;
cluster.stat2 = stat2;

%% Plot results & significant channels 
% 
allStats = {'stat1'};
allEffects = {'all_trl_effEv'};
effectNames = {'|\Delta \Psi|'};

if strcmp(par.averageChan, 'no')
    for stats = 1%:4%:length(allStats)
    f = figure;
    thisStat = eval(allStats{stats});
    thisEffect = eval(allEffects{stats}); 
    
    cfg = [];
    cfg.channel   = cfg_dat.label;%'all';
    cfg.latency   = 'all';
    cfg.parameter = 'values';
    GA_Effect = ft_timelockgrandaverage(cfg, thisEffect{:});
    GA_Null = ft_timelockgrandaverage(cfg, all_trl_null{:});

	cfg = [];
    cfg.operation = 'subtract';
    cfg.parameter = 'avg';
    diffEffect = ft_math(cfg, GA_Effect, GA_Null);

    % Make a vector of all p-values associated with the clusters from ft_timelockstatistics.
    if ~isempty(thisStat.posclusters(:))
    pos_cluster_pvals = [thisStat.posclusters(:).prob];
    
    % Then, find which clusters are deemed interesting to visualize, here we use a cutoff criterion based on the
    % cluster-associated p-value, and take a 5% two-sided cutoff (i.e. 0.025 for the positive and negative clusters,
    % respectively
        pos_clust = find(pos_cluster_pvals < 0.05);
    else
        pos_clust = [];
    end
    pos       = ismember(thisStat.posclusterslabelmat, pos_clust);
    
    % and now for the negative clusters...close all 
    if ~isempty(thisStat.negclusters(:))
        neg_cluster_pvals = [thisStat.negclusters(:).prob];
        neg_clust         = find(neg_cluster_pvals < 0.05);
    else
        neg_clust = [];
    end
    neg               = ismember(thisStat.negclusterslabelmat, neg_clust);

    %Plot sequence of effect topographies
    timestep      = 0.05; % timestep between time windows for each subplot (in seconds)
    sampling_rate = 200; % Data has a temporal resolution of 300 Hz
    sample_count  = length(thisStat.time);
    % number of temporal samples in the statistics object
    j = [0:timestep:1]; % Temporal endpoints (in seconds) of the ERP average computed in each subplot
    m = [1:timestep*sampling_rate:sample_count]; % temporal endpoints in M/EEG samples
    
%     % First ensure the channels to have the same order in the average and in the statistical output.
%     % This might not be the case, because ft_math might shuffle the order
    [i1,i2] = match_str(diffEffect.label, thisStat.label);
%     
    for k = 6%1:19
        subplot(1,3,1);
%         subplot(4,5,k);

        cfg = [];
        cfg.xlim = [j(k) j(k+1)];   % time interval of the subplot
        cfg.zlim = [-2.5e-13 2.5e-13];
        % If a channel is in a to-be-plotted cluster, then
        % the element of pos_int with an index equal to that channel
        % number will be set to 1 (otherwise 0).
        
        % Next, check which channels are in the clusters over the
        % entire time interval of interest.
%         pos_int = zeros(numel(diffEffect.label),1);
%         neg_int = zeros(numel(diffEffect.label),1);
        
        pos_int(i1) = all(pos(i2, m(k):m(k+1)), 2);
        neg_int(i1) = all(neg(i2, m(k):m(k+1)), 2);
        
%         pos_int = repmat(0,128,1);
%         pos_int([1,3,19]) = 1;
%                 pos_int([1,3,19]) = 1;

%         diffEffect.avg(:) = 0;
        
        cfg.highlight   = 'on';
        % Get the index of the to-be-highlighted channel
        cfg.highlightchannel = find(pos_int | neg_int);  %exp.frontBackCluster;%
        cfg.comment     = 'xlim';
        cfg.commentpos  = 'title';
        cfg.layout      = 'biosemi128.lay';
        cfg.interactive = 'no';
        cfg.figure      = 'gca'; % plots in the current axes, here in a subplot
        cfg.parameter   = 'avg';
        cfg.markersize = 0.5;
        cfg.headrad = 0.5;
        cfg.colormap = 'jet';
       
        thisTopoTimes = nanmean(diffEffect.avg(:,j(k)*200:j(k+1)*200),2);

          VT_topoplot(thisTopoTimes,chanlocs(1:128), ...
                'electrodes', 'off', ...
                'emarker', {'.','k',[],5},...
                'emarker2', {exp.centroParietal, '*', 'k',2,0.5},...)
                'numcontour', 0, ...
                'plotrad',0.65, ...
                'headrad', 0.5)

            colorbar
%         caxis([-0.075 0.075]);
    end
    sgtitle(effectNames{stats}); 
    f.Units = 'centimeters';
    f.OuterPosition = [0.25 0.25 13 16];

    end
end