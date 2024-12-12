% Run cluster-based permutation tests
% Permutation parameters 
function [stat1] = VT_clusterPerm_rSquared(exp,permData, chanLab, modelLab)
cfg_dat = struct;
% close all; 

load('C:\Users\elisa\Desktop\VolatilityTask\Analysis\functions\chanlocsBioSemi_128_noext.mat')

% load('chanlocsBioSemi_128_noext.mat') %Load channel locations
chans = 1:128;
% figure; topoplot([], chanlocs(chans)); 

chans = 1; %Just lateralisation

% load('chanlocsBioSemi_128_noext.mat') %Load channel locations
cfg_dat.label       =  {'RsqDiff'};%Channel locations; (exp.CPPCluster_small)

cfg = struct();
% cfg.feedback = 'yes';

%Triangulation method
cfg.layout           = 'biosemi128.lay';
cfg_neighb           = [];
cfg_neighb.layout    = ft_prepare_layout(cfg);%ft_prepare_layout(cfg_layout);
cfg_neighb.method    = 'triangulation'; %triangulation
cfg_neighb.feedback  = 'no';
cfg.neighbours       = ft_prepare_neighbours(cfg_neighb); %ALL channels


cfg_dat.fsample = 200;  % sampling rate
cfg_dat.time = [0.005:1/cfg_dat.fsample:1]; 
cfg_dat.dimord = 'chan_time';

cfg.minnbchan   = [0];
cfg.latency     = 'all'; %'all'; %Analysis period
cfg.channel     = 1;  %
cfg.avgoverchan = 'no';
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
 
cfg.numrandomization = 1000;

%Design matrix 
Nsub = 20
cfg.design(1,1:2*Nsub)  = [ones(1,Nsub) 2*ones(1,Nsub)];
cfg.design(2,1:2*Nsub)  = [1:Nsub 1:Nsub];
cfg.ivar                = 1; % the 1st row in cfg.design contains the independent variable
cfg.uvar                = 2; % the 2nd row in cfg.design contains the subject number

%Null condition 
for s = 1:20
    all_trl_diff{s} = cfg_dat; all_trl_diff{s}.values = permData.deltaPsi_minus_LLR.(modelLab).(chanLab){s}(1,:);%(exp.CPPCluster_small,:);

    all_trl_null{s} = cfg_dat;
    all_trl_null{s}.values = zeros(size(all_trl_diff{1}.values,1),size(all_trl_diff{1}.values,2))

end

stat1 = ft_timelockstatistics(cfg,all_trl_diff{:},all_trl_null{:});

stat1.channels = exp.centroParietal;

