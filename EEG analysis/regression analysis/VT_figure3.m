function [] = VT_figure3(exp)

%% Figure 3A/ Figure S3
par.hp = 1; %HP filter >1Hz for Fig. S3
par.lp = 0; %LP filter <6Hz for Viz. 
VT_figure3A(exp,par)

%% Figure 3B - CPP ~ LLR + LLR*surprise
par.chans = exp.centroParietal;%Which channels to include in cluster analysis?
par.averageChan = 'yes'; %For statistics, average over centroparietal cluster
par.clusterTime = [0 1]; %Time window to test

VT_figure3B(exp,par)
%% Figure 3C - CPP ~ deltaPsi
par.chans = exp.centroParietal; %Which channels to include in cluster analysis?
par.averageChan = 'yes'; %For statistics, average over centroparietal cluster
par.clusterTime = [0 1]; %Time window to test
VT_figure3C(exp,par)

%% Figure 3D, median split evoked potentials of equal |LLR| by high/low deltaPsi
par.nGroups = 4; %How many groups to split by?
par.dataType = 'ERP';
VT_figure3D(exp, par)
