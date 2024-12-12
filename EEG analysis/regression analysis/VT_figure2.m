function [] = VT_figure2(exp)
%% Figure 2 panels 
%% A - Plot lateralisation from trial onset 
VT_Figure2A(exp)
%% B - lateralisation sorted by CP time 
par.plotPsi = 1 ; %Overlay estimated Psi for each CP
par.bothResp = 1; %Plot L/R response separately (1) or as an index (0)
VT_Figure2B(exp,par) %Plot lateralisation sorted by CP time 
%% C & D - regression results
adjR_C = VT_Figure2C(exp) %Full regression: Psi + LLR + LLR*surprise 
adjR_D = VT_Figure2D(exp) %Brief regression: Psi + deltaPsi 
%Compare fit
[~,p,~,stats] = ttest(adjR_C, adjR_D); 

%% Figure 2C, quartile split MBL 
par.nGroups = 4; %How many groups to split by?
par.dataType = 'TFA';
VT_figure3D(exp, par)
%% Plot LEFT/RIGHT HEMI data - Figure S2
VT_FigureS2 %Left/Right hemispheres regression results 
