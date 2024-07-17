function [] = VT_Figure2B(exp,par)
%% Lateralisation, sorted by CP
id = 0;
figure
for sub = 1:20
    id = id +1;
    subj = exp.sub_id{sub}
    %Load full trials
    load(['D:\Elisabeth EEG Data\VolatilityTask\Data\P' subj '\S2\EEG\VT_TFA_sumData_P' subj '_2.mat'])
    
     %Load normative model to get PSI traces
    load(['C:\Users\elisa\Desktop\VolatilityTask\Analysis\ProcStim_fitted_glaze_P' subj '_S2.mat']);
    
     
    if size(dat.fitted_glaze.acc_full,1) ~= size(dat.fitted_glaze.fCPpos_full)
        isfull = ~isnan(dat.fitted_glaze.LLR(:,10));
        dat.fitted_glaze.acc_full = dat.fitted_glaze.acc_full(isfull,:);
        clear isfull;
    end

    
    if strcmp(subj, '01'); idcs = [225:235];
        dat.fitted_glaze.psi_full(idcs,:) = [];
        dat.fitted_glaze.LLR_full(idcs,:) = [];
        dat.fitted_glaze.choices_full(idcs,:) = [];
        dat.fitted_glaze.nCP_full(idcs,:) = [];
        dat.fitted_glaze.fCPpos_full(idcs,:) = [];
        dat.fitted_glaze.acc_full(idcs,:) = [];
    end
    

    if size(dat.fitted_glaze.nCP_full,1) ~= size(dat.fitted_glaze.fCPpos_full)
        isfull = ~isnan(dat.fitted_glaze.LLR(:,10)); 
        dat.fitted_glaze.nCP_full = dat.fitted_glaze.nCP(isfull,:);
        dat.fitted_glaze.fCPpos_full = dat.fitted_glaze.fCPpos(isfull,:); 
        clear isfull; 
    end

    choices = dat.fitted_glaze.choices_full;
    acc = dat.fitted_glaze.acc_full; 
    
    %Subtract baseline
    clear baselines;
    baselines = nanmean(trl_data(:,:,0.4*200:0.5*200),3);
    dbData = trl_data-baselines;
    
    %Compute lateralisation: IPSI-CONTRA
    ldata = dbData(:,exp.movCluster(1,:),:);
    rdata = dbData(:,exp.movCluster(2,:),:);
  
    leftChoice_b(id,:,:) = squeeze(nanmean(nanmean(ldata(choices == 0 & acc == 1,:,:) - rdata(choices == 0& acc == 1,:,:))));
    rightChoice_b(id,:,:) = squeeze(nanmean(nanmean(ldata(choices == 1 & acc == 1,:,:) - rdata(choices == 1& acc == 1,:,:))));
    
    topoData_left(id,:,:,:) = nanmean(dbData(choices == 0,:,:));
    topoData_right(id,:,:,:) = nanmean(dbData(choices == 1,:,:));
    
    %Lateralisation sorted by CP
    earlyCP = dat.fitted_glaze.fCPpos_full <= 3 & dat.fitted_glaze.nCP_full == 1;
    midCP = dat.fitted_glaze.nCP_full == 1 & dat.fitted_glaze.fCPpos_full > 3 & dat.fitted_glaze.fCPpos_full <= 7;
    lateCP = dat.fitted_glaze.nCP_full == 1 &  dat.fitted_glaze.fCPpos_full >= 8;

    indices = [{earlyCP}, {midCP}, {lateCP}];
    labs = {'EarlyCP', 'MidCP', 'LateCP'};
   
    assert(length(choices) == size(dat.fitted_glaze.psi_full,1));
    
    nexttile
    for cpTime = 1:3
        lat_leftChoice.(labs{cpTime})(id,:,:) = squeeze(nanmean(nanmean(ldata(choices == 0 & acc == 1& indices{cpTime},:,:) - rdata(choices == 0 & acc == 1& indices{cpTime},:,:))));
        lat_rightChoice.(labs{cpTime})(id,:,:) = squeeze(nanmean(nanmean(rdata(choices == 1 & acc == 1& indices{cpTime},:,:) - ldata(choices == 1 & acc == 1& indices{cpTime},:,:))));
        
        singleLat.(labs{cpTime})(id,:,:) = lat_leftChoice.(labs{cpTime})(id,:,:) + lat_rightChoice.(labs{cpTime})(id,:,:);
        
        normativeLeft.(labs{cpTime})(id,:) = nanmean(dat.fitted_glaze.psi_full(choices == 0 & acc == 1& indices{cpTime},2:end));
        normativeRight.(labs{cpTime})(id,:) = nanmean(dat.fitted_glaze.psi_full(choices == 1& acc == 1& indices{cpTime},2:end));
        
        if cpTime == 3; plot(normativeLeft.(labs{cpTime})(id,:), 'LineStyle', '-'); hold on;  plot(normativeRight.(labs{cpTime})(id,:), 'LineStyle', '--'); end
    end
    clear indices; clear earlyCP; clear midCP; clear lateCP;
end


%% Figure
colors = {'k', 'b', 'r'};
f = figure;
tiledlayout(3,1, 'TileSpacing', 'compact')
c = 0;
for i = 1:10
        sampleTimes(i) = [(0.1*200)+(0.4*i*200)]
end
scale = 7;
colors = hot(6);% {'k', 'r','b','r', 'g'};
DVcol = colors(3,:);
labels = {'Early CP (1-3)', 'Mid CP (4-7)', 'Late CP (8-10)'}
cpTimes = {[0.5,0.5,(0.4)*3+0.5, (0.4)*3+0.5]*200,...
    [0.4*3+0.5,0.4*3+0.5, 0.4*7+0.5,0.4*7+0.5]*200,...
    [0.4*7+0.5,0.4*7+0.5, 0.4*10+0.5,0.4*10+0.5]*200};

for cpTime = 1:3
    c = c+1;

    if par.bothResp == 0; %Plot lateralisation index (left-right resp lateralisation)
        nexttile
        
        singleSignal = lat_leftChoice.(labs{cpTime})+lat_rightChoice.(labs{cpTime});
        stdshade(singleSignal(:,timeLim(1):timeLim(2)), 0.05, 'k', '-',[],1); hold on;
        if par.plotPsi == 1
            normativeLat =  -((normativeRight.(labs{cpTime}) - normativeLeft.(labs{cpTime}))');
            plot([0,sampleTimes]+0.4*200,[0,-nanmean((normativeLat)')./scale], 'linestyle', '-','color', DVcol, 'linewidth', 1.5)
            clear normativeLat;
        end
    else %Plot Left/Right responses separately
        nexttile
        for r= 1:2
        timeLim = [0.4,(4.5+0.4)]*200;

            
            if r == 1
                stdshade(lat_leftChoice.(labs{cpTime})(:,timeLim(1):timeLim(2)), 0.05, 'k', '-',[],1); hold on;
                if par.plotPsi == 1; plot([0,sampleTimes]+0.4*200,-[0,nanmean((normativeLeft.(labs{cpTime})')')./scale], 'linestyle', '-','color', DVcol, 'linewidth', 1.5); end
                hold on
            elseif r == 2
                stdshade(-lat_rightChoice.(labs{cpTime})(:,timeLim(1):timeLim(2)), 0.05, 'k', '--',[],1); hold on;
                if par.plotPsi == 1; plot([0,sampleTimes]+0.4*200,-[0,nanmean((normativeRight.(labs{cpTime})')')./scale], 'linestyle', '--','color', DVcol, 'linewidth', 1.5); end
            end
        end
    end
  
    a = fill([cpTimes{cpTime}],[-1.5,1.5,1.5,-1.5], 'k')
    a.FaceAlpha = 0.05; a.EdgeAlpha = 0;
    xline(sampleTimes, '--k'); xlim([0,900]); xline([0.1*200]); 
   
    if c == 2,
ylabel({'\leftarrow  Error               MBL (dB)             Correct \rightarrow'})
    end
    yline(0, '--k'); 
    ylim([-1,1]); 
    set(gca, 'fontsize', 10);
    if c == 3
        xticks(([0.1,0.5:0.4:4.5]*200));
        xticklabels(([0,0.4,0.8,1.2,1.6,2,2.4,2.8,3.2,3.6,4,4.4]));
         xlabel('Time (s)');
    else
        xticks([]); xlabel([]);
        xticklabels([]);
    end
    set(gca, 'fontsize', 10);
    title(labels{cpTime}, 'fontsize', 10);
    if cpTime == 1
        legend({'','MBL', '\Psi'}, 'Location', 'se');
    end
    
end

set(gca, 'fontsize', 10);
f.Units = 'centimeters';

f.OuterPosition = [0.25 0.35 8 17];

exportgraphics(f,[exp.figPath, 'Fig2B.tiff'], 'Resolution', 300)
