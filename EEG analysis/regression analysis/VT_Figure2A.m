function [] = VT_Figure2A(exp)
%% Load for all participants and save in the same file 
id = 0; clear leftChoice; clear rightChoice; clear allLeft; clear allRight; 

for sub = 1:20
    id = id +1;
    subj = exp.sub_id{sub}
    %Load full trials
    load(['D:\Elisabeth EEG Data\VolatilityTask\Data\P' subj '\S2\EEG\VT_TFA_sumData_P' subj '_2.mat'])

    %Load and concatenate behaviour
    load([exp.procBehPath, 'PreprocData_' subj '-.mat']) 
    
    choices = allBeh.allChoices(logical(allBeh.fullTrials));
  
   
    clear baselines; 
    baselines = nanmean(trl_data(:,:,0.4*200:0.5*200),3);

    dbdata = trl_data - baselines;
    allPbaselines(id,:) = nanmean(baselines);
   
    leftChoice_db = squeeze(nanmean(nanmean(dbdata(choices == 0,exp.movCluster(1,:),:) - dbdata(choices == 0,exp.movCluster(2,:),:))));
    rightChoice_db = squeeze(nanmean(nanmean(dbdata(choices ==1,exp.movCluster(1,:),:) - dbdata(choices== 1,exp.movCluster(2,:),:))));

    topoData_left(id,:,:,:) = nanmean(dbdata(choices == 0,:,:));
    topoData_right(id,:,:,:) = nanmean(dbdata(choices == 1,:,:));
    
    % Realign samples correctly 
    for i = 1:10
        offset_sampleTimes(i) = round([(0.1*200)+(0.386*i*200)]);
        real_sampleTimes(i) = round([(0.1*200)+(0.4*i*200)]);
    end
    
    %Realign: subjects 13-20 had a technical issue where stimuli lasted 7ms less than they should. 
    %Here, we are "filling-in" the missing 7ms (3 samples) with NaN values, and shifting data to the originally intended time
    %for whole-trial ERP plots.
    if sub > 12 
         ldata = leftChoice_db(:,:,:); 
        rdata = rightChoice_db(:,:,:);

        ldata(:,:,0.1*200:end) = NaN;
        rdata(:,:,0.1*200:end) = NaN;

        for i = 1:10
            toNan = 3; %Add 3 nan samples to align with correctly timed data
            ldata(real_sampleTimes(i):real_sampleTimes(i)+80) = leftChoice_db(offset_sampleTimes(i):offset_sampleTimes(i)+80); 
            rdata(real_sampleTimes(i):real_sampleTimes(i)+80) = rightChoice_db(offset_sampleTimes(i):offset_sampleTimes(i)+80);
        end
        clear leftChoice_db; clear rightChoice_db;  leftChoice_db = ldata; rightChoice_db = rdata; 
        for i = 1:10
            leftChoice_db(real_sampleTimes(i)-toNan:real_sampleTimes(i)-1) = NaN;
            rightChoice_db(real_sampleTimes(i)-toNan:real_sampleTimes(i)-1) = NaN;

        end
    else
        for i = 1:10
            toNan = 3;
            leftChoice_db(real_sampleTimes(i)-toNan:real_sampleTimes(i)-1) = NaN;
            rightChoice_db(real_sampleTimes(i)-toNan:real_sampleTimes(i)-1) = NaN;
        end
     end
    allLeft(id,:) = leftChoice_db; 
    allRight(id,:) = rightChoice_db; 
end


%% Plot
f = figure;

timeLim = [0.4,(4.5+0.4)]*200;
normLeft = nanmean(allLeft(:,timeLim(1):timeLim(1)+0.2*200),2);
stdshade(allLeft(:,timeLim(1):timeLim(2))./normLeft, 0.15, 'k', '-',[],1) %hex2rgb(plt.brightPalettes{3}([2])
hold on;
normRight= nanmean(allRight(:,timeLim(1):timeLim(1)+0.2*200),2);
stdshade(allRight(:,timeLim(1):timeLim(2))./normRight, 0.15, 'k', '--',[],1); hold on; 

for i = 1:11
    sampleTimes(i) = [(0.1*200)+(0.4*(i-1)*200)]
end

resps = {'allLeft', 'allRight'}
for var = 1:2
    thisResp = eval(resps{var});
    amean=squeeze(nanmean(thisResp(:,timeLim(1):timeLim(2))));
    astd=nanstd(thisResp(:,timeLim(1):timeLim(2)))/sqrt(20); % to get sem shading
    
    errorup = [amean+astd];
    errordown = [(amean-astd)];
    for s = 1:11
        x = round(sampleTimes(s)+1:1:sampleTimes(s)+76);
        thisError = [errorup(x), fliplr(errordown(x))];
        fill([x fliplr(x)],thisError, 'k', 'FaceAlpha', 0.15,'linestyle','none'); hold on;
        
    end
    hold on;
end
    
xline(sampleTimes, '--k');
xlim([0,900])
xline([0.1*200])
xlabel('Time (s)')
ylabel({'MBL (dB) [Left-Right]'})

xticks(([0.1,0.5:0.4:4.5]*200));
xticklabels({'Mask','S1', 'S2', 'S3', 'S4', 'S5', 'S6', 'S7', 'S8', 'S9', 'S10',''});

legend({'','Left response','', 'Right response'})
set(gca, 'fontsize', 11);

% ax1 = gca();
% % ADD second X-axis for sample labels
% ax2 = axes('Position', get(ax1,'Position'), ...
%     'XAxisLocation','top', ...
%     'Color','none', ...
%     'XColor','k',...
%     'XLim', [ax1.XLim]);
% ax2.YAxis.Visible = 'off';
% ax2.XAxis.Visible = 'on';
% ax2.XAxis.TickLength = [0,0];
% 
% ax2.XTick = [([0.1,0.5:0.4:4.5]*200)];
% ax2.XTickLabel = {['Trial' ,'\newline', ' start'], 'S1', 'S2', 'S3', 'S4', 'S5', 'S6', 'S7', 'S8', 'S9', 'S10',''}
% xticks(ax2.XTick);
% xticklabels(ax2.XTickLabel)

set(gca, 'fontsize', 10);

% ADD topography 
% Inset topographies 
% Create smaller axes in top left, and plot on it

% axes('Position',[.08 .65 .22 .22])
% box off
% load('C:\Users\elisa\Desktop\VolatilityTask\Analysis\functions\chanlocsBioSemi_128_noext.mat')
% 
% ttopoData = squeeze(nanmean(topoData_left - topoData_right));
% ttopoData = nanmean(ttopoData(:,[(4.2+0.4)*200]:[(4.4+0.4)*200]),2);
% 
% VT_topoplot(ttopoData,chanlocs, ...
%     'electrodes', 'off', ...
%     'emarker', {'.','k',[],5},...
%     'emarker2', {[exp.movCluster(1,:), exp.movCluster(2,:)], '*', 'k',3,0.5},...)
%     'numcontour', 0, ...
%     'plotrad',0.65, ...
%     'headrad', 0.5)

% caxis([-0.6 0.6])

f.Units = 'centimeters';
f.OuterPosition = [0.25 0.35 19 11];
% exportgraphics(f,[exp.figPath, 'Fig2A.tiff'], 'Resolution', 600)

