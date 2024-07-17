%% Behavioural fits
% Load behavioural data and variety of model fits, and estimate desired
% behavioural measures (accuracies, psychophysical kernels, etc.).

function [GA] = VT_pullBehFits_optimised(nsub, exper, par,GA)
Htrue = exper.H;

thisModel = par.thisModel;
sub = [exper.sub_id{nsub}];
subj = ['P' exper.sub_id{nsub}];
thisPath = [exper.dataPath, subj];
fprintf(['Processing ' subj])
sdirs = dir([thisPath]);

dat = struct();
for s = 1:2 %Loop through sessions
    
    % Load behavioural data & calculate model-derived belief updates
    behPath = [thisPath, '\S' num2str(s) '\Behaviour\'];
    samplePath = [thisPath, '\S' num2str(s) '\Sample_seqs\'];
    fdirs = dir([behPath]);
    
    for b = 1:(length(fdirs)-2)  % looping through each block within this meg dataset
        
        load([behPath,sub,'_',num2str(s),'_',num2str(b),'.mat'])
        load([samplePath,sub,'_',num2str(s),'_',num2str(b),'.mat'])
        
        % Converting sample and choice values to appropriate signs for choice regressions
        stimIn = round(stimIn.*-1); %Round and change sign
        choices = Behav(:,2)-1;
        acc = Behav(:,3);
        
        % Calculate LLR
        %  Log(P(GenDist = 2)/P(GenDist = 1))
        %  Y = normpdf(X,MU,SIGMA) returns the pdf of the normal distribution with
        %  mean MU and standard deviation SIGMA, evaluated at the values in X.
        %  This effectively transforms the raw stimulus values (stimIn)into
        %  probability values, based on the probability distributions
        %  mean & SD.
        % Convert stimulus values to LLRs & calculate sample-wise surprise
        oLLR = log(normpdf(stimIn,gen.mu(2)*-1,gen.sigma(2))./normpdf(stimIn,gen.mu(1)*-1,gen.sigma(1)));
        
        if ~isfield(dat, thisModel), dat.(thisModel) = struct(); end %Initialise model field
        for f = 1:length(par.datFields)%Initialise data fields
            thisField = par.datFields{f};
            if ~isfield(dat.(thisModel), thisField), dat.(thisModel).(thisField) = []; end
        end
        switch thisModel
            case {'data', 'last_sample'}
                pIn = cat(3,normpdf(stimIn,26,29),normpdf(stimIn,-26,29)); %pCP used as surprise metric
                [LPR, surprise, psi] = accGlaze_fast(oLLR,gen.H,0,'pCP',pIn);%accPerf_fast
%                 [LPR, surprise, psi] = accPerf_fast(oLLR,0,'pCP',pIn, Htrue);%accPerf_fast
%                   LPR= oLLR; %This will yield surprise & psi like ideal obs, but real LLRs. %non transformed; use function only to compute surprise & psi
                thisLLR = oLLR;
            case 'ideal' % norm_glaze
                %% NORMATIVE GLAZE - OBJECTIVE model parameters, without fitting to subject.
                %                     [LPR,surprise,psi] = accGlaze_fast(oLLR,gen.H,0,par.stype,[]); %4th input argument defines how to calculate surprise.
                pIn = cat(3,normpdf(stimIn,26,29),normpdf(stimIn,-26,29)); %pCP used as surprise metric
                [LPR,surprise,psi] = accGlaze_fast(oLLR,gen.H,0,'pCP', pIn);
                thisLLR = oLLR; 
            case 'fitted_glaze'
                %% FITTED GLAZE - transform to SUBJECTIVE model values, based on subject fits.
                load([exper.modPath,'/Fits/',thisModel '/' subj,'_fixed_fit.mat'])  % load model fit data
                GA.(thisModel).H(nsub,1) = pm_fit(1);
                GA.(thisModel).gainB(nsub,1) = pm_fit(2);
                GA.(thisModel).noise(nsub,1) = pm_fit(3);
                
                LLRinMlin = oLLR.*pm_fit(2);  % calculate model-derived sample-wise LLRs
                
                c = 2; while stimIn(1,c)==0, c=c+1; end  % get sample ~= zero
                cLLR = LLRinMlin(1,c); cStim = stimIn(1,c);
                
                sigma = sqrt((-(cStim-26).^2 + (cStim+26).^2)./(2.*cLLR));
                pIn = cat(3,normpdf(stimIn,26,sigma),normpdf(stimIn,-26,sigma));
                H = pm_fit(1);
                
                [LPR,surprise,psi] = accGlaze_fast(LLRinMlin,H,0,'pCP',pIn);
                thisLLR = LLRinMlin;
            case 'fitted_perf'
                %% PERFECT ACCUMULATION
                % --> all evidence has the same weight
                load([exper.modPath,'/Fits/fitted_perf/' subj,'_fixed_fit.mat'])  % load model fit data
                
                GA.(thisModel).gainB(nsub,1) = pm_fit(1);
                GA.(thisModel).noise(nsub,1) = pm_fit(2);
                
                LLRinM_perf = oLLR.*pm_fit(1);
                c = 2; while stimIn(1,c)==0, c=c+1; end  % get sample ~= zero
                cLLR = LLRinM_perf(1,c); cStim = stimIn(1,c);
               
                sigma = sqrt((-(cStim-26).^2 + (cStim+26).^2)./(2.*cLLR));
                pIn_perf= cat(3,normpdf(stimIn,26,sigma),normpdf(stimIn,-26,sigma));
                
                [LPR,surprise,psi] = accPerf_fast(LLRinM_perf,0,'pCP',pIn_perf,Htrue); %Fitting with H true now; NOT subjective H; this would be> GA.fitted_glaze.H(nsub)
                thisLLR = LLRinM_perf;
                
            case 'perf'
                %% NOISELESS PERFECT ACCUMULATION
                
                GA.(thisModel).gainB(nsub,1) = 1;%pm_fit(1);
                GA.(thisModel).noise(nsub,1) = 0;%pm_fit(2);
                
                LLRinM_perf = oLLR;
                c = 2; while stimIn(1,c)==0, c=c+1; end  % get sample ~= zero
                cLLR = LLRinM_perf(1,c); cStim = stimIn(1,c);
                
%                 cLLR = LLRinM_perf(1,5); cStim = stimIn(1,5);
                sigma = sqrt((-(cStim-26).^2 + (cStim+26).^2)./(2.*cLLR));
                pIn_perf= cat(3,normpdf(stimIn,26,sigma),normpdf(stimIn,-26,sigma));
                
                [LPR,surprise,psi] = accPerf_fast(LLRinM_perf,0,'pCP',pIn_perf,Htrue); %Fitting with H true now; NOT subjective H; this would be> GA.fitted_glaze.H(nsub)
                thisLLR = LLRinM_perf;

            case 'fitted_perf_naBounds'
                %% PERFECT ACCUMULATION
                % --> all evidence has the same weight
                load([exper.modPath,'/Fits/fitted_perf_naBounds/' subj,'_fixed_fit.mat'])  % load model fit data
                
                GA.(thisModel).A(nsub,1) = pm_fit(1);
                GA.(thisModel).gainB(nsub,1) = pm_fit(2);
                GA.(thisModel).noise(nsub,1) = pm_fit(3);

                LLRinM_perf = oLLR.*pm_fit(2);
                c = 2; while stimIn(1,c)==0, c=c+1; end  % get sample ~= zero
                cLLR = LLRinM_perf(1,c); cStim = stimIn(1,c);
                
                sigma = sqrt((-(cStim-26).^2 + (cStim+26).^2)./(2.*cLLR));
                pIn_perf= cat(3,normpdf(stimIn,26,sigma),normpdf(stimIn,-26,sigma));
                
                [LPR,surprise,psi] = accBound_fast(LLRinM_perf,pm_fit(1),0,'pCP',pIn_perf,Htrue); %Perfect accumulation towards non-absorbing bounds
                thisLLR = LLRinM_perf;

            case 'fitted_leak'
                %% LEAKY ACCUMULATION
                % --> all evidence has the same weight
                load([exper.modPath,'/Fits/fitted_leak/' subj,'_fixed_fit.mat'])  % load model fit data
                
                GA.(thisModel).leak(nsub,1) = pm_fit(1);
                GA.(thisModel).gainB(nsub,1) = pm_fit(2);
                GA.(thisModel).noise(nsub,1) = pm_fit(3);

                LLRin_leak= oLLR.*pm_fit(2);
                c = 2; while stimIn(1,c)==0, c=c+1; end  % get sample ~= zero
                cLLR = LLRin_leak(1,c); cStim = stimIn(1,c);
                sigma = sqrt((-(cStim-26).^2 + (cStim+26).^2)./(2.*cLLR));
                pIn_leak= cat(3,normpdf(stimIn,26,sigma),normpdf(stimIn,-26,sigma));
                thisLLR = LLRin_leak;

                [LPR,surprise,psi] = accLeak_fast(LLRin_leak,pm_fit(1),0,'pCP',pIn_leak,Htrue); %Perfect accumulation towards non-absorbing bounds
        end
        % Pulling samples per trial and position of final change-point
        nsamps=[]; fCPpos=[];
        for t = 1:size(oLLR,1)
            nsamps(t,1) = length(find(~isnan(oLLR(t,:))));
            if isempty(find(pswitch(t,:)==1, 1))
                fCPpos(t,1) = 1;
            else fCPpos(t,1) = find(pswitch(t,:)==1,1,'last');
            end
        end
        
        % Pulling # samples from last change-point
        fsmps = (exper.nSamples-fCPpos)+1;  % full-length no-CP trials will take value of 11 %changed from P, those would be 13 in his task
        fsmps(fsmps==11) = 10;
        
        % Isolating useable full-sequence trials
        %ts: Trial Selection - full sequence
        %tsCE: Trial Selection Correct Error - if response value
        %>=2, wrong key
        ts=[]; tsCE=[];
        for t = 1:length(choices)
            if sum(isnan(stimIn(t,:)))==0 && choices(t)<2, ts(end+1) = t; end   % trials with full sequence (ts)
            if choices(t)<2, tsCE(end+1) = t; end   % all trials, regardless of sequence length
        end
        fCPpos_full = fCPpos(ts,:);
        
       
        subinds = sub2ind(size(LPR),1:size(LPR,1),nsamps');
        LPRfinal = LPR(subinds)'; %LPR at end of trial
        LLRfinal = oLLR(subinds)'; %LLR at enf of trial
        subinds = sub2ind(size(psi),1:size(psi,1),(nsamps+1)');
        psifinal = psi(subinds)'; %Psi at end of trial
        
        dat.(thisModel).oLLR= [dat.(thisModel).oLLR; oLLR]; %OBJECTIVE LLR; same for all models
        dat.(thisModel).oLLR_full = [dat.(thisModel).oLLR_full; oLLR(ts,:)]; %OBJECTIVE LLR; same for all models
        
        dat.(thisModel).LLR = [dat.(thisModel).LLR; thisLLR]; %SUBJECTIVE LLR oLLR*gain
        dat.(thisModel).LPR = [dat.(thisModel).LPR; LPR]; %SUBJECTIVE LPR (posterior)
        dat.(thisModel).surprise = [dat.(thisModel).surprise; surprise];
        dat.(thisModel).psi = [dat.(thisModel).psi; psi];
        
        dat.(thisModel).nCP_full = [dat.(thisModel).nCP_full; nansum(pswitch(ts,:),2)];
        dat.(thisModel).nCP = [dat.(thisModel).nCP; nansum(pswitch(:,:),2)];

        dat.(thisModel).fCPpos = [dat.(thisModel).fCPpos; fCPpos]; %Final change point
        dat.(thisModel).fCPpos_full = [dat.(thisModel).fCPpos_full; fCPpos_full]; %Final change point
        
        dat.(thisModel).psi_final_fsmps= [dat.(thisModel).psi_final_fsmps; psifinal(tsCE)];
        dat.(thisModel).LPR_final_fsmps = [dat.(thisModel).LPR_final_fsmps; LPRfinal(tsCE)]; %
        dat.(thisModel).LLR_final_fsmps = [dat.(thisModel).LLR_final_fsmps; LLRfinal(tsCE)];
        dat.(thisModel).final_fsmps = [dat.(thisModel).final_fsmps; fsmps(tsCE,:)]; %Samples from last change point
        
        dat.(thisModel).psi_final= [dat.(thisModel).psi_final; psifinal]; %subjective LPR at end of trial, full trials
        dat.(thisModel).LPR_final= [dat.(thisModel).LPR_final; LPRfinal]; %subjective LPR at end of trial, full trials
        dat.(thisModel).LPR_final_full= [dat.(thisModel).LPR_final_full; LPRfinal(ts)]; %subjective LPR at end of trial, full trials
        dat.(thisModel).LLR_full = [dat.(thisModel).LLR_full; thisLLR(ts,:)]; %subjective LLR at the end of trial
        dat.(thisModel).LPR_full = [dat.(thisModel).LPR_full; LPR(ts,:)]; %subjective LPR at the end of trial
        dat.(thisModel).surprise_full = [dat.(thisModel).surprise_full; surprise(ts,:)];
        dat.(thisModel).psi_full = [dat.(thisModel).psi_full; psi(ts,:)];
        
        dat.(thisModel).choices = [dat.(thisModel).choices; choices];
        dat.(thisModel).choices_full = [dat.(thisModel).choices_full; choices(ts)];

        dat.(thisModel).fdist_full = [dat.(thisModel).fdist_full; Behav(tsCE,1)-1];
        
        dat.(thisModel).acc = [dat.(thisModel).acc; acc];
        dat.(thisModel).acc_full = [dat.(thisModel).acc_full; acc(tsCE)];
%         dat.(thisModel).acc_full = [dat.(thisModel).acc_full; acc(ts)];

        dat.(thisModel).fsmps= [dat.(thisModel).fsmps; fsmps];
        dat.(thisModel).fsmps_full = [dat.(thisModel).fsmps_full; fsmps(tsCE)];
        
        %Flip sign
        dat.(thisModel).fdistc_full = dat.(thisModel).fdist_full;
        dat.(thisModel).fdistc_full(dat.(thisModel).fdistc_full == 0) = -1;
        
        dat.(thisModel).finalDist_full = [dat.(thisModel).finalDist_full; Behav(ts,1)-1];
        
    end
end

% %Overwrite if any full trials missing 
fullTrials = ~isnan(dat.(thisModel).oLLR(:,10));
if sum(fullTrials) ~= size(dat.(thisModel).LPR_final_full,1)
    dat.(thisModel).fCPpos_full = [dat.(thisModel).fCPpos(fullTrials,:)];%subjective LLR at the end of trial
    dat.(thisModel).nCP_full = [dat.(thisModel).nCP(fullTrials,:)];%subjective LLR at the end of trial
    dat.(thisModel).oLLR_full = [dat.(thisModel).oLLR(fullTrials,:)];%subjective LLR at the end of trial
    dat.(thisModel).LLR_full = [dat.(thisModel).LLR(fullTrials,:)];%subjective LLR at the end of trial
    dat.(thisModel).surprise_full = [dat.(thisModel).surprise(fullTrials,:)];
    dat.(thisModel).psi_full = [dat.(thisModel).psi(fullTrials,:)];
    dat.(thisModel).LPR_full = [dat.(thisModel).LPR(fullTrials,:)]; %subjective LPR at the end of trial
    dat.(thisModel).LPR_final_full = [dat.(thisModel).LPR_final(fullTrials,:)]; %subjective LPR at the end of trial
    dat.(thisModel).choices_full = [dat.(thisModel).choices(fullTrials,:)]; %subjective LPR at the end of trial
    dat.(thisModel).acc_full = [dat.(thisModel).acc(fullTrials)];
    dat.(thisModel).fsmps_full = [dat.(thisModel).fsmps(fullTrials)];
end
assert(sum(fullTrials) == size(dat.(thisModel).LPR_final_full,1))

if strcmp(thisModel, 'ideal')
    GA.Normative_surprise_full{nsub}= dat.(thisModel).surprise_full; 
    GA.Normative_psi_full{nsub} = dat.(thisModel).psi_full; 
end
save(['C:\Users\elisa\Desktop\VolatilityTask\Analysis\ProcStim_' thisModel '_P' num2str(sub) '_S' num2str(s) '.mat'], 'dat');

if ~strcmp(thisModel, 'data') & ~strcmp(thisModel, 'perf') & ~strcmp(thisModel, 'last_sample') & ~strcmp(thisModel, 'ideal')
    
    %Compute BIC
    k = length(pm_fit); % number of free parameters
    n = size(dat.(thisModel).oLLR,1);% number of trials
    
    BIC = 2*err + k*log(n);
    
    dat.(thisModel).BIC = BIC;
end


% % Append sample of last change
% % Load processed data info
% load(['C:/Users/elisa/Desktop/VolatilityTask/Modelling/PreprocessedStim/PreprocData_' sub '-.mat'])

%Change time indices 1 - early, 2-mid, 3-late.
changeTime = dat.(thisModel).fCPpos_full;
changeTime(changeTime <=5) = 1; %Early
changeTime(changeTime>5 & changeTime<8) = 2; %Mid
changeTime(changeTime>=8) = 3; %Late;

if length(dat.(thisModel).choices_full) == length(changeTime)
    dat.(thisModel).changeTime = repmat(changeTime,1,10);
    dat.(thisModel).lastChange = repmat(dat.(thisModel).fCPpos_full,1,10);
    dat.(thisModel).finalDist = repmat(dat.(thisModel).finalDist_full,1,10);
else warning('Mismatch!')
end

%% Running LLR/SURPRISE/UNCERTAINTY regressions of choice, for all models
par.ss = nsub;

% Run full regression:
% Choice ~ LLR + LLR*Surprise + LLR*Uncertainty (phi)
if ~strcmp(par.thisModel, 'last_sample')
    [GA, par] = VT_runBehRegression(exper,par,dat,GA)
end

% Calculating accuracy per subject & as a function of # samples from final change-point
% Accuracy --
GA = VT_calculateAccuracy(exper,par,dat,GA)


%Add BIC for each model 
if ~strcmp(thisModel, 'data') & ~strcmp(thisModel, 'perf') & ~strcmp(thisModel, 'last_sample') & ~strcmp(thisModel, 'ideal')
    GA.(thisModel).BIC(nsub,1) = dat.(thisModel).BIC;
end

if strcmp(par.thisModel, 'data')
GA.(thisModel).nCP(nsub,1) = sum(dat.data.nCP== 0)/size(dat.data.nCP,1); 
GA.(thisModel).nCP(nsub,2) = sum(dat.data.nCP== 1)/size(dat.data.nCP,1); 
GA.(thisModel).nCP(nsub,3) = sum(dat.data.nCP> 1)/size(dat.data.nCP,1);  
end

end

