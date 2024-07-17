function [GA] = VT_eeg2behRegression(exp, par)

smpwin = [0 0.4];  % window for sample-wise analyses
trlwin = [-0.1 5.4];  % window for full-trial-wise analyses

clear GA;

thisPath = ['P01\S2\EEG\'];
fullPath = ([exp.dataPath, thisPath]);
cd(fullPath);
EEG = pop_loadset( [fullPath '\csd_DSacICASL_f01fc_P01_S2.set'] );

if strcmp(par.thisModel, 'CPPRes_choice_timeRes'), par.res= {'timeResolved_frontPostCluster', 'timeResolved_frontalCluster', 'timeResolved_posteriorCluster'};
else par.res = {'tfa_res'};
end

figure;
tiledlayout(4,5)
for sub = [1:20]
    for r = 1:length(par.res)
        thisRes = par.res{r};
        
        subj = exp.sub_id{sub};
        sess_data = [];
        disp(['Processing sub ' num2str(sub) '...']);
        
        if strcmp(par.thisModel, 'BetaRes_choice_timeRes')
            ResLab = 'BetaRes';
            
            %Get indices of bad trials
            %Load sample EEG dataset for timepopints
            %Load standard regression results (EEG ~ Prior + LLR + LLR*Surprise + LLR*Uncertainty + session)
            times = EEG.times/1000;
            
            semifullPath = ([exp.dataPath, 'P' subj '\']);
            %Load results of psi_deltaPsi regression of motor beta
            %lateralisation 
            load([semifullPath '\VT_regression_psi_effEv_TFA_allT_subtr_lateralisationrsq0_P' subj '_2.mat']);
            
            %Average residuals over sample time
            times = times(times>=trlwin(1) & times<=trlwin(2));
            onsets = 0.4:0.4:0.4*10;  onsets=round(onsets,1);  % vector of all sample onset times relative to pre-mask; epoched around pre-mask on
            smptimes = times(times>=smpwin(1) & times<=smpwin(2));  % getting vector of sample times relative to dot onset
            
            smp_data = squeeze(nanmean(dat.prior.r(:,(0.2*200:1*200),:,:),2))'; %Smp_data dimensions are [trials x channels x timepoints x samples]
            
            smp_data_res = dat.prior.r; %Neural residuals
            smp_data_res = squeeze(permute(smp_data_res,[3,2,4,1]));
            
            %Get behaviour
            psi_all = []; surprise_all = []; llr_all =[];
            psi_full = []; llr_full = []; surprise_full = [];
            choices_all = []; fullTrials_all = [];
            
            load(['C:\Users\elisa\Desktop\VolatilityTask\Analysis\ProcStim_fitted_glaze_P' subj '_S2.mat']) % This MDLVARS. comes from here_
            
            choices_full = dat.fitted_glaze.choices_full;
            psi_full = dat.fitted_glaze.psi_full;
            llr_full = dat.fitted_glaze.oLLR_full;
            surprise_full = dat.fitted_glaze.surprise_full;
            fCP_full = dat.fitted_glaze.fCPpos_full;
            
        elseif strcmp(par.thisModel, 'CPPRes_choice_timeRes') 
            ResLab = 'LLRxCPPRes';
            times = EEG.times/1000;
            semifullPath = ([exp.dataPath, 'P' subj '\']);

            %LOAD RESIDUALS
            load([semifullPath '\VT_CPPresiduals_deltaPsi.mat']); %Output of VT_getResiduals "deltaPsi"; 
            smp_data_res = res.(thisRes); %Neural residuals

        end
        
        %Get behaviour
        load(['C:\Users\elisa\Desktop\VolatilityTask\Analysis\ProcStim_fitted_glaze_P' subj '_S2.mat']) % This MDLVARS. comes from here_
        
        choices_full = dat.fitted_glaze.choices_full;
        psi_full = dat.fitted_glaze.psi_full;
        deltaPsi = diff(dat.fitted_glaze.psi_full,1,2);
        llr_full = dat.fitted_glaze.oLLR_full; %Try with fitted values
        surprise_full = dat.fitted_glaze.surprise_full;
        fCP_full = dat.fitted_glaze.fCPpos_full;
    end
    
    if sub == 1;idcs = [225:235];  psi_full(idcs,:) = []; deltaPsi(idcs,:) = []; llr_full(idcs,:) = []; surprise_full(idcs,:) = []; choices_full(idcs,:)=[];fCP_full(idcs,:) = [];
    end
    idcs = [];
    effEv = diff(psi_full,1,2);
    
    assert(length(surprise_full)==size(smp_data_res,3),'ERROR: Trial counts in eeg/behaviour are unequal')
   
    
    %LOGIT TRANSFORM SURPRISE! 
    surprise_full = log(surprise_full./(1-surprise_full)); 
    
    %Loop through samples
    for t = 1:200
        if strcmp(par.thisModel, 'BetaRes_choice_timeRes')
            
            idx = 1:size(llr_full,1);
            
            %% Time-resolved regression, as in Wyart et al. 2015
            X = [nanzscore(llr_full(:,1:10)),... %LLR
                nanzscore(nanzscore(llr_full(:,2:10)).*nanzscore(surprise_full(:,2:10)),0,1),... %SURPRISE
                nanzscore(nanzscore(llr_full(:,2:10)).*nanzscore(-abs(psi_full(:,2:10))),0,1),...;%,...   %UNCERTAINTY
                nanzscore(squeeze(smp_data_res(1:9,t,:))')];
            
        elseif strcmp(par.thisModel, 'CPPRes_choice_timeRes')
            idx = 1:size(llr_full,1);
            
            X = [nanzscore(llr_full(:,1:10)),... %LLR
                nanzscore(nanzscore(llr_full(:,2:10)).*nanzscore(surprise_full(:,2:10)),0,1),... %SURPRISE
                nanzscore(nanzscore(llr_full(:,2:10)).*nanzscore(-abs(psi_full(:,2:10))),0,1),...;%,...   %UNCERTAINTY
                nanzscore(nanzscore(llr_full(:,1:10)).*nanzscore(squeeze(smp_data_res(1:10,t,:))'))];%,...
        end
        
        Y = choices_full(idx,:);
        
        %Exclude any wrong keypresses
        toExclude = find(Y~= 1 & Y~=0);
        Y(toExclude,:) = [];
        X(toExclude,:) = [];
                
        X = nanzscore(X);
        
        [B,dev,stats] = glmfit(X,Y,'binomial');
        
        if strcmp(par.thisModel, 'CPPRes_choice_timeRes')
            GA.(par.thisModel).(thisRes).LLR(sub,:,t) = B(2:exp.nSamples+1);          
            GA.(par.thisModel).(thisRes).LLRxSurprise(sub,:,t) =  B(exp.nSamples+2:(2*exp.nSamples));           
            GA.(par.thisModel).(thisRes).LLRxUncertainty(sub,:,t) = B((2*exp.nSamples)+1:((3*exp.nSamples))-1);       
            GA.(par.thisModel).(thisRes).Res(sub,:,t) = B((3*exp.nSamples):((4*exp.nSamples))-1); 
        elseif  strcmp(par.thisModel, 'BetaRes_choice_timeRes')
            GA.(par.thisModel).(thisRes).LLR(sub,:,t) = B(2:exp.nSamples+1); 
            GA.(par.thisModel).(thisRes).LLRxSurprise(sub,:,t) =  B(exp.nSamples+2:(2*exp.nSamples)); 
            GA.(par.thisModel).(thisRes).LLRxUncertainty(sub,:,t) = B((2*exp.nSamples)+1:((3*exp.nSamples))-1);
            GA.(par.thisModel).(thisRes).Res(sub,:,t) = B((3*exp.nSamples):((4*exp.nSamples)-2)); 
        end
    end
    
    clear smp_data; clear toExclude; clear sess_data; clear trl_data; clear fullTrials;
    clear X; clear Y;
    clear smp_data; clear llr_full; clear surprise_full; clear psi_full;
    clear fsmp_data_res;
end


save([exp.dataPath, 'GA_' par.thisModel '_regression.mat'], 'GA')
