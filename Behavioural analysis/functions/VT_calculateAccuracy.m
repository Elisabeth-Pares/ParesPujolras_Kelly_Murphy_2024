function [GA] = VT_calculateAccuracy(exp,par,dat,GA)

thisModel = par.thisModel;
subj = par.ss


switch thisModel
    
    case 'data' 
        GA.(thisModel).acc(subj,1) = nanmean(dat.(thisModel).acc_full); 
        cacc = (dat.(thisModel).acc_full); 
          %Sort accuracy by time from final change point, in full trials only
        for b = 1:size(exp.sampBounds,1)
            c = find(dat.(thisModel).fsmps_full>=exp.sampBounds(b,1) & dat.(thisModel).fsmps_full<=exp.sampBounds(b,2));
            GA.(thisModel).acc_fsmp(subj,b) = nanmean(dat.(thisModel).acc_full(c));
        end
        
      case 'ideal' % = normative glaze 
        % Posterior (LPR) at end of trial, following ideal normative computation
        cacc = dat.(thisModel).LPR_final_fsmps; 

        % If final distribution = left, flip LPR sign; all correct
        % decisions (based on LPR) are positive. 
        cacc(sign(dat.(thisModel).fdistc_full)==-1) = cacc(sign(dat.(thisModel).fdistc_full)==-1).*-1; % conditioning dat.(thisModel).LPRs on accuracy, not choice
        
        %ERF gives a probability of choice, based on how extreme the LPR
        %value is. The computation below forces values to be between 0-1,
        %where 0.5 is chance probability of correct/incorrect, 0 = sure
        %error, 1 = sure correct. 
        cacc = 0.5+0.5.*erf(cacc);
        GA.(thisModel).acc(subj,1) = nanmean(cacc);
        %% The acc_binary computation just compares the sign, without
        % considering probability.
        cacc = sign(dat.(thisModel).LPR_final_fsmps);
        GA.(thisModel).acc_binary(subj,1) = length(find(cacc-dat.(thisModel).fdistc_full==0))./length(cacc);
       
        %Sort accuracy by time from final change point, in full trials only
        for b = 1:size(exp.sampBounds,1)
            c = find(dat.(thisModel).fsmps_full>=exp.sampBounds(b,1) & dat.(thisModel).fsmps_full<=exp.sampBounds(b,2));
            GA.(thisModel).acc_fsmp(subj,b) = length(find(cacc(c)-dat.(thisModel).fdistc_full(c)==0))./length(c);
        end
        
    case 'fitted_glaze'
        cacc = dat.(thisModel).LPR_final_fsmps./GA.(thisModel).noise(subj);  %%% Normative model fit

        cacc(sign(dat.(thisModel).fdistc_full)==-1) = cacc(sign(dat.(thisModel).fdistc_full)==-1).*-1; % conditioning dat.(thisModel).LPRs on accuracy, not choice
        cacc = 0.5+0.5.*erf(cacc);
        GA.(thisModel).acc(subj,1) = nanmean(cacc);
       
        %Sort accuracy by time from final change point, in full trials only
        for b = 1:size(exp.sampBounds,1)
            GA.(thisModel).acc_fsmp(subj,b) = nanmean(cacc(dat.(thisModel).fsmps_full>=exp.sampBounds(b,1) & dat.(thisModel).fsmps_full <=exp.sampBounds(b,2)));
        end
        
        
        cacc = sign(dat.(thisModel).LPR_final_fsmps);
        GA.(thisModel).acc_binary(subj,1) = length(find(cacc-dat.(thisModel).fdistc_full==0))./length(cacc);
        
    case 'fitted_leak'
        cacc = dat.(thisModel).LPR_final_fsmps./GA.(thisModel).noise(subj);  %%% Normative model fit
        cacc(sign(dat.(thisModel).fdistc_full)==-1) = cacc(sign(dat.(thisModel).fdistc_full)==-1).*-1; % conditioning dat.(thisModel).LPRs on accuracy, not choice
        cacc = 0.5+0.5.*erf(cacc);
        GA.(thisModel).acc(subj,1) = nanmean(cacc);
       
        %Sort accuracy by time from final change point, in full trials only
        for b = 1:size(exp.sampBounds,1)
            GA.(thisModel).acc_fsmp(subj,b) = nanmean(cacc(dat.(thisModel).fsmps_full>=exp.sampBounds(b,1) & dat.(thisModel).fsmps_full <=exp.sampBounds(b,2)));
        end
        
        
        cacc = sign(dat.(thisModel).LPR_final_fsmps);
        GA.(thisModel).acc_binary(subj,1) = length(find(cacc-dat.(thisModel).fdistc_full==0))./length(cacc);
%         
    case 'perf'
        
        cacc = dat.(thisModel).LPR_final_fsmps;  %%% Normative model fit
        cacc(sign(dat.(thisModel).fdistc_full)==-1) = cacc(sign(dat.(thisModel).fdistc_full)==-1).*-1; % conditioning dat.(thisModel).LPRs on accuracy, not choice
        cacc = 0.5+0.5.*erf(cacc);
        GA.(thisModel).acc(subj,1) = nanmean(cacc);
       
        %Sort accuracy by time from final change point, in full trials only
        for b = 1:size(exp.sampBounds,1)
            GA.(thisModel).acc_fsmp(subj,b) = nanmean(cacc(dat.(thisModel).fsmps_full>=exp.sampBounds(b,1) & dat.(thisModel).fsmps_full <=exp.sampBounds(b,2)));
        end
        
        cacc = sign(dat.(thisModel).LPR_final_fsmps);
        GA.(thisModel).acc_binary(subj,1) = length(find(cacc-dat.(thisModel).fdistc_full==0))./length(cacc);

    case 'last_sample'
        cacc = dat.(thisModel).LLR_final_fsmps;       %%% Last-sample only (no noise)
        cacc(sign(dat.(thisModel).fdistc_full)==-1) = cacc(sign(dat.(thisModel).fdistc_full)==-1).*-1; % conditioning dat.(thisModel).LPRs on accuracy, not choice
        cacc = 0.5+0.5.*erf(cacc);
        GA.(thisModel).acc(subj,1) = nanmean(cacc);
        
        cacc = sign(dat.(thisModel).LLR_final_fsmps);
        GA.(thisModel).acc_binary(subj,1) = length(find(cacc-dat.(thisModel).fdistc_full==0))./length(cacc);
        
        %Sort accuracy by time from final change point, in full trials only
        for b = 1:size(exp.sampBounds,1)
            c = find(dat.(thisModel).fsmps_full>=exp.sampBounds(b,1) & dat.(thisModel).fsmps_full<=exp.sampBounds(b,2));
            GA.(thisModel).acc_fsmp(subj,b) = length(find(cacc(c)-dat.(thisModel).fdistc_full(c)==0))./length(c);
        end
        
    case 'fitted_perf'
        cacc = dat.(thisModel).LPR_final_fsmps;      %%% Perfect accumulation (no noise)
%         cacc = dat.(thisModel).psi_final_fsmps;  %%% Perfect accumulation w/non-absorbing bounds

        cacc(sign(dat.(thisModel).fdistc_full)==-1) = cacc(sign(dat.(thisModel).fdistc_full)==-1).*-1; % conditioning dat.(thisModel).LPRs on accuracy, not choice
        cacc = 0.5+0.5.*erf(cacc);
        GA.(thisModel).acc(subj,1) = nanmean(cacc);
        
        cacc = sign(dat.(thisModel).LPR_final_fsmps);
        GA.(thisModel).acc_binary(subj,1) = length(find(cacc-dat.(thisModel).fdistc_full==0))./length(cacc);
  
        %Sort accuracy by time from final change point, in full trials only
        for b = 1:size(exp.sampBounds,1)
            c = find(dat.(thisModel).fsmps_full>=exp.sampBounds(b,1) & dat.(thisModel).fsmps_full<=exp.sampBounds(b,2));
            GA.(thisModel).acc_fsmp(subj,b) = length(find(cacc(c)-dat.(thisModel).fdistc_full(c)==0))./length(c);
        end
        
        %% Perfect accumulation with Non-absorbing & human noise 
    case 'fitted_perf_naBounds'
        cacc = dat.(thisModel).LPR_final_fsmps./GA.(thisModel).noise(subj);  %%% Perfect accumulation w/non-absorbing bounds

        cacc(sign(dat.(thisModel).fdistc_full)==-1) = cacc(sign(dat.(thisModel).fdistc_full)==-1).*-1; % conditioning dat.(thisModel).LPRs on accuracy, not choice
        cacc = 0.5+0.5.*erf(cacc);
        GA.(thisModel).acc(subj,1) = nanmean(cacc);
        
        cacc = sign(dat.(thisModel).LPR_final_fsmps);

        GA.(thisModel).acc_binary(subj,1) = length(find(cacc-dat.(thisModel).fdistc_full==0))./length(cacc);
  
        %Sort accuracy by time from final change point, in full trials only
        for b = 1:size(exp.sampBounds,1)
            c = find(dat.(thisModel).fsmps_full>=exp.sampBounds(b,1) & dat.(thisModel).fsmps_full<=exp.sampBounds(b,2));
            GA.(thisModel).acc_fsmp(subj,b) = length(find(cacc(c)-dat.(thisModel).fdistc_full(c)==0))./length(c);
        end
        disp('check')
end