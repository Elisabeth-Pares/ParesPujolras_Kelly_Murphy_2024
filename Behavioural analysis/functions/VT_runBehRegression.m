function [GA, par] = VT_runBehRegression(exp,par,dat,GA)

surprise = log(dat.(par.thisModel).surprise_full./(1-dat.(par.thisModel).surprise_full)); 

X = [nanzscore(dat.(par.thisModel).oLLR_full(:,1:end),0,1),... %LLR
    nanzscore(dat.(par.thisModel).oLLR_full(:,2:end).*nanzscore(surprise(:,2:end)),0,1),... %SURPRISE
    nanzscore(dat.(par.thisModel).oLLR_full(:,2:end).*nanzscore(-abs(dat.(par.thisModel).psi_full(:,2:10))),0,1)];   %UNCERTAINTY

%Define DV
if strcmp(par.thisModel, 'ideal') || strcmp(par.thisModel, 'data') || strcmp(par.thisModel, 'perf') %Ideal == norm_glaze  strcmp(par.thisModel, 'norm_glaze') ||
    Y = [dat.(par.thisModel).choices_full ones(length(dat.(par.thisModel).choices_full),1)]; %Second (optional) column indicates that each row is one trial; could be eliminated.
    Y(Y~=0& Y~=1,1) = NaN; %If wrong key, replace with NaN; 
else
    %Any other fitted models
    % Compute kernels/LLR regresss for model fits
    % In methods: We then fit variants of the normative model to participants’ behavior, apssuming
    % that choices were based on the log-posterior odds Ln,trl [LRP_full] for the observed stimulus
    % sequence on each trial. Ln,trl was corrupted by a noise term ν, such that choice
    % probability ˆr was computed as follows.
    % erf(x) = 2/sqrt(pi) * integral from 0 to x of exp(-t^2) dt.
    
    CP = 0.5+0.5.*erf(dat.(par.thisModel).LPR_final_full./GA.(par.thisModel).noise(par.ss));  % calculate choice probabilities scaled by noise

    Y = [CP ones(length(dat.(par.thisModel).choices_full),1)];
end

%Run regression 
[B,~,stats] = glmfit(X,Y,'binomial');

%Save results
%B returns 29 elements:
if strcmp(par.regression, 'Murphy')
    GA.(par.thisModel).B_LLR(par.ss,:) = B(2:exp.nSamples+1); %         2-11 LLR for samples 1-10; old name = par.ssB_regU
    GA.(par.thisModel).B_surprise(par.ss,:) =  B(exp.nSamples+2:(2*exp.nSamples)); % %         2-11 LLR for samples 1-10; old name = GA_par.ssB_surpU
    GA.(par.thisModel).B_uncertainty(par.ss,:) = B(end-length(2:exp.nSamples)+1:end); %         2-11 LLR for samples 1-10; old name = GA_par.ssB_uncertU
end

% %Plot single subject beta scores (EPP)
% figure; plot(B(2:exp.nSamples+1), 'Color', 'b');
% hold on; plot(B(exp.nSamples+2:end-length(2:exp.nSamples)), 'Color', 'r');
% hold on; plot(B(end-length(2:exp.nSamples)+1:end), 'Color', 'g');
% ylabel('Beta'); yline([0 0]); xlabel('Sample position'); title(par.thisModel)
% legend({'LLR','surprise','uncertainty','',''}, 'Location', 'northwest')

par.Y = Y; 
        
end