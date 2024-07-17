% Creates sequence of updated beliefs given some input sequence of log
% likelihood ratios, according the Glaze et al. normative model. H is
% hazard rate of changes in generative distribution, startpoint is
% starting belief in bans (LLR units).

function [LPRout,surprise,scaled_prior] = accGlaze_fast(LLRin,H,startpoint,Stype,pIn)

% Initializing LPR and surprise vectors
if H>0
    LPRout = LLRin(:,1)+startpoint+log(((1-H)/H)+exp(-startpoint))-log(((1-H)/H)+exp(startpoint));
else
    LPRout = LLRin(:,1)+startpoint;
end
surprise = nan(size(LPRout));
scaled_prior = zeros(size(LPRout));


% Run belief updating
if H>0
    for s = 2:size(LLRin,2)
        scaled_prior(:,end+1) = LPRout(:,end)+log(((1-H)/H)+exp(-LPRout(:,end)))-log(((1-H)/H)+exp(LPRout(:,end)));
        LPRout(:,end+1) = LLRin(:,s)+scaled_prior(:,end);
    end
    if s == size(LLRin,2)
        scaled_prior(:,end+1) = LPRout(:,end)+log(((1-H)/H)+exp(-LPRout(:,end)))-log(((1-H)/H)+exp(LPRout(:,end)));
    end

else
    for s = 2:size(LLRin,2)
        scaled_prior(:,end+1) = LPRout(:,end);
        LPRout(:,end+1) = LLRin(:,s)+LPRout(:,end);
    end
end

% Calculate sample-wise surprise
if strcmp(Stype,'pCP')  % analytically-derived change-point probability
    
% CPP was the posterior probability that a change in
% generative task state has just occurred, given the expected H, the evidence
% carried by xn and the observer’s belief before encountering that sample Ln−1 

% First, the numerator computes a weighted sum of the likelihood of the new sample under both generative states assuming that a change has occurred, 
% with each weight determined by the strength of the observer’s existing belief in the opposing state. 
% This means that a new sample of evidence that is inconsistent with the observer’s belief (i.e.
% sign(LLRn)≠sign(Ln-1)) will yield a larger CPP than a sample that is consistent. 
% Second, if the new sample carries no information about the current generative state (i.e. LLRn=0), eq. S3 evaluates to H.
% In other words, when a new sample is ambiguous, the observer must rely more on their base expected rate of state change as an estimate for CPP. 
% Similarly, if the observer is agnostic as to the task state (i.e. Ln-1=0), eq. S3 again evaluates to H. That is, a belief about a change-point having occurred over and above the base expected rate of change can only form if the observer has some level
% of belief in the task state before encountering the new sample.
% (From supplementary note 1 in Murphy et al. 2022, @https://static-content.springer.com/esm/art%3A10.1038%2Fs41593-021-00839-z/MediaObjects/41593_2021_839_MOESM1_ESM.pdf)

    pR = l2p([zeros(size(LPRout,1),1) LPRout],'p'); %Probability sample belongs to Right distribution
    pL = l2p([zeros(size(LPRout,1),1) LPRout],'n'); %Probability sample belongs to Left distribution
    
    for s = 1:size(LLRin,2)
        surprise(:,s) = (H.* ((squeeze(pIn(:,s,1)).*pL(:,s)) + (squeeze(pIn(:,s,2)).*pR(:,s)))) ./ ...
            ((H.* ((squeeze(pIn(:,s,1)).*pL(:,s)) + (squeeze(pIn(:,s,2)).*pR(:,s)))) + ((1-H).*((squeeze(pIn(:,s,2)).*pL(:,s)) + (squeeze(pIn(:,s,1)).*pR(:,s)))));
    end
end
    
end