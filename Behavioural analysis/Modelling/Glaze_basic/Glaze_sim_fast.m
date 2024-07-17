% Script for simulating behaviour generated by Glaze et al. normative model
% during two-alternative sequential evidence accumulation task where
% samples are drawn from one of 2 distributions with same std but different
% means, and mean can switch at fixed hazard rate.

% H = hazard rate
% B = gain term applied to sample-wise LLRs
% noise = divisive scaling factor applied to final belief


function CP = Glaze_sim_fast(LLRin,nsamps,H,B,noise)

% Apply gain term to LLRs (represents subjective component of generative variance)
LLRin = LLRin.*B;

% Increment belief for each sample
LPRout = zeros(size(LLRin,1),1);
if H>0
    for s = 1:size(LLRin,2)
        LPRout(:,end+1) = LLRin(:,s)+LPRout(:,end)+log(((1-H)/H)+exp(-LPRout(:,end)))-log(((1-H)/H)+exp(LPRout(:,end)));
    end
else
    for s = 1:size(LLRin,2)
        LPRout(:,end+1) = LLRin(:,s)+LPRout(:,end);
    end
end

% Retrieve final LPRs for each trial (accounting for variable sequence lengths)
subinds = sub2ind(size(LPRout),1:size(LPRout,1),nsamps'+1);
LPRfinal = LPRout(subinds)';

% Calculate choice probabilities based on LPR @ end of each sequence, scaled by noise
CP = 0.5+0.5.*erf(LPRfinal./noise);