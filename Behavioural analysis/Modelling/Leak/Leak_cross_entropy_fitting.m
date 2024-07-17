% Takes sequences of samples, computes trial-wise choice probabilities
% using bounded accumulator model with specified non-absorbing bound (A),
% generative distribution  gain factor (B) and noise, and computes
% cross-entropy error term (e) with respect to observed responses.

% pm(1:3) = [leak, B, noise]

function e = Leak_cross_entropy_fitting(pm)

% Retrieving global variables from initializing script
global LLRin choices nsamps

% Looping through PSO particles
for p = 1:size(pm,1);
    % Generate choice probabilities given current parameter set
    CPs = Leak_sim_fast(LLRin,nsamps,pm(p,1),pm(p,2),pm(p,3));
    CPs(CPs==1) = 0.9999; CPs(CPs==0) = 0.0001;  % adjusting extreme CPs that lead to inf cross-entropy values
    
    % Calculate cross-entropy between model choice probabilities and observed responses
    e(p,1) = -sum((ones(size(choices))-choices).*log(1-CPs)+(choices.*log(CPs)));
    
end