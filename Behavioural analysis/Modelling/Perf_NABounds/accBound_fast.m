% Creates sequence of updated beliefs given some input sequence of log
% likelihood ratios, according the Glaze et al. normative model. H is
% hazard rate of changes in generative distribution, startpoint is
% starting belief in bans (LLR units).

function [LPRout,surprise,prior] = accBound_fast(LLRin,A,startpoint,Stype,pIn,H)

% Initializing LPR and surprise vectors
LPRout = LLRin(:,1)+startpoint;
surprise = nan(size(LPRout));
prior = ones(size(LPRout,1),1).*startpoint;

% Run belief updating %
for s = 2:size(LLRin,2)
    prior(:,end+1) = LPRout(:,end);
    prior(LPRout(:,end)>A,end) = A; prior(LPRout(:,end)<-A,end) = -A;   % implementing bound transformation *before* accumulation so final LPR is equivalent to Glaze (i.e. untransformed)    
    LPRout(:,end+1) = prior(:,end)+LLRin(:,s);

end
%Test
if s == size(LLRin,2)
    prior(:,end+1) = LPRout(:,end);
    prior(LPRout(:,end)>A,end) = A; prior(LPRout(:,end)<-A,end) = -A;   % implementing bound transformation *before* accumulation so final LPR is equivalent to Glaze (i.e. untransformed)
end

% Calculate sample-wise surprise
if strcmp(Stype,'pCP')  % analytically-derived change-point probability
    pR = l2p([zeros(size(LPRout,1),1) LPRout],'p');
    pL = l2p([zeros(size(LPRout,1),1) LPRout],'n');
    
    for s = 1:size(LLRin,2)
        surprise(:,s) = (H.* ((squeeze(pIn(:,s,1)).*pL(:,s)) + (squeeze(pIn(:,s,2)).*pR(:,s)))) ./ ...
            ((H.* ((squeeze(pIn(:,s,1)).*pL(:,s)) + (squeeze(pIn(:,s,2)).*pR(:,s)))) + ((1-H).*((squeeze(pIn(:,s,2)).*pL(:,s)) + (squeeze(pIn(:,s,1)).*pR(:,s)))));
    end

end
