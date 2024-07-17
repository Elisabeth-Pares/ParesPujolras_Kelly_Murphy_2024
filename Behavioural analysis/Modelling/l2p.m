function pout = l2p(Lin,dir)
%LLR to direction probability 

if strcmp(dir,'n') %Negative --> LEFT
    pout = 1./(exp(Lin)+1);
elseif strcmp(dir,'p') %Positive --> RIGHT 
    pout = exp(Lin)./(exp(Lin)+1);
end