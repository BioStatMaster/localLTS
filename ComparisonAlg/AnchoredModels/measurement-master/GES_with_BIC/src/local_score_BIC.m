function score = local_score_BIC(Data,i,PAi) % calculate the local score with BIC for the linear Gaussian case
T=size(Data,1);
X = Data(:,i);

if(~isempty(PAi))
    PA = Data(:,PAi);    
    D = size(PA, 2);
    %% derive the parameters by maximum likelihood
    H = PA*pdinv(PA'*PA)*PA';
    E = X-H*X;
    sigma2 = sum(E.^2)/T;
    %% BIC
    score = T*log(sigma2)+D*log(T);
else
    
    sigma2 = sum(X.^2)/T;
    %% BIC
    score = T*log(sigma2);
    
end








