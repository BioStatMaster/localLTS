function [A, sigma, w, loglAll, posMeanE] = ifaEM_bk1(X, numGauss, parsEM)

% Em algorithm for independent factor analysis, e as latent variable, no update for mu
% Inputs:
%       X = data points
%       numGauss = number of Gaussian for each component, p
%       parsEM = parameters for EM

% Outputs:
%       A = mixing matrix
%       sigma = mixture of gaussian stds for all components
%       mu = mixture of gaussian means
%       w = mixture of gaussian weights
% (c) Code written by Mingming Gong.
%     Only to be used for academic purposes.

[dim, N]=size(X);
thres = parsEM.thres;

% initialize parameters
mu = parsEM.mu;
sigma = parsEM.sigma;
w = parsEM.w;
Dim = size(w,1);
Lambda = parsEM.noise*eye(dim);

A = parsEM.A;
%% E step
qidMat = qidMatrix(numGauss, 1, Dim);
priorq = priorLatentQ(w, Dim, numGauss, 1, qidMat);

% posterior of q and e
posq = zeros(N, numGauss.^(Dim));
posMeanE = zeros(Dim, N);
posCovE = zeros(Dim, Dim, N);
posMeanEqi = zeros(Dim, numGauss, N);
posCovEqi = zeros(Dim, numGauss, N);
posCovEq = zeros(Dim, Dim, numGauss^Dim, N);
marginalx = zeros(N, 1);
for t = 1:N  %can be paralled here, use parfor
    [posq(t,:), marginalx(t)] = evaluateQ(X(:,t), A, mu, sigma, Dim, numGauss, qidMat, priorq, Lambda);
    [posMeanE(:,t), posCovE(:,:,t), posMeanEqi(:,:,t), posCovEqi(:,:,t), posCovEq(:,:,:,t)] = ...
        evaluateE(X(:,t), posq(t,:), A, Lambda,...
        mu, sigma, Dim, numGauss, qidMat);
end
% fVal = stockEmLowerBound(X, Atilde(:,3:4), posq, posMeanE, posCovE, posCovEqi, posCovEq, w,...
% sigmas, numLag, numGauss, Lambda, qidMat);

% update parameters
logl = -sum(log(marginalx));
logl_prev = logl;
iter = 0;
fprintf('iter%d: negative loglik: %.10f\n', iter, logl);
A
mu
sigma
w
iter = iter + 1;
loglAll = [];
loglAll = [loglAll, logl];

eta = 1;
if parsEM.fast == 0
    alpha = 1;
else
    alpha = 1.2;
end
A_Prev = A;
while iter == 1 || iter < parsEM.maxIter && abs(logl_prev-logl) >= abs(logl_prev)*thres
    
    delta = (logl_prev-logl)/abs(logl_prev);
    
    if delta < thres && iter~=1 && parsEM.fast == 1
        eta = 1;
        A = A_EM;
        qidMat = qidMatrix(numGauss, 1, Dim);
        priorq = priorLatentQ(w, Dim, numGauss, 1, qidMat);
        
        % posterior of q and e
        for t = 1:N
            [posq(t,:), marginalx(t)] = evaluateQ(X(:,t), A, mu, sigma, Dim, numGauss, qidMat, priorq, Lambda);
            [posMeanE(:,t), posCovE(:,:,t), posMeanEqi(:,:,t), posCovEqi(:,:,t), posCovEq(:,:,:,t)] = ...
                evaluateE(X(:,t), posq(t,:), A, Lambda,...
                mu, sigma, Dim, numGauss, qidMat);
        end
        logl = -sum(log(marginalx));
        fprintf('iter%d: negative loglik: %.10f\n', iter, logl);
        A
        mu
        sigma
        w
        loglAll(end) = logl;
        logl_prev = logl;
    else
        eta = eta * alpha;
        logl_prev = logl;
    end
    %% M step
    % update A
    A_EM = (X*posMeanE')/(sum(posCovE,3));
    A = A_Prev + eta*(A_EM-A_Prev);
    A_Prev = A;
    
    if parsEM.updatePrior == 1
        posqi = marginalPosQ(posq, Dim, 1, numGauss, qidMat, 0);
        % update mu
        if mu~=zeros(Dim,numGauss)
            mu = updateMu(posqi, posMeanEqi);
        end
        % update sigma
        sigma = updateSigma(mu, posqi, posCovEqi);
        
        % update w
        w = updateW(posqi);
    end
    
    %% E step
    qidMat = qidMatrix(numGauss, 1, Dim);
    priorq = priorLatentQ(w, Dim, numGauss, 1, qidMat);
    
    % posterior of q and e, can be paralled here, use parfor
    for t = 1:N
        [posq(t,:), marginalx(t)] = evaluateQ(X(:,t), A, mu, sigma, Dim, numGauss, qidMat, priorq, Lambda);
        [posMeanE(:,t), posCovE(:,:,t), posMeanEqi(:,:,t), posCovEqi(:,:,t), posCovEq(:,:,:,t)] = ...
            evaluateE(X(:,t), posq(t,:), A, Lambda,...
            mu, sigma, Dim, numGauss, qidMat);
    end
    logl = -sum(log(marginalx));
    fprintf('iter%d: negative loglik: %.10f\n', iter, logl);
    A
    mu
    sigma
    w
    iter = iter + 1;
    loglAll = [loglAll, logl];
end


function w = updateW(posqi)
% update w

% the posterior of qi
w = sum(posqi,3)/size(posqi,3);


function sigma = updateSigma(mu, posqi, posCovEqi)
% update sigma

% the posterior of qi
posCovEqi = sum(posCovEqi,3);
posqi = sum(posqi,3);
sigma = posCovEqi./(posqi+eps^20)-mu.^2;
sigma = sigma.^0.5;


function mu = updateMu(posqi, posMeanEqi)
% update mu

% the posterior of qi
posMeanEqi = sum(posMeanEqi,3);
posqi = sum(posqi,3);
mu = posMeanEqi./(posqi+eps^20);


function [posMeanE, posCovE, posMeanEqi, posCovEqi, posCovEq1] = evaluateE(yt, posq, A, Lambda, ...
    mu, sigma, Dim, numGauss, qidMat)
% Calculate the posterior of E
posMeanE = zeros(Dim, numGauss^(Dim));
posCovEq = zeros(Dim, Dim, numGauss^(Dim));
posCovEq1 = zeros(Dim, Dim, numGauss^(Dim));
posCovEqi = zeros(Dim, numGauss);
posMeanEqi = zeros(Dim, numGauss);

for q = 1:size(posMeanE,2)
    muq = zeros(1,size(qidMat,2));
    sigmaq = zeros(1,size(qidMat,2));
    for j = 1:size(qidMat,2)
        muq(j) = mu(j,qidMat(q,j));
        sigmaq(j) = sigma(j,qidMat(q,j));
    end
    sigmaq = sigmaq.^2;
    posMeanE(:,q) = muq' + diag(sigmaq) * A' * ...
        inv(A*diag(sigmaq)*A'+Lambda) * (yt-A*muq');
    posCovEq(:,:,q) = (diag(sigmaq) - diag(sigmaq)*A'*...
        inv(A*diag(sigmaq)*A'+Lambda) * A * diag(sigmaq)) + posMeanE(:,q)*posMeanE(:,q)';
    posCovEq1(:,:,q) = (diag(sigmaq) - diag(sigmaq)*A'*...
        inv(A*diag(sigmaq)*A'+Lambda) * A * diag(sigmaq));
    posMeanE(:,q) = posq(q) * posMeanE(:,q);
end

for i = 1:size(posMeanEqi,1)
    for j = 1:size(posMeanEqi,2)
        ind = qidMat(:,i)==j;
        posMeanEqi(i,j) = sum(ind.*posMeanE(i,:)');
    end
end

posMeanE = sum(posMeanE,2);
for i = 1:size(posCovEqi,1)
    for j = 1:size(posCovEqi,2)
        ind = qidMat(:,i)==j;
        posCovEqi(i,j) = sum(posq(ind)'.*squeeze(posCovEq(i,i,ind)));
    end
end

for q = 1:size(posCovEq,3)
    posCovEq(:,:,q) = posq(q) * posCovEq(:,:,q);
end
posCovE = sum(posCovEq,3);


function [posq, marginalx] = evaluateQ(yt, A, mu, sigma, Dim, numGauss, qidMat, priorq, Lambda)
%  Posterior distribution of latent variable q

% the conditional distribution of q
condProb = conditionalProb(yt, A, mu, sigma, Dim, numGauss, qidMat, Lambda);

% the posterior
marginalx = sum(condProb .* priorq);
posq = condProb .* priorq / marginalx;


function condProb = conditionalProb(yt, A, mu, sigma, Dim, numGauss, qidMat, Lambda)
% Conditional probability p(x_tp1|xtp,q)

condProb = zeros(1, numGauss^Dim);
for q = 1:length(condProb)
    muq = zeros(1,size(qidMat,2));
    sigmaq = zeros(1,size(qidMat,2));
    for j = 1:size(qidMat,2)
        muq(j) = mu(j,qidMat(q,j));
        sigmaq(j) = sigma(j,qidMat(q,j));
    end
    muq = A * muq';
    U = A * diag(sigmaq.^2) * A' + Lambda;
    condProb(q) = gaussian(yt, muq, U);
end

function priorq = priorLatentQ(w, dim, numGauss, numLag, qidMat)
% Calculate the prior distribution of q

% repeat the parameters for all time lags
ws = repmat(w, [numLag,1]);

priorq = zeros(1, numGauss^(numLag*dim));

for q = 1:length(priorq)
    wq = zeros(1,size(qidMat,2));
    for j = 1:size(qidMat,2)
        wq(j) = ws(j,qidMat(q,j));
    end
    priorq(q) = prod(wq);
end

function qidMat = qidMatrix(numGauss, numLag, dim)
% Get all the qs

qidMat = zeros(numGauss^(numLag*dim), numLag*dim);
for i = size(qidMat,2):-1:1
    id = size(qidMat,2) - i;
    qidMatTemp = [];
    for j = 1:numGauss
        qidMatTemp = [qidMatTemp; j*ones(numGauss^id,1)];
    end
    qidMat(:,i) = repmat(qidMatTemp, [numGauss^(numLag*dim-id-1),1]);
end

function prob = gaussian(x, mu, Sigma)
% probability of Gaussian distribution

prob = (2*pi)^(-size(x,1)*0.5) * det(Sigma)^-0.5 * exp(-0.5*(x-mu)'*inv(Sigma)*(x-mu));
