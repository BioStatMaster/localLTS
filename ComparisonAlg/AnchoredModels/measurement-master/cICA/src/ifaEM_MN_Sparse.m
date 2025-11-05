function [B, sigma, w, loglAll, posMeanE] = ifaEM_MN_Sparse(X, numGauss, parsEM)

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
B = parsEM.B;
B0 = B;
B(logical(eye(dim))) = diag(zeros(dim));
mu = parsEM.mu;
sigma = parsEM.sigma;
w = parsEM.w;
dimG = parsEM.dimG;
nEye = parsEM.nEye;
Dim = size(w,1);
Lambda = parsEM.noise*eye(dim);
lambda = parsEM.lambda;

iter = 0;
loglAll = [];
logl_prev = inf;
logl = inf;
while iter < parsEM.minIter || iter < parsEM.maxIter && abs(logl_prev-logl) >= abs(logl_prev)*thres
    % M step
    % update A
    logl_prev = logl;
    if iter > 0
        % update B, A=(I-B)^-1
%         A0 = (X*posMeanE')/(sum(posCovE,3));
%         B0 = eye(dimG) - inv(A0(1:dimG,1:dimG));
%         B = B0;
%         B = updateB_SCAD(B, X, posMeanE, posCovE, Lambda, B0);
        if parsEM.refine
            B = B0;
            B = B.*parsEM.mask;
        else
            B = updateB_ALasso(B, X, posMeanE, posCovE, Lambda, B0, lambda);
        end
        
        B(1:dimG+1:end) = zeros(1,dimG);
        if parsEM.updatePrior == 1
            posqi = marginalPosQ(posq, Dim, dimG, numGauss, qidMat);
            % update mu
            if mu~=zeros(Dim,numGauss)
                mu = updateMu(posqi, posMeanEqi);
            end
            % update sigma
            sigma = updateSigma(mu, posqi, posCovEqi);
            
            % update w
            w = updateW(posqi);
        end
    end
    
    % E step
    qidMat = qidMatrix(numGauss, Dim, dimG);
    priorq = priorLatentQ(w, Dim, dimG, numGauss, qidMat);
    
    % posterior of q and e
    A = [inv(eye(dimG)-B),eye(nEye)];
    [posq, marginalx] = evaluateQ(X, A, mu, sigma, Dim, dimG, numGauss, qidMat, priorq,Lambda);
    [posMeanE, posCovE, posMeanEqi, posCovEqi] = evaluateE(X, posq, A, Lambda, ...
        mu, sigma, Dim, dimG, numGauss, qidMat);
    logl = -sum(log(marginalx));
    fprintf('iter%d: negative loglik: %.10f\n', iter, logl);
    B
    mu
    sigma
    w
    iter = iter + 1;
    loglAll = [loglAll, logl];
    if parsEM.logl == 1
        return;
    end
end
