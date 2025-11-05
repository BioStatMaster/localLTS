% Experiment: EM algorithm for constrained ICA
clear;

% load data

% set EM parameters
parsEM.dimG = D;
parsEM.nEye = D;
parsEM.thres = 1e-6;
parsEM.minIter = 500;
parsEM.maxIter = 1000;
initNoise = 0.8;
parsEM.B = B+initNoise*(rand(D,D)-1/2);
parsEM.mu = [mu(1:r,:);mu(D+1:end,:)];
parsEM.sigma = [sigma(1:r,:);repmat(sqrt(sum(w(D+1:end,:).*sigma(D+1:end,:).^2,2)),[1 2])];
% parsEM.sigma = [sigma(1:r,:);repmat(sqrt(sum(w(D+1:end,:).*sigma(D+1:end,:).^2,2)),[1 2])].*(1+2*(rand(n,2)-0.5));
parsEM.w = [w(1:r,:);w(D+1:end,:)];
parsEM.updatePrior = 1;
parsEM.zeroMean = 0;
parsEM.noise = 1e-2;
parsEM.fast = 0;
parsEM.logl = 0;
tic;
[B_Hat, sigmaHat, wHat, loglAll, posMeanE] = runIFA_MN_Sparse(x, 2, parsEM);
toc;
plot(loglAll);