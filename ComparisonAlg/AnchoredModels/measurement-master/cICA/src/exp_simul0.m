% Experiment: EM algorithm for constrained ICA
clear;
addpath('../../ifa/src');
addpath('../../GIST_package_1.0/GIST');
addpath(genpath('../../FISTA'));

% load data
dataset.D = 5;
dataset.T = 500;
dataset.n_var = 0.1;
dataset.id = 1;
dataset.load_func = @dataset_simul0;
[X, B, mu, sigma, w] = dataset.load_func(dataset);

% set EM parameters
D = dataset.D;
parsEM.refine = false;
parsEM.dimG = D;
parsEM.nEye = D;
parsEM.thres = 1e-6;
parsEM.minIter = 500;
parsEM.maxIter = 1000;
initNoise = 0.8;
parsEM.B = B+initNoise*(rand(D,D)-1/2);
parsEM.mu = mu;
parsEM.sigma = sigma;
parsEM.sigma = [sigma(1:D,:);repmat(sqrt(sum(w(D+1:end,:).*sigma(D+1:end,:).^2,2)),[1 2])];
% parsEM.sigma = [sigma(1:D,:);repmat(sqrt(sum(w(D+1:end,:).*sigma(D+1:end,:).^2,2)),[1 2])].*(1+2*(rand(n,2)-0.5));
parsEM.w = w;
parsEM.updatePrior = 1;
parsEM.zeroMean = 0;
parsEM.noise = 1e-2;
parsEM.fast = 0;
parsEM.logl = 0;

% start training
parsEM.lambda = 0;
[B_Hat, sigmaHat, wHat, loglAll, posMeanE] = ifaEM_MN_Sparse(X, 2, parsEM);
parsEM.B = B_Hat;

parsEM.lambda = 1/2*log(size(X,2));
[B_Hat, sigmaHat, wHat, loglAll, posMeanE] = ifaEM_MN_Sparse(X, 2, parsEM);

% plot(loglAll);

% refie results
% B_Hat(B_Hat<0.1) = 0;
% parsEM.refine = true;
% parsEM.mask = B_Hat~=0;
% [B_HatRef, sigmaHat, wHat, loglAll, posMeanE] = ifaEM_MN_Sparse(X, 2, parsEM);

% display results
B
B_Hat
% B_HatRef
mse(B,B_Hat)
% mse(B,B_HatRef)