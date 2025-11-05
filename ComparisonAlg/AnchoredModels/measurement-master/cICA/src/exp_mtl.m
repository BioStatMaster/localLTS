% Experiment: EM algorithm for MoG based ICA, stock
clear;
addpath('../ifa');
x = [];
for i = 14:15
fileName = sprintf('../data/hippocampus/sub0%d_mtlmean.txt', i);
x = [x;load(fileName)];
end
X = x';
X = X(2:7,:);
% generate data
rng('default');
r = 5;
D = size(X,1);
n = D + r;

A = rand(D,n);
A(:,r+1:end) = eye(D);
mu = repmat([0 0],[n,1]);
w = repmat([0.8 0.2],[n 1]);
sigma = repmat([1e-1 2],[n 1]);
w(n-D+1:n,:) = repmat([1 0],[D 1]);
sigma(n-D+1:n,:) = 1;

% experiments
parsEM.dimG = D;
parsEM.nEye = D;
parsEM.thres = 1e-7;
parsEM.minIter = 5000;
parsEM.maxIter = 10000;
parsEM.A = A;
parsEM.mu = mu;
parsEM.sigma = sigma;
parsEM.w = w;
parsEM.updatePrior = 1;
parsEM.zeroMean = 0;
parsEM.noise = 1e-2;
parsEM.fast = 0;
tic;
[A_Hat, sigmaHat, wHat, loglAll, posMeanE] = ifaEM_MN(X, 2, parsEM);
toc;
%figure(5), plot(loglAll(2:end));
save(sprintf('mtl%d.mat',r)); 
exit
