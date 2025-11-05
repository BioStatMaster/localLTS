% Experiment: EM algorithm for MoG based ICA, stock
clear;
addpath('ifa');
fileName = 'data/simul_kun/data_simul3.mat';
load(fileName);
% generate data
rng('default');
r = sum(sum(B)~=0);
D = size(X,1);
n = D + r;

A = inv(eye(D)-B);
A = [A(:,1:r) eye(D)];
mu = repmat([0 0],[n,1]);
w = repmat([0.8 0.2],[n 1]);
sigma = repmat([1e-1 1.5],[n 1]);
w(n-D+1:n,:) = repmat([1 0],[D 1]);
sigma(n-D+1:n,:) = 2;

% experiments
parsEM.dimG = D;
parsEM.nEye = D;
parsEM.thres = 1e-6;
parsEM.minIter = 500;
parsEM.maxIter = 5000;
initNoise = 0.8;
parsEM.A = A+initNoise*(rand(D,n)-1/2);
parsEM.A(:,end-parsEM.nEye+1:end) = eye(parsEM.nEye);
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
scale = sqrt(sum(wHat(1:r,:).*sigmaHat(1:r,:).^2,2))./std(Err(1:r,:),0,2);
A_Hat(:,1:r) = A_Hat(:,1:r).*repmat(scale',[D 1]);
mse0 = mse(parsEM.A(:,1:r), A(:,1:r));
mse1 = mse(A_Hat(:,1:r), A(:,1:r));
fprintf('initial MSE is %f, and estimated MSE is %f.\n', mse0, mse1);
figure(5), plot(loglAll(2:end));
save(fileName, '-append');
