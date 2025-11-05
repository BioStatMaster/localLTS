function exp_simulKunCV()
addpath('../ifa');
resDir0 = '../result/simulKunCV/';
dataDir = '../result/simulKun/';
setting = 1;
if ~exist(resDir0, 'dir')
    mkdir(resDir0);
end
load([dataDir,'result_simul' num2str(setting) '.mat']);
rng('default');
if setting == 1 % r = 3, l = 5; chain;
    rng(1);
    aa = sign(rand(1,11)-0.3).*(.7*rand(1,11) + .3);
    B = [0 0 0 0 0 0 0 0; aa(1) 0 0 0 0 0 0 0; 0 aa(2) 0 0 0 0 0 0; 0 0 aa(3) 0 0 0 0 0; aa(4) aa(5) 0 0 0 0 0 0; 0 aa(6) aa(7) 0 0 0 0 0; 0 0 aa(8) aa(9) 0 0 0 0;...
        aa(10) 0 0 aa(11) 0 0 0 0];
elseif setting ==2 %
    B = [0 0 0 0 0 0 0 0; 0.6 0 0 0 0 0 0 0; 0 .7 0 0 0 0 0 0; .7 .4 0 0 0 0 0 0; .5 .4 0 0 0 0 0 0; 0 .4 .7 0 0 0 0 0; 0.6 0 0.5 -.4 0 0 0 0;...
        0.4 -0.4 0.5 0 .6 0 0 0];
end
X = X_all_500(:,:,4);
[rOpt, loglik] = crossValidation(X, N, 5, A_m);
save([resDir0,'result_simul' num2str(setting) '.mat']);
end

function [rOpt loglik] = crossValidation(X, N, folds, A0)
T = size(X,2);
idPerm = randperm(T);
X = X(:,idPerm);
r = 3:6;
loglik = zeros(numel(r),folds);
for i = 1:length(r)
    for j = 1:folds
        foldSize = T/folds;
        X_val = X(:,(j-1)*foldSize+1:j*foldSize);
        X_train = X;
        X_train(:,(j-1)*foldSize+1:j*foldSize) = [];
        loglik(i,j) = loglikelihood(X_train, X_val, r(i), N, A0);
    end
end
loglik = sum(loglik,2);
[~, i] = min(loglik);
rOpt = r(i);
end

function loglik = loglikelihood(X, X_val, r, N, A0)
D = N;
n = D + r;

A = [A0(:,1:r) eye(D)];
mu = repmat([0 0],[n,1]);
w = repmat([0.8 0.2],[n 1]);
sigma = repmat([1e-1 1.5],[n 1]);
w(n-D+1:n,:) = repmat([1 0],[D 1]);
sigma(n-D+1:n,:) = 2;

% experiments
parsEM.dimG = D;
parsEM.nEye = D;
parsEM.thres = 1e-6;
parsEM.minIter = 50;
parsEM.maxIter = 1000;
initNoise = 0;
parsEM.A = A+initNoise*(rand(D,n)-1/2);
parsEM.A(:,end-parsEM.nEye+1:end) = eye(parsEM.nEye);
parsEM.mu = mu;
parsEM.sigma = sigma;
parsEM.w = w;
parsEM.updatePrior = 1;
parsEM.zeroMean = 0;
parsEM.noise = 1e-2;
parsEM.fast = 0;
parsEM.logl = 0;
[A_Hat, sigmaHat, wHat, loglAll, posMeanE] = ifaEM_MN(X, 2, parsEM);

% validation
parsEM.logl = 1;
parsEM.A = A_Hat;
parsEM.sigma = sigmaHat;
parsEM.w = wHat;
[A_Hat, sigmaHat, wHat, loglAll, posMeanE] = ifaEM_MN(X_val, 2, parsEM);
loglik = loglAll(end);
end
