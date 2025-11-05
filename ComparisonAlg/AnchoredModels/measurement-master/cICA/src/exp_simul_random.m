% Experiment: EM algorithm for constrained ICA
clear;
addpath('../../ifa/src');
addpath('../../GIST_package_1.0/GIST');
addpath(genpath('../../FISTA'));

% load data
dataset.filePath = '../../datasets/simul_random';
dataset.id = 0;
dataset.T = 500;
dataset.load_func = @dataset_simul_random;
[X, B, G] = dataset.load_func(dataset);

exp.filePath = '../outputs/simul_random';
if ~exist(exp.filePath, 'dir')
    mkdir(exp.filePath);
end
exp.noise_std = 0.4;

% 10 random experiments
for i = 1:50
    X_i = X{i}';
    B_i = B{i}';
    X_i = X_i + exp.noise_std*randn(size(X_i));
    % set EM parameters
    D = size(X_i,1);
    parsEM.refine = false;
    parsEM.dimG = D;
    parsEM.nEye = D;
    parsEM.thres = 1e-6;
    parsEM.minIter = 500;
    parsEM.maxIter = 1000;
    initNoise = 0.8;
    parsEM.B = B_i+initNoise*(rand(D,D)-1/2);
    parsEM.mu = repmat([0 0],[2*D,1]);
    parsEM.w = repmat([0.8 0.2],[2*D 1]);
    parsEM.sigma = repmat([1e-2 1],[2*D 1]);
    parsEM.sigma(D+1:2*D,:) = 0.1;
    parsEM.updatePrior = 1;
    parsEM.zeroMean = 0;
    parsEM.noise = 1e-2;
    parsEM.fast = 0;
    parsEM.logl = 0;
    
    % start training
    [B_Hat_i, sigmaHat_i, wHat_i, loglAll_i, posMeanE_i] = ifaEM_MN_Sparse(X_i, 2, parsEM);
    % plot(loglAll);
    
    % refie results
%     B_Hat_i(B_Hat_i<0.1) = 0;
    parsEM.refine = true;
    parsEM.mask = B_Hat_i~=0;
    [B_HatRef_i, sigmaHat_i, wHat_i, loglAll_i, posMeanE_i] = ifaEM_MN_Sparse(X_i, 2, parsEM);
    
%     % display results
%     B_i
%     B_Hat_i
%     B_HatRef_i
%     mse(B_i,B_Hat_i)
%     mse(B_i,B_HatRef_i)

    Bs{i} = B_HatRef_i;
end
save(fullfile(exp.filePath,sprintf('Bmat_id%d_T%d_std%.1f.mat',dataset.id,dataset.T,exp.noise_std)), 'Bs');
