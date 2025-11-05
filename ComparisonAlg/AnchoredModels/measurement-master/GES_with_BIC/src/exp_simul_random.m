% Experiment: EM algorithm for constrained ICA
clear;
addpath('../../ifa/src');
addpath('../../GIST_package_1.0/GIST');
addpath(genpath('../../FISTA'));
addpath('../../cICA/src')

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
exp.noise_std = 0;
exp.maxP = 2; % maximum number of parents when searching the graph
% 10 random experiments
for i = 1:50
    X_i = X{i};
    X_i = X_i + exp.noise_std*randn(size(X_i));
    X_i=X_i-repmat(mean(X_i),size(X_i,1),1);
    X_i=X_i*diag(1./std(X_i));
    Record_i = GES(X_i,exp.maxP);
    Bs{i} = Record_i.G;
end
save(fullfile(exp.filePath,sprintf('Bmat_id%d_T%d_std%.1f.mat',dataset.id,dataset.T,exp.noise_std)), 'Bs');
