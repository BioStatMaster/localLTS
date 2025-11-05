% Experiment: EM algorithm for MoG based ICA, stock
clear;

% generate data
rng('default');
n = 3; % number of sources
D = 2; % number of mixtures
T = 500;

% A = [.5 .0 ; .3 .4];
% w = [0.8 0.2; 0.9, 0.1];
% sigma = [1e-2, 1; 1e-2, 1];
% A = [0.9 0 0; 0.8 0.7 0; 0.9 0.7 0.8];
mu = repmat([0 0],[n,1]);
w = repmat([0.8 0.2],[n 1]);
sigma = repmat([1e-2 1],[n 1]);
% w = [0.8 0.2; 0.8, 0.2; 0.8 0.2];
% sigma = [1e-2 1; 1e-2 1; 1e-2 1];

% random experiments
for exps_id = 1
    rng(exps_id);
    A = rand(D,n)-0.5;
    p = zeros(n, T);
    E = zeros(n, T);
%     for i = 1:n
%         p(i,:) = random('binomial', 1, w(i,2), 1, size(E,2)) + 1;
%         for j = 1:2
%             EE = mu(i,j) + sigma(i,j) * randn(1,sum(p(i,:)==j));
%             E(i,p(i,:)==j) = EE;
%         end
%     end
    
    for i = 1:n
        MU = mu(i,:)';
        SIGMA = reshape(sigma(i,:),1,1,size(w,2));
        P = w(i,:);
        obj = gmdistribution(MU,SIGMA,P);
        E(i,:) = random(obj, T);
    end
    % estimate by MoG ICA
    x = A*E;
        
%     figure, plot(E(1,:), E(2,:), '.');
%     figure, plot(x(1,:), x(2,:), '.');
    
    parsEM.thres = 1e-7;
    parsEM.maxIter = 1000;
    initNoise = 0.8;
    parsEM.A = A+initNoise*rand(D,n)-initNoise/2;
    parsEM.mu = mu;
    parsEM.sigma = sigma;
    parsEM.w = w;
    parsEM.updatePrior = 1;
    parsEM.zeroMean = 0;
    parsEM.noise = 1e-2;
    parsEM.fast = 0;
    tic;
    [A_Hat, sigmaHat, wHat, loglAll, posMeanE] = ifaEM(x, 2, parsEM);
    toc;
%     AA_Hat = repmat(std(posMeanE')./std(E'), D, 1) .* A_Hat;
    A_Hat = A_Hat(:, [1, 3, 2])
    AA_Hat = repmat(sqrt(sum(wHat.*sigmaHat,2)')./sqrt(sum(w.*sigma,2)'), D, 1) .* A_Hat
    mse(AA_Hat, A)
    mse(parsEM.A, A)
    figure(5), plot(loglAll(2:end));
end