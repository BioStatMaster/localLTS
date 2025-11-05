% Experiment: EM algorithm for MoG based ICA, stock
clear;
addpath('../ifa');

% generate data
rng('default');
r = 2;
D = 3;
n = r + D;
T = 500;

n0 = 2 * D;
mu = repmat([0 0],[n0,1]);
w = repmat([0.8 0.2],[n0 1]);
sigma = repmat([1e-2 1],[n0 1]);
% w(D+1:n0,:) = repmat([1 0],[D 1]);
% sigma(D+1:n0,:) = 0.5;
sigma(D+1:n0,:) = repmat([.1 0.5],[D 1]);
% experiments
rng(1);
B = [0 0 0; 0.3 0 0.4; 0 0 0];
A = [1 0 1 0 0;0.3 0.4 0 1 0;0 1 0 0 1];
p = zeros(n0, T);
E = zeros(n0, T);
for i = 1:n0
    p(i,:) = random('binomial', 1, w(i,2), 1, size(E,2)) + 1;
    for j = 1:2
        EE = mu(i,j) + sigma(i,j) * randn(1,sum(p(i,:)==j));
        E(i,p(i,:)==j) = EE;
    end
end

% estimate by MoG ICA
x = (eye(D)-B)\E(1:D,:) + E(D+1:end,:);

%     figure, plot(E(1,:), E(2,:), '.');
%     figure, plot(x(1,:), x(2,:), '.');
parsEM.dimG = D;
parsEM.nEye = D;
parsEM.thres = 1e-6;
parsEM.maxIter = 1000;
initNoise = 0.8;
parsEM.A = A+initNoise*(rand(D,n)-1/2);
parsEM.A(:,end-parsEM.nEye+1:end) = eye(parsEM.nEye);
parsEM.mu = [mu(1:r,:);mu(D+1:end,:)];
parsEM.sigma = [sigma(1:r,:);repmat(sqrt(sum(w(D+1:end,:).*sigma(D+1:end,:).^2,2)),[1 2])];
% parsEM.sigma = [sigma(1:r,:);repmat(sqrt(sum(w(D+1:end,:).*sigma(D+1:end,:).^2,2)),[1 2])].*(1+2*(rand(n,2)-0.5));
parsEM.w = [w(1:r,:);w(D+1:end,:)];
parsEM.updatePrior = 1;
parsEM.zeroMean = 0;
parsEM.noise = 1e-2;
parsEM.fast = 0;
tic;
[A_Hat, sigmaHat, wHat, loglAll, posMeanE] = ifaEM_MN(x, 2, parsEM);
toc;
scale = sqrt(sum(wHat(1:r,:).*sigmaHat(1:r,:).^2,2))./sqrt(sum(w(1:r,:).*sigma(1:r,:).^2,2));
A_Hat(:,1:r) = A_Hat(:,1:r).*repmat(scale',[D 1]);
mse(A_Hat(:,1:r), A(:,1:r))
mse(parsEM.A(:,1:r), A(:,1:r))
figure(5), plot(loglAll(2:end));
