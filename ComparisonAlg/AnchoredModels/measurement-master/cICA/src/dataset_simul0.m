function [X, B, mu, sigma, w] = dataset_simul0(pars)
% simulation dataset x1->x2<-x3
D = pars.D;
T = pars.T;
n_var = pars.n_var;
id = pars.id;
rng('default');
rng(id);
n = 2 * D;
mu = repmat([0 0],[n,1]);
w = repmat([0.8 0.2],[n 1]);
sigma = repmat([1e-2 1],[n 1]);
sigma(D+1:n,:) = n_var;
% sigma(D+1:n,:) = repmat([.1 0.5],[D 1]);

% B = [0 0 0; 0.3 0 0.4; 0 0 0];
% with cycle
B = [0 0 0 0 0; 0.6 0 0 -0.3 0; 0 0.8 0 0 0; 0 0 -.7 0 0; 0 0.7 0 0 0];
% without cycle
% B = [0 0 0 0 0; 0.6 0 0 0 0; 0 0.8 0 0 0; 0 0 -.7 0 0; 0 0.7 0 0 0];

% A = [1 0 1 0 0;0.3 0.4 0 1 0;0 1 0 0 1];
E = zeros(n, T);
for i = 1:n
    MU = mu(i,:)';
    SIGMA = reshape(sigma(i,:),1,1,size(w,2));
    P = w(i,:);
    obj = gmdistribution(MU,SIGMA,P);
    E(i,:) = random(obj, T);
end

% estimate by MoG ICA
X = (eye(D)-B)\E(1:D,:) + E(D+1:end,:);