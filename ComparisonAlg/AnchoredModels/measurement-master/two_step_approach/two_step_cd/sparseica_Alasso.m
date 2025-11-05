function [y, W, Score] = sparseica_Alasso(lambda, x)
% ICA with SCAD penalized entries of the de-mixing matrix

[N,T] = size(x);
xx = x - mean(x')'*ones(1,T);
% % To avoid instability
% xx = diag(1./std(xx')) * xx;

% learning rate
mu = 1E-6;
m = 60; % for approximate the derivative of SCAD
% a = 3.7;
itmax = 3000;
iter_M = 2000;
delta_H = 0;
Tol = 5e-4;
w11_back = [];
w12_back = [];

% initiliazation
fprintf('Initialization....\n')
 WW = diag(1./std(xx'));
 WW = natural_grad(xx, WW);
 save WW_temp.mat WW;
% load WW_temp.mat; WW = WW * diag(std(x'));
omega = 1./abs(WW);
% to avoid instability
Upper = 2 * mean(mean(omega));
omega = (omega>Upper)*Upper + omega.*(omega<=Upper)

omega 
W = WW;

fprintf('Starting penalization...\n');
for iter = 1:itmax
    if iter>iter_M
        mu = mu*0.99;
    end
    fprintf('.');
    if ~mod(iter,100)
        fprintf('\n');
    end
    y = W * xx;

    % update W: linear ICA with marginal score function estimated from data...
    y_psi = [];
    for i = 1:N
        tem = estim_beta_pham(y(i,:));
        y_psi =[y_psi; tem(1,:)];
    end

    dev = omega .* tanh(m*W);
    G = y_psi * y'/T - lambda * dev * W';
    yy = y*y'/T;
    I_N = eye(N);
    H = G - diag(diag(G)) + I_N - diag(diag(yy));
    W = (I_N + mu*H) * W;
    delta_H = norm(H),
    if delta_H < Tol
        break;
    end
end
Score = omega .* abs(W);