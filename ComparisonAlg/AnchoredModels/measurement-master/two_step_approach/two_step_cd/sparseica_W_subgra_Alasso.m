function [y, W, Score] = sparseica_W_adasize_Alasso(lambda, x)
% ICA with SCAD penalized entries of the de-mixing matrix

[N,T] = size(x);
xx = x - mean(x')'*ones(1,T);
% % To avoid instability
% xx = diag(1./std(xx')) * xx;

% learning rate
mu = 1E-2; % 1E-6
mu_B = 1E-3;
beta = 0;
% m = 60; % for approximate the derivative of SCAD
% a = 3.7;
itmax = 2100;
iter_M = 1400;
delta_H = 0;
Tol = 5e-5;
w11_back = [];
w12_back = [];
W_backup = zeros(N,N,1);
eta_backup = zeros(N,N,1);
z_backup = zeros(N,N,1);
grad_backup = zeros(N,N,1);

% initiliazation
fprintf('Initialization....\n');
 WW = diag(1./std(xx'));
 WW = natural_grad(xx, WW);
 save WW_temp.mat WW;
% load WW_temp.mat; WW = WW * diag(std(x'));
omega = 1./abs(WW);
% to avoid instability
Upper = 4 * mean(mean(omega));
omega = (omega>Upper)*Upper + omega.*(omega<=Upper);

% omega 
W = WW;

z = zeros(N,N);
eta = mu * ones(size(W));
W_old = W + eye(N);

fprintf('Starting penalization...\n');
for iter = 1:itmax
    %     if iter>iter_M
    %         mu = mu*0.99;
    %     end
    fprintf('.');
    if ~mod(iter,100)
        fprintf('\n');
    end
    y = W * xx;
    
    if any(any(isnan(W)))
        pause;
    end

    % normalization
    W = diag(1./std(y')) * W;
    y = diag(1./std(y')) * y;
    
    if sum(sum(abs(W - W_old)))<Tol
        break;
    end
    W_old = W;

    % update W: linear ICA with marginal score function estimated from data...
    y_psi = [];
    for i = 1:N
        tem = estim_beta_pham(y(i,:));
        y_psi =[y_psi; tem(1,:)];
    end

    % dev = omega .* tanh(m*W);
    dev = omega .* sign(W);
%     dev = zeros(size(W));
%     grad_new = y_psi*x'/T + inv(W') - dev*lambda/T;
    % modified objective function
    grad_orig = y_psi*x'/T + inv(W') -4*beta* (diag(diag(y*y'/T)) - eye(N)) * (y*x'/T);
    % obtain the subgradient
      
    for i = 1:N
        for j = 1:N
            if abs(W(i,j)) > 0
                grad_new(i,j) = grad_orig(i,j) - dev(i,j)*lambda/T;
            else if abs(grad_orig(i,j)) > omega(i,j)*lambda/T
                    grad_new(i,j) = sign(grad_orig(i,j)) * (abs(grad_orig(i,j)) - omega(i,j)*lambda/T);
                else
                    grad_new(i,j) = 0;
                end
            end            
        end
    end
    
    % grad_new = y_psi*x'/T + inv(W') -4*beta* (diag(diag(y*y'/T)) - eye(N)) * (y*x'/T) - dev*lambda/T;
    if iter==1
        grad_old = grad_new;
    end
    %     G = y_psi * y'/T - lambda * dev * W';
    %     yy = y*y'/T;
    %     I_N = eye(N);
    %     H = G - diag(diag(G)) + I_N - diag(diag(yy));
    %     W = (I_N + mu*H) * W;
    
    % adaptive size
    [eta, z] = adaptive_size(grad_new, grad_old, eta, z);
    delta_W = eta.*z;
    %         delta_A = mu2 * grad_new;
    W = W + delta_W;
    for i=1:N
        for j=1:N
            if eta(i,j) < mu_B & grad_new(i,j)*grad_old(i,j)<0
                W(i,j) = 0;
            end
        end
    end
    if min(sum(abs(W),2)) == 0
        pause;
    end
%     delta_H = norm(H),
%     if delta_H < Tol
%         break;
%     end
    grad_old = grad_new;
    W_backup(:,:,iter) = W;
    z_backup(:,:,iter) = z;
    eta_backup(:,:,iter) = eta;
    grad_backup(:,:,iter) = grad_new;
%     if sum(sum(abs(delta_W)))<Tol
%         break;
%     end
end

figure, for i=1:4 for j=1:4 subplot(4,4,(i-1)*4 + j), plot(squeeze(W_backup(i,j,:))); end; end
figure, for i=1:4 for j=1:4 subplot(4,4,(i-1)*4 + j), plot(squeeze(eta_backup(i,j,:))); end; end
figure, for i=1:4 for j=1:4 subplot(4,4,(i-1)*4 + j), plot(squeeze(z_backup(i,j,:))); end; end
Score = omega .* abs(W);