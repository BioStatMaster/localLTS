function [y, W, W_back, sal_back] = sparseica_path(x)
% ICA with SCAD penalized entries of the de-mixing matrix

[N,T] = size(x);
xx = x - mean(x')'*ones(1,T);
% To avoid instability
xx = diag(1./std(xx')) * xx;

% % learning rate
% mu = 1E-6;
% m = 60; % for approximate the derivative of SCAD
% % a = 3.7;
% itmax = 3000;
% iter_M = 2000;
% delta_H = 0;
% Tol = 5e-4;
% w11_back = [];
% w12_back = [];

% initiliazation
fprintf('Initialization....\n')
% WW = diag(1./std(xx'));
% WW = natural_grad(xx, WW);
% save WW_temp.mat WW;
load WW_temp.mat; %WW = WW * diag(std(x'));
% omega = 1./abs(WW);
% % to avoid instability
% Upper = 2 * mean(mean(omega));
% omega = (omega>Upper)*Upper + omega.*(omega<=Upper)
% 
% omega
W = WW;

W_back = W;
sal_back = 0;
Ind_nonzero = 1:N*N;
mu1 = 2E-2;
epoch = 1;
W_col = reshape(W, 1, N*N);
EE = xx * xx'/T;
while 1
    epoch = epoch + 1; fprintf('.');
    W_inv_T = inv(W');
    grad_scaling = -4*diag(diag( W*EE*W' - eye(N) )) * (W*EE);

    y = W*xx;
    y_psi = [];
    for i = 1:N
        tem = estim_beta_pham(y(i,:));
        y_psi =[y_psi; tem(1,:)];
    end

    for i=1:N
        for j=1:N
            %             gx(i,j,:) = repmat(W_inv_T(i,j) + grad_scaling(i,j), 1, T) + y_psi(i,:).*xx(j,:);
            gx(i+(j-1)*N,:) = repmat(W_inv_T(i,j) + grad_scaling(i,j), 1, T) + y_psi(i,:).*xx(j,:);
        end
    end
    % calculate the Hessian matrix

    H = [];
    %     for I = 1:length(Ind_nonzero)
    %         for J = 1:length(Ind_nonzero)
    %             Ind1 = Ind_nonzero(I); Ind2 = Ind_nonzero(J);
    %             H(I,J) = mean(squeeze( gx(mod(Ind1,N)+(~mod(Ind1,N))*N, floor((Ind1-1)/N)+1, :) ).*...
    %                 squeeze( gx(mod(Ind2,N)+(~mod(Ind2,N))*N,  floor((Ind2-1)/N)+1, :) ));
    %         end
    %     end
    H = gx(Ind_nonzero,:) * gx(Ind_nonzero,:)'/T;
    H_inv = inv(H);
    % calculate the saliency and the corresponding lambda
%     W_col = reshape(W, 1, N*N);
    W_col_nonzero = W_col(Ind_nonzero);
    sal = .5 * W_col_nonzero.^2 ./ diag(H_inv)';
    [sal_small, Ind_small] = min(sal * T);
    sal_back = [sal_back sal_small];
    eq = zeros(length(Ind_nonzero),1); eq(Ind_small) = 1;
    W_col_nonzero = W_col_nonzero - W_col_nonzero(Ind_small)/H_inv(Ind_small, Ind_small) * (H_inv * eq)';
    W_col_nonzero = W_col_nonzero([1:Ind_small-1 Ind_small+1:length(Ind_nonzero)]);
%     Ind_temp = find(Ind_nonzero==Ind_small);
    Ind_nonzero = Ind_nonzero([1:Ind_small-1 Ind_small+1:length(Ind_nonzero)] );
    H_new = H([1:Ind_small-1 Ind_small+1:length(Ind_nonzero)],:);
    H_new = H_new(:,[1:Ind_small-1 Ind_small+1:length(Ind_nonzero)]);
    W_col = zeros(1,N*N);
    W_col(Ind_nonzero) = W_col_nonzero;
    W = reshape(W_col, N, N);
    % new W
    for iter = 1:8
        y = W * xx;
        y_psi = [];
        for i = 1:N
            tem = estim_beta_pham(y(i,:));
            y_psi =[y_psi; tem(1,:)];
        end
        W = W + mu1 * ( inv(W') + y_psi*xx'/T );
        W_col_temp = reshape(W, 1, N*N);
        W_col(Ind_nonzero) = W_col_temp(Ind_nonzero);
        W = reshape(W_col, N, N);
    end
    W_back(:,:,epoch) = W;
    % elaborate adjustment
    if sal_small > 20 * log(T)/4 % lambda_IC/2
        break;
    end
end
fprintf('\n');

figure, hold on;
for k=1:N for i=1:N
%         if abs(W(k,i)) < mean(mean(abs(W)))
            plot(sal_back, squeeze(W_back(k,i,:)),'.-');
%         end
    end; end;
set(gca,'XDir','reverse');