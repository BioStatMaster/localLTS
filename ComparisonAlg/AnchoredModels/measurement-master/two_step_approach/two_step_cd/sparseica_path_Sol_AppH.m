function [y, W, W_back, sal_back] = sparseica_path_Sol(x)
% function [y, W, W_back, sal_back] = sparseica_path(x)
% ICA with sparse entries of the de-mixing matrix:
% INPUT: 
%   x: the data (each row corresponds to a variable);
% OUTPUTS:
%   y: the final separated signals;
%   W: the final separation matrix;
%   W_back: contains all separation matrix with different number of non-zero
%       entries (if you want W with k zero entries, you should use W_back(:,:,k-1));
%   sal_back: the saliency of each removed entry.
% NOTE: Cri_stop in line 35 determines when to stop the pruning process.
% The default value is 4*\lambda_{BIC}.

[N,T] = size(x);
xx = x - mean(x')'*ones(1,T);
% To avoid instability
xx = diag(1./std(xx')) * xx;

% initiliazation
fprintf('Initialization....\n')
% WW = diag(1./std(xx'));
% WW = natural_grad(xx, WW);
% save WW_temp.mat WW;
load WW_temp_RET_NORM.mat; %WW = WW * diag(std(x'));
% load WW_temp.mat; %WW = WW * diag(std(x'));
W = WW;
fprintf('\n');

% Settings
W_back = W;
sal_back = 0;
Ind_nonzero = 1:N*N;
H_change_norm = [];
H_new = zeros(N*N, N*N);
Cri_stop = 4 * log(T)/4;
mu1 = 0.02;
epoch = 1;
W_col = reshape(W, 1, N*N);
EE = xx * xx'/T;

% find the solution path
while 1
    epoch = epoch + 1; fprintf('*removing the %d-th entry\n', epoch-1);
    % for finding the Hessian matrix
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
            gx(i+(j-1)*N,:) = repmat(W_inv_T(i,j) + grad_scaling(i,j), 1, T) + y_psi(i,:).*xx(j,:);
        end
    end
    
    % calculate the Hessian matrix
    H = [];
    H = gx(Ind_nonzero,:) * gx(Ind_nonzero,:)'/T;
    H_change_norm = [ H_change_norm norm(H-H_new)];
    % To reduce the complexity, we use approximations
    current_length = length(Ind_nonzero);
    if mod(epoch,5) == 2
        H_inv = inv(H);
    else
        Temp = (H-H_new)*H_new_inv;
        H_inv = H_new_inv * ( eye(current_length) - Temp *(eye(current_length) - Temp));
        H_inv1 = inv(H); norm(H_inv1 - H_inv),
    end
%     H_inv = inv(H);

    % calculate the saliency and the corresponding lambda, and prune
    %     W_col = reshape(W, 1, N*N);
    W_col_nonzero = W_col(Ind_nonzero);
    sal = .5 * W_col_nonzero.^2 ./ diag(H_inv)';
    [sal_small, Ind_small] = min(sal * T);
    if sal_small > Cri_stop % lambda_IC/2
        break;
    else
        sal_back = [sal_back sal_small];
    end

    % update the new weights accoding to the current Hessian matrix
    eq = zeros(current_length,1); eq(Ind_small) = 1;
    W_col_nonzero = W_col_nonzero - W_col_nonzero(Ind_small)/H_inv(Ind_small, Ind_small) * (H_inv * eq)';
    W_col_nonzero = W_col_nonzero([1:Ind_small-1 Ind_small+1:current_length]);
    %     Ind_temp = find(Ind_nonzero==Ind_small);
    Ind_nonzero = Ind_nonzero([1:Ind_small-1 Ind_small+1:current_length] );
    H_new = H([1:Ind_small-1 Ind_small+1:current_length],:);
    H_new = H_new(:,[1:Ind_small-1 Ind_small+1:current_length]);
    H_new_inv = inv(H_new); %%% To find a easy and efficient way to calculate H_new_inv...
    W_col = zeros(1,N*N);
    W_col(Ind_nonzero) = W_col_nonzero;
    W = reshape(W_col, N, N);
    % elaborate adjustment of the separation matrix
    for iter = 1: ( 30 + ~mod(epoch,5)*150)
% for iter = 1: 1
        fprintf('.');
        y = W * xx;
        y_psi = [];
        for i = 1:N
            tem = estim_beta_pham(y(i,:));
            y_psi =[y_psi; tem(1,:)];
        end

        grad_scaling = -4*diag(diag( W*EE*W' - eye(N) )) * (W*EE);

        %         grad_temp = reshape(inv(W') + y_psi*xx'/T, 1, N*N);
        grad_temp = reshape(inv(W') + y_psi*xx'/T + grad_scaling, 1, N*N);
        grad_temp = grad_temp(Ind_nonzero);
        %                 W_col_nonzero = W_col_nonzero + .5 * grad_temp * H_new_inv;
%         W_col_nonzero = W_col_nonzero + mu1 * grad_temp;
        W_col_nonzero = W_col_nonzero + .1 * grad_temp * H_new_inv';

        scal_grad = sum(sum(abs(grad_temp ))); fprintf('%f',scal_grad);
        % check if further iterations are needed
        %         if norm(grad_temp)< 7E-4*N
        if sum(sum(abs(scal_grad))) < 1E-4 * current_length
            break;
        end

        W_col(Ind_nonzero) = W_col_nonzero;
        W = reshape(W_col, N, N);
    end
    W_back(:,:,epoch) = W;
end
fprintf('\n');

% plot the results
figure, hold on;
for k=1:N for i=1:N
        plot(sal_back, squeeze(W_back(k,i,:)),'.-');
    end; end;
xlabel('\lambda = 1/2*\beta_q^2/[H^{-1}]_{q,q}');
ylabel('Weight'); title('Evolution of all weights');
set(gca,'XDir','reverse'); 
line([log(T)/4, log(T)/4], [min(min(W))-.3, max(max(W))+.3],...
    'Color', 'red', 'LineWidth',2,'MarkerFaceColor','r'); text(log(T)/4-.02, min(min(W))-.2, 'BIC')

% zero weights
figure, hold on;
for k=1:N for i=1:N
        if ~abs(W(k,i))
            plot(sal_back, squeeze(W_back(k,i,:)),'.-');
        end
    end; end;
xlabel('\lambda = 1/2*\beta_q^2/[H^{-1}]_{q,q}');
ylabel('Weight'); title('Evolution of the weights set to zero');
set(gca,'XDir','reverse');
% line([log(T)/4, log(T)/4], [min(min(W))-.3, max(max(W))+.3],...
%     'Color', 'red', 'LineWidth',2,'MarkerFaceColor','r'); text(log(T)/4-.02, min(min(W))-.2, 'BIC')

% figure, plot(H_change_norm);