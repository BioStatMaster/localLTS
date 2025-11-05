%given x

% lambda = 0.04; %.15
[N,T] = size(x);
xx = x - mean(x')'*ones(1,T);

%%%%%%%%%%%%%%%%%%%%%%%%%
xx = xx/mean(std(x'));

lambda = log(T)/2;
T1 = 200;
delta_A = ones(N,N);
mu2 = 4E-2;
mu1 = 0.03;
% m = 1000; %8
% a = 3.4;
itmax = 3000;
iter_M = 1700;
delta_H = 0;
Tol = 5e-5;
Th = 0.09/mean(std(x'));
% lambda_in = 1.6 * lambda;
A_backup = zeros(N,N,itmax);
eta_backup = zeros(N,N,itmax);
z_backup = zeros(N,N,itmax);
grad_backup = zeros(N,N,itmax);
% Regular = 1;
Init_fastica = 1;
%%%%%%%%%%%%%%%%%%%%%%%%%%%% since A has been obtained
% if Init_fastica
%     figure, [ica_fastica, AA, WW] = fastica(xx, 'g', 'tanh'); %'approach', 'symm',
% end
% % % else
% %     AA = diag(std(xx'));
% % % end
% % %
% % % initialization
 load A_temp_fin.mat;
 AA= A;
A = AA;


% % the weights
% omega = 1./abs(AA);
% % to avoid instability
% Upper = 1 * mean(mean(omega)); % 2
% omega = (omega>Upper)*Upper + omega.*(omega<=Upper)
% omega,

for iter = 1:itmax
    if iter>iter_M
        mu = mu*0.99;
        if ~mod(iter,10)
            %             m = m + 3;
        end
    end
    fprintf('.');
    if ~mod(iter,100)
        fprintf('\n');
    end
    y = inv(A) * xx;

    % normalization
    A = A * diag(std(y'));
    y = diag(1./std(y')) * y;

    % update W: linear ICA with marginal score function...
    y_psi = []; %%%
    for i = 1:N %%%
        tem = estim_beta_pham(y(i,:)); %%%
        y_psi =[y_psi; tem(1,:)]; %%%
    end %%%

    %     dev = zeros(N,N);
    %     if Regular
    %         for row=1:N
    %             for column = 1:N
    %                 if abs(A(row, column)) <= lambda_in
    %                     dev(row, column) = tanh(m*A(row, column));
    %                 else if abs(A(row, column)) <= a * lambda_in
    %                         dev(row, column) = sign(A(row, column))*tanh(m*lambda_in)*(a*lambda_in - abs(A(row, column)))/(a-1)/lambda_in;
    %                     end
    %                 end
    %             end
    %         end
    if iter > T1
        %         dev = omega .* tanh(m*A);
        dev = omega .* sign(A);
        %         G = y_psi * y'/T + lambda * A' *dev; %%%
        grad_old = grad_new;
        %         grad_new = inv(A') * (-y_psi*y'/T - eye(N) ) - dev * lambda/T;
        grad_new = A * ( (-y_psi*y'/T - eye(N) ) - A' * dev * lambda/T);
        %         mu = mu2;
        %         %         G = y_psi * y'/T + lamda * A' * tanh(m*A); %%%
    else
        G = y_psi * y'/T;
        mu = mu1;
        z = zeros(N,N);
        eta = mu2 * ones(size(A));
        if iter == T1
            omega = 1./abs(A);
            % to avoid instability
            Upper = 1.5 * mean(mean(omega)); % 2
            omega = (omega>Upper)*Upper + omega.*(omega<=Upper)
            omega,            %             dev = omega .* tanh(m*A);
            dev = omega .* sign(A);
            %             grad_new = inv(A') * (-y_psi*y'/T - eye(N) ) - dev * lambda/T;
            grad_new = A * ( (-y_psi*y'/T - eye(N) ) - A' * dev * lambda/T);
            AA = A; save AA_temp.mat AA;
            % the weights

        else
            %             grad_new = inv(A') * (-y_psi*y'/T - eye(N) );
            grad_new = A * (-y_psi*y'/T - eye(N) ) ;
        end
    end
    yy = y*y'/T;
    I_N = eye(N);
    % H = G - diag(diag(G)) + I_N - diag(diag(UU)); %%%
    if iter <= T1
        H = G - diag(diag(G)) + I_N - diag(diag(yy));
        A = A * (I_N - mu*H);  %%%
        delta_H = norm(H);
    else
        [eta, z] = adaptive_size(grad_new, grad_old, eta, z);
        delta_A = eta.*z;
        %         delta_A = mu2 * grad_new;
        A = A + delta_A;
    end
    A_backup(:,:,iter) = A;
    z_backup(:,:,iter) = z;
    eta_backup(:,:,iter) = eta;
    grad_backup(:,:,iter) = grad_new;
    if delta_H < Tol | sum(sum(abs(delta_A)))<Tol
        break;
    end
end

% if Regular
%     A = (abs(A)>Th).*A;
% end

fprintf('iter = %d, delta_H = %6.4d\n', iter, delta_H);
% % A = A * diag(std(y')) * mean(std(x')),

figure, for i=1:4 for j=1:10 subplot(4,10,(i-1)*10 + j), plot(ss(i,:), icasig(j,:), '.'); end; end; title('icasig')
figure, for i=1:4 for j=1:10 subplot(4,10,(i-1)*10 + j), plot(ss(i,:), y(j,:), '.'); end; end; title('y')