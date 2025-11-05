function [y, W] = sparseica_SCAD_FUN(lambda_in, x, Init,WW)
% ICA with SCAD penalized entries of the de-mixing matrix
% If you have initialization value for W, let Init = 1, and WW is the
% initialization value; otherwise use 0 for Init.
[N,T] = size(x);
xx = x - mean(x')'*ones(1,T);

% learning rate
mu = 0.01;
m = 180; % for approximate the derivative of SCAD
a = 3.7; 
itmax = 1500;
iter_M = 1400;
delta_H = 0;
Tol = 5e-4;
w11_back = [];
w12_back = [];
Regular = 1,
Init_fastica = 1;

if ~Init
    if Init_fastica
        figure(7), [ica_fastica, AA, WW] = fastica(xx, 'approach', 'symm', 'g', 'tanh');
    else
        WW = diag(1./std(xx'));
    end
end

% initialization
W = WW;

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
    
    dev = zeros(N,N);
    if Regular & (Init | iter > 50)
        % calculates the derivatives of the SCAD term
        for row=1:N
            for column = 1:N
                if abs(W(row, column)) <= lambda_in
                    dev(row, column) = tanh(m*W(row, column));
                else if abs(W(row, column)) <= a * lambda_in
                        dev(row, column) = sign(W(row, column))*tanh(m*lambda_in)*(a*lambda_in - abs(W(row, column)))/(a-1)/lambda_in;
                    end
                end
            end
        end
        G = y_psi * y'/T - lambda_in * dev * W'; 
    else
        G = y_psi * y'/T;    
    end
    yy = y*y'/T;  
    I_N = eye(N);
    H = G - diag(diag(G)) + I_N - diag(diag(yy));
    W = (I_N + mu*H) * W;  
    delta_H = norm(H);
    if delta_H < Tol 
        break;
    end
end