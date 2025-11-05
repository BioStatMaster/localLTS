function B = updateB_ALasso(B, X, posMeanE, posCovE, Lambda, B0, lambda)
% update B, with sparse penalty

% input parameters
% lambda = 1e-1*abs(randn);
% lambda = 1/2*log(size(X,2));
% lambda = 0;
B = B(:);
B0 = B0(:);
lambda = lambda./abs(B0+eps);
opts.lambda = lambda;
B = fista_backtracking(@(x)loss_mse(x, X, posMeanE, posCovE, Lambda), ...
    @(x)loss_mse_grad(x, X, posMeanE, posCovE, Lambda), B, opts, ...
    @(x)loss_mse_L1(x, X, posMeanE, posCovE, Lambda, opts.lambda));  
B = reshape(B, sqrt(length(B)), sqrt(length(B)));


    
