function [x_final B_orig W_orig e] = forward_generate(ninputs,T_length,nonlinear_level)
% In our simulation study, ninputs = 8;
% T_length = 1000;
% nonlinear_level = 0.1;
% generating B_orig
B_orig = rand(ninputs,ninputs) - .5;
B_orig = sign(B_orig) .* (.05 + .95*abs(B_orig));
for i=1:ninputs
    for j=i : ninputs
        B_orig(i,j) = 0;
    end
end
% zero entries in the stricy lower triangular part
Ind1 = [4 5 6 6 7 7 8 8 8]; Ind2 = [2 4 3 5 1 2 2 4 6];
for i=1:length(Ind1)
    B_orig(Ind1(i), Ind2(i)) = 0;
end
W_orig = eye(ninputs) - B_orig;
A_orig = inv(W_orig);
% varance of e
var_e = 0.8*rand(ninputs,1) + 0.2;
% generate e
e = zeros(ninputs, T_length);
for i=1:ninputs
    tt = normrnd(0,1,1,T_length);
    expon = rand * .5 + 1.5;
    e(i,:) = sign(tt).*(abs(tt).^expon);
    e(i,:) = e(i,:) - mean(e(i,:));
    e(i,:) = sqrt(var_e(i)) * e(i,:)/std(e(i,:));
end

% generate linear mixtures
x_linear = A_orig * e;
% generate nonlinear mixtures using MLP
[x_nonlinear, MSE, A_fit] = nonlinear_generate(ninputs, e);

% generate x with different levels of nonlinear distortion
x_final = x_linear + sqrt(nonlinear_level) * diag(std(x_linear')./std(x_nonlinear')) * x_nonlinear;