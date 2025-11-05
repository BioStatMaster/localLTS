function mse_loss = loss_mse(B, X, posMeanE, posCovE, Lambda)
% EM negative lower bound
Lambda = inv(Lambda);
d = size(B,1)^0.5;
n = size(X,2);
B = reshape(B,d,d);
A = [inv(eye(d)-B), eye(d)];
mse_loss = 0.5*trace(A'*Lambda*A*sum(posCovE,3)) - trace(A'*Lambda*X*posMeanE');