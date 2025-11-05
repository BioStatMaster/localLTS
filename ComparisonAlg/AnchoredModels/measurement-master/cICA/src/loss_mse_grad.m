function grad_mse_loss = loss_mse_grad(B, X, posMeanE, posCovE, Lambda)
% EM negative lower bound gradient
d = size(B,1)^0.5;
B = reshape(B,d,d);
Lambda = inv(Lambda);
A = [inv(eye(d)-B), eye(d)];
A1 = A(1:d,1:d);
grad1 = Lambda * A * sum(posCovE,3);
grad2 = Lambda * X * posMeanE';
grad_loss = A1'*(grad1(1:d,1:d)-grad2(1:d,1:d))*A1';
grad_mse_loss = grad_loss(:);
grad_mse_loss(1:d+1:end) = zeros(1,d);