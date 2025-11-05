function mse_loss = loss_mse_L1(B, X, posMeanE, posCovE, Lambda, lambda)
% EM negative lower bound + L1 norm
mse_loss = loss_mse(B, X, posMeanE, posCovE, Lambda);
mse_loss = mse_loss + sum(lambda.*abs(B(:)));