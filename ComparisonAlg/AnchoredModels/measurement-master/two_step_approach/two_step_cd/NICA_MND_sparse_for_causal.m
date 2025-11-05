% to generate nonlinear mixtures with a given nonlinear distortion level
% and with the linear part following LiNGAM, and then to identify their 
% linear causal relations

% nonlinear distortion level in the mixing system
Non_level = 0.01; % 0.03?
% dimentionality & sample size
n = 8; T = 1000;
[x_final B_orig W_orig e_s] = forward_generate(n,T,Non_level);

trpattern = x_final;

netpar_MND;
% value of lambda, lambda_in, and a
lambda = 0.05 / (sum(diag(trpattern*trpattern'/ntrain))/ninputs);
lambda_in = 0.06; a = 3.7; % 0.10?
IF_direct = 1;
Linear_init = 1;
nepochs = 800;

% net initialization
netinit_MND;
% network training & testing for causality
train_MND_sparse;

if IF_causality
    fprintf('Linear causal relations successfully identified. Causal relations implied by B = \n');
    disp(Best);
else
    fprintf('Linear causal relations identification: failed...\n');
end