% to generate nonlinear mixtures with a given nonlinear distortion level
% and with the linear part following LiNGAM, and then to identify their
% linear causal relations

% nonlinear distortion level in the mixing system
Non_level = 0.02;
lambda_in = 0.2;
% dimentionality & sample size
n = 8; T = 1000;
[x_final B_orig W_orig e_s] = forward_generate(n,T,Non_level);

% Run ICA with the SCAD penalty on each entry of the de-mixing matrix
% [y, W, Score] = sparseica_Alasso(1.2 * log(T)/2, x_final);
[y, W, Score] = sparseica_Alasso(log(T), x_final);
%check causality
[Wp,rowp] = nzdiaglinprog( W );
fprintf('Done!\n');
% threshold
Tol = 0.01;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Divide each row of Wp by the diagonal element
estdisturbancestd = 1./diag(abs(Wp));
Wp = Wp./(diag(Wp)*ones(1,ninputs)),
Wp = (abs(Wp)>Tol) .* Wp;

% Compute corresponding B
Best = eye(ninputs)-Wp;

% Check if B is lower-triangular after permutation
[IF_causality, Bestcausal,causalperm] = check_causality( Best );
if IF_causality
    fprintf('Linear causal relations successfully identified. Causal relations implied by B = \n');
    disp(Best);
    icausal = iperm(causalperm);
    B = Bestcausal(icausal, icausal);
    k = causalperm;
else
    fprintf('Linear causal relations identification with linear ICA: failed...\n');
end