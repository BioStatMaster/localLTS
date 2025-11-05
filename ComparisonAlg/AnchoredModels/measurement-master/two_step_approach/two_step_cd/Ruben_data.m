X = importdata('BOLDnoise_1.txt', '\t', 1);
X = (X.data)';

[N,T] = size(X);

% estimate the mask
Mask = zeros(N,N);
for i=1:N
    [beta_al, beta_new_n, beta2_al, beta2_new_n] = betaAlasso_grad_2step(X([1:i-1 i+1:N],:), X(i,:), 0.65^2, log(T)/2); % 0.7^2
    Mask(i,1:i-1) = abs(beta2_al(1:i-1)) >0.01;
    Mask(i,i+1:N) = abs(beta2_al(i:N-1)) >0.01;
end
Mask = Mask + Mask';
Mask = (Mask~=0);

% perform constained_ICA
[y_m, W_m, WW_m, Score] = sparseica_W_adasize_Alasso_mask(log(T)/2, Mask , X);

figure, subplot(1,3,1), imagesc(WW_m); colorbar; subplot(1,3,2), imagesc(W_m); colorbar; subplot(1,3,3), imagesc(eye(N) - W_m); colorbar;