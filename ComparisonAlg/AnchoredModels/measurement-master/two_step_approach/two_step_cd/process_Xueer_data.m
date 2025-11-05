%% load 'data_initial_process.Xueer.mat';
TT = 3683;

Cor_I_c = corr(I_data_txt, C_data_txt);
sum_total_Cor = sum(Cor_I_c.^2);
[Cor_sort, Ind_sort] = sort(sum_total_Cor, 'descend');

Num_sel = 500;
C_selected = C_data_txt(:,Ind_sort(1:Num_sel));
I_data_pro = I_data_txt - repmat(mean(I_data_txt), TT, 1);
C_data_pro = C_selected - repmat(mean(C_selected), TT, 1);
figure, imagesc(corr(I_data_pro, C_data_pro)); colorbar

%
I_data_pro = I_data_pro * diag(1./std(I_data_pro));
C_data_pro = C_data_pro * diag(1./std(C_data_pro));

X = [I_data_pro C_data_pro]';
[N,T] = size(X);

% estimate the mask
Mask = zeros(N,N);
for i=1:N
    [beta_al, beta_new_n, beta2_al, beta2_new_n] = betaAlasso_grad_2step(X([1:i-1 i+1:N],:), X(i,:), 0.75^2, log(T)/2); % 0.7^2
    Mask(i,1:i-1) = abs(beta2_al(1:i-1)) >0.01;
    Mask(i,i+1:N) = abs(beta2_al(i:N-1)) >0.01;
end
Mask = Mask + Mask';
Mask = (Mask~=0);

% perform constained_ICA
[y_m, W_m, WW_m, Score] = sparseica_W_adasize_Alasso_mask(log(3704)/2, Mask , X);

figure, subplot(1,2,1), imagesc(WW_m); colorbar; subplot(1,2,2), imagesc(W_m); colorbar;