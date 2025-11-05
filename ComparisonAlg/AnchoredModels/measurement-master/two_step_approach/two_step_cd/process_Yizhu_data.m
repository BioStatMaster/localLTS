% process Yizhu's data

data = load('scer_spardata_yizhu_12042016.txt');

data(:,1) = log(data(:,1) + 80);
data(:,2) = log(data(:,2) + 50);
data(:,5) = log(data(:,5) + 1);
data(:,7) = log(data(:,7) + 0.5);
data(:,8) = log(data(:,8) + 15);
data(:,9) = -log(-data(:,9) + 3);
data(:,10) = -log(-data(:,10) + 3);
data(:,10) = -log(-data(:,10) + 5);
data(:,11) = -log(-data(:,11) + 8);
data(:,12) = -log(-data(:,12) + 8);

figure, for i=1:12 subplot(2,6,i), hist(data(:,i)); end

figure, 
for i=1:12
    for j=1:12  
        subplot(12,12,(i-1)*12+j), 
        plot(data(:,i), data(:,j), '.'); 
    end;
end

X = data';
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
[y_m, W_m, WW_m, Score] = sparseica_W_adasize_Alasso_mask(log(T)/2, Mask , X);

figure, subplot(1,2,1), imagesc(WW_m); colorbar; subplot(1,2,2), imagesc(W_m); colorbar;
