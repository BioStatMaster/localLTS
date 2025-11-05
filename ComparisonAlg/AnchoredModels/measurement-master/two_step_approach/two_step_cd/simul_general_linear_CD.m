% SImulation for the two-step method for linear causal discovery with
% confounders and feedbacks

N = 10;
T = 500;
Num_conf = ceil(N/10); 
s = 2 * rand(N,T) - 1;
conf = 2* rand(10,T) - 1;

BB = .8*(2*rand(N,N) - 1);
BB_nonzero = (rand(N,N)>0.85); % there is an edge with prob. 0,15
for i=1:N 
    for j=i:N
        BB_nonzero(i,j) = 0;
    end
end
BB = BB .* BB_nonzero;

BB(5,7) = 0.3;
% BB(15,18) = -0.3;
BB(9,10) = 0.3;
% BB(16,19) = 0.4;

Ordering = randperm(N);
% B = BB(Ordering, :);
B = BB;

X = inv( eye(N)-B) * s;

Ind_conf = randperm(N);
Ind_conf = Ind_conf(1:Num_conf);

for i=1:floor(Num_conf/2)
    X(Ind_conf(2*i-1),:) = X(Ind_conf(2*i-1),:) + conf(i);
    X(Ind_conf(2*i),:) = X(Ind_conf(2*i),:) + 1.3* conf(i);
end

% 
% X(1,:) = s(1,:);
% for i = 2:N
%     X(i,:) = 0.7 * X(i-1,:) + s(i,:);
% end

% % confounder
% X(4,:) = X(4,:) + .6*conf;
% X(6,:) = X(6,:) + .5*conf;
% 
% % feedback
% X(8,:) = 0.7*X(8,:) + 0.5 * X(9,:);

% True Mask
Mask_t = zeros(N,N);
Mask_t = ((B+B') ~=0);
Mask_t(Ind_conf,Ind_conf) = 1;
Mask_t(1:N+1:end) = 0;
% for i=1:N
%     if i~=1
%         Mask(i,i-1) = 1;
%         Mask(i-1,i) = 1;
%     elseif i~=N
%         Mask(i,i+1) = 1;
%         Mask(i+1,i) = 1;
%     end
% end
% Mask(4,6) = 1;
% Mask(6,4) = 1;
% Mask(7,9) = 1;
% Mask(9,7) = 1;

% standardization
X = diag(1./std(X')) * X;

% estimate the mask
Mask = zeros(N,N);
for i=1:N
    [beta_al, beta_new_n, beta2_al, beta2_new_n] = betaAlasso_grad_2step(X([1:i-1 i+1:N],:), X(i,:), 0.5^2, log(T)/2); % 0.7^2
    Mask(i,1:i-1) = abs(beta2_al(1:i-1)) >0.01;
    Mask(i,i+1:N) = abs(beta2_al(i:N-1)) >0.01;
end
Mask = Mask + Mask';
Mask = (Mask~=0);