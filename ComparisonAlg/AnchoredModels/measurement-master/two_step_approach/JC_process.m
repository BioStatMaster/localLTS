% Yizhu's data
clear all,clc,close all
addpath(genpath(pwd))
Data = csvread('JC_processed_WS6_header_NoCodon.csv',2,2);
Data = Data(:,1:17);
[T,N] = size(Data);

Data = Data - repmat(mean(Data), T, 1);
Data = Data * diag(1./std(Data));
[B,W,WW,y,Mask] = two_step_CD(Data');

Mask = Mask - diag(diag(Mask));
% % refine the initialization value
W_init = natural_grad_Adasize_Mask(Data',Mask,WW);

lambda1 = 0.2;
lambda2 = 1;

x_all{1} = Data';
[y_final, W_final] = W_adasize_groupwise_sparse(lambda1 * log(T)/2 ,lambda2 * log(T)/2 , Mask, W_init, x_all);
BB = eye(N) - W_final{1};
save JC BB Mask

