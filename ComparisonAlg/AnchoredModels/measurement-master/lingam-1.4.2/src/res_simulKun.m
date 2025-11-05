clear;
resDir0 = '../result/simulKun_Lingam/';
setting = 1;
load([resDir0,'result_simul' num2str(setting) '.mat']);
Precision = zeros(2,2);
Recall = zeros(2,2);
for t_i = 1:2
    B_est_i = B_est(:,:,:,t_i);
    TP = sum(sum(sum((B_est_i~=0) .* repmat(B~=0,[1 1 size(B_est_i,3)]))));
    Precision(t_i,1) = TP/sum(sum(sum((B_est_i~=0))));
    Recall(t_i,1) = TP/sum(sum(sum(repmat(B~=0,[1 1 size(B_est_i,3)]))));
    B_est_i = B_est_til(:,:,:,t_i);
    TP = sum(sum(sum((B_est_i~=0) .* repmat(B~=0,[1 1 size(B_est_i,3)]))));
    Precision(t_i,2) = TP/sum(sum(sum((B_est_i~=0))));
    Recall(t_i,2) = TP/sum(sum(sum(repmat(B~=0,[1 1 size(B_est_i,3)]))));
end
fprintf('directed graph\n')
Precision
Recall

Precision = zeros(2,2);
Recall = zeros(2,2);
B = B + B';
for t_i = 1:2
    B_est_i = B_est(:,:,:,t_i);
    for i = 1:size(B_est_i,3)
        B_est_i(:,:,i) = B_est_i(:,:,i) + B_est_i(:,:,i)';
    end
    TP = sum(sum(sum((B_est_i~=0) .* repmat((B~=0),[1 1 size(B_est_i,3)]))));
    Precision(t_i,1) = TP/sum(sum(sum((B_est_i~=0))));
    Recall(t_i,1) = TP/sum(sum(sum(repmat(B~=0,[1 1 size(B_est_i,3)]))));
    B_est_i = B_est_til(:,:,:,t_i);
    for i = 1:size(B_est_i,3)
        B_est_i(:,:,i) = B_est_i(:,:,i) + B_est_i(:,:,i)';
    end
    TP = sum(sum(sum((B_est_i~=0) .* repmat(B~=0,[1 1 size(B_est_i,3)]))));
    Precision(t_i,2) = TP/sum(sum(sum((B_est_i~=0))));
    Recall(t_i,2) = TP/sum(sum(sum(repmat(B~=0,[1 1 size(B_est_i,3)]))));
end
fprintf('undirected graph\n')
Precision
Recall