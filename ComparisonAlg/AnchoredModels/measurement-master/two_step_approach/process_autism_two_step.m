% process the autism data

Data_index = [272:276 278 282 284 285 291:295 297 298 300:302 304 310 312 314 315  318:321 324 325 327 329 330:333 336:345 347:370 372:377 379 380 381];
% In total there are 79 (30 autistic and 49 typical) subject, for each of
% which there are 295 data points.
Data= [];
for i=1:79
    if i<= 30
    fname = sprintf('UM_1_autism_0050%3d_rois_aal.1D', Data_index(i));
    else
        fname = sprintf('UM_1_typical_0050%3d_rois_aal.1D', Data_index(i));
    end
    data1 = importdata(fname);
    TT = length(data1.data);
    Data = [Data; data1.data - repmat(mean(data1.data), TT, 1)];
end

[T_total,N] = size(Data);

%Data = Data - repmat(mean(Data), T_total, 1);
Data = Data * diag(1./std(Data));
[B,W_m,y_m] = two_step_CD(Data');

% now we have the values for initialization in init_autism.mat

for i = 1:79
    i,
    load init_autism.mat; Mask = Mask - diag(diag(Mask));
    X = Data((i-1)*TT+1:i*TT,:)';
%     [y_f{i}, W_f{i}, WW_f{i}, Score] = sparseica_W_adasize_Alasso_mask_init(log(TT)/2/2, Mask, WW_m_sparse4, X);
% B{i} = eye(N) - W_f{i};
%     B{i} = B{i} .* (abs(B{i}) > 0.04);
    W_noSparse{i} = natural_grad_Adasize_Mask(X,Mask,WW_m_sparse4);
    B_noSparse{i} = eye(N) - W_noSparse{i};
    B_noSparse{i} = B_noSparse{i} .* (abs(B_noSparse{i}) > 0.04);
end

% refine the initialization value
W_init = natural_grad_Adasize_Mask(Data',Mask,WW_m_sparse4);
save W_init_autism.mat W_init

% group analysis
for i = 1:79
    % subject 47 is ill-posed; removed for now
    if i<47
    x_all{i} = Data((i-1)*TT+1:i*TT,:)';
    elseif i>47
        x_all{i-1} = Data((i-1)*TT+1:i*TT,:)';
    end
end

Mask = Mask - diag(diag(Msk));
[y_group, W_group] = W_adasize_groupwise(1 * log(TT)/2 , Mask, W_init, x_all);

% Data_index(47) = []; for i=1:78 B_less_similar{i} = eye(116) -
% W_group{i};  fname1 = sprintf('B_0050%3d_rois_aal.txt', Data_index(i));
% Bt = B_less_similar{i}; save(fname1, 'Bt', '-ASCII'); end

% %%% for similation of the group analysis
% ss = 2*(rand(900,10)) - 1;
% W1 = rand(10,10) - .5; W1 = W1 -diag(diag(W1)) + eye(10);
% W2 = W1; W3 = W1;
% W2(2,3) = W1(2,3) + .3;
% W3(5,4) = W1(5,4) - .3;
% xx(1:300,:) = ss(1:300,:) * inv(W1);
% xx(301:600,:) = ss(301:600,:) * inv(W2);
% xx(601:900,:) = ss(601:900,:) * inv(W3);
% x_all{1} = xx(1:300,:)'; x_all{2} = xx(301:600,:)'; x_all{3} = xx(601:900,:)';
% W_init = natural_grad_Adasize_Mask(xx', abs(ones(10,10) - eye(10))>1E-5);
% [y_group, W_group] = W_adasize_groupwise(1 * log(300)/2 , abs(ones(10,10) - eye(10))>1E-5, W_init, x_all);