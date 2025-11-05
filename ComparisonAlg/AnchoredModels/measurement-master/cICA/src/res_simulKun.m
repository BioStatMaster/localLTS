clear;
resDir = '../result/simulKun/';
N = 8;
mseA = zeros(20,2);
T = [500 2000];
setting = 2;
switch setting
    case 1
        K = 4;
    case 2
        K = 5;
end
A_est = zeros(N,K,20,2);
X_all_500 = zeros(N,500,20);
X_all_2000 = zeros(N,2000,20);
X_til_all_500 = zeros(N,500,20);
X_til_all_2000 = zeros(N,2000,20);
for t_i = 1:2
    for expId = 1:20
        if setting == 1
            fileName = sprintf('%sN%d_K%d_T%d_expId%d_set%d.mat', resDir, N, K, T(t_i), expId, setting);
        else
            fileName = sprintf('%sN%d_K%d_T%d_expId%d.mat', resDir, N, K, T(t_i), expId);
        end
        load(fileName);
        if t_i==1
            X_all_500(:,:,expId) = X;
            X_til_all_500(:,:,expId) = X - Noise;
        elseif t_i==2
            X_all_2000(:,:,expId) = X;
            X_til_all_2000(:,:,expId) = X - Noise;
        end
        A_mn = A_m(:,1:K).*repmat(std(Err(1:K,:),0,2)',N,1);
        scale = sqrt(sum(wHat(1:K,:).*sigmaHat(1:K,:).^2,2));
        A_Hat(:,1:K) = A_Hat(:,1:K).*repmat(scale',[N 1]);
        A_est(:,:,expId,t_i) = A_Hat(:,1:K);
        mseA(expId,t_i) = mse(A_mn(:,1:K),A_Hat(:,1:K));
    end
end
save([resDir,'result_simul' num2str(setting) '.mat']);
boxplot(mseA([1:17,19:20],:));