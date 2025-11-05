addpath('code');
addpath('FastICA_23/');
resDir0 = '../result/simulKun_Lingam/';
dataDir = '../result/simulKun/';
setting = 1;
if ~exist(resDir0, 'dir')
    mkdir(resDir0);
end
load([dataDir,'result_simul' num2str(setting) '.mat']);
rng('default');
if setting == 1 % r = 3, l = 5; chain;
    rng(1);
    aa = sign(rand(1,11)-0.3).*(.7*rand(1,11) + .3);
    B = [0 0 0 0 0 0 0 0; aa(1) 0 0 0 0 0 0 0; 0 aa(2) 0 0 0 0 0 0; 0 0 aa(3) 0 0 0 0 0; aa(4) aa(5) 0 0 0 0 0 0; 0 aa(6) aa(7) 0 0 0 0 0; 0 0 aa(8) aa(9) 0 0 0 0;...
        aa(10) 0 0 aa(11) 0 0 0 0];
elseif setting ==2 %
    B = [0 0 0 0 0 0 0 0; 0.6 0 0 0 0 0 0 0; 0 .7 0 0 0 0 0 0; .7 .4 0 0 0 0 0 0; .5 .4 0 0 0 0 0 0; 0 .4 .7 0 0 0 0 0; 0.6 0 0.5 -.4 0 0 0 0;...
        0.4 -0.4 0.5 0 .6 0 0 0];
end
B_est_til = zeros(N,N,20,2);        
B_est = zeros(N,N,20,2);
for t_i = 1:2
    for expId = 1:20
        if t_i == 1
            [B0 stde ci k] = estimate(X_all_500(:,:,expId));
            B_est(:,:,expId,t_i) = prune(X_all_500(:,:,expId),k);
            [B0 stde ci k] = estimate(X_til_all_500(:,:,expId));
            B_est_til(:,:,expId,t_i) = prune(X_til_all_500(:,:,expId),k);
        elseif t_i == 2
            [B0 stde ci k] = estimate(X_all_2000(:,:,expId));
            B_est(:,:,expId,t_i) = prune(X_all_2000(:,:,expId),k);
            [B0 stde ci k] = estimate(X_til_all_2000(:,:,expId));
            B_est_til(:,:,expId,t_i) = prune(X_til_all_2000(:,:,expId),k);
        end
    end
end
save([resDir0,'result_simul' num2str(setting) '.mat']);

