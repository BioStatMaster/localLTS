% script to generate data with measurement noise and do causal analysis

%% generate noisy data
% first generate
addpath('../ifa');
resDir = '../result/simulKun/';
if ~exist(resDir, 'dir')
    mkdir(resDir);
end
Small = 1;
Real = 0;
Setting = 1;
if Small
    N = 8;
    for T = [500 2000]
        if Setting == 1 % r = 3, l = 5; chain;
            %         % aa = .7*rand(1,6) + .2;
            %          aa(1) = -.2;
            %         B = [0 0 0 0 0 0 0 0; aa(1) 0 0 0 0 0 0 0; 0 aa(2) 0 0 0 0 0 0; .6 aa(3) 0 0 0 0 0 0; .4 0 .6 0 0 0 0 0; 0.6 aa(4) -.75 0 0 0 0 0; 0.6 0 aa(5) 0 0 0 0 0;...
            %         0 aa(6) 0.6 0 0 0 0 0];
%             B = [0 0 0 0 0 0 0 0; -0.2 0 0 0 0 0 0 0; 0 .5454 0 0 0 0 0 0; .6 .4374 0 0 0 0 0 0; .4 0 .6 0 0 0 0 0; 0.6 .7015 -.75 0 0 0 0 0; 0.6 0 .6171 0 0 0 0 0;...
%                 0 .5160 0.6 0 0 0 0 0];
            rng(1);
            aa = sign(rand(1,11)-0.3).*(.7*rand(1,11) + .3);
            B = [0 0 0 0 0 0 0 0; aa(1) 0 0 0 0 0 0 0; 0 aa(2) 0 0 0 0 0 0; 0 0 aa(3) 0 0 0 0 0; aa(4) aa(5) 0 0 0 0 0 0; 0 aa(6) aa(7) 0 0 0 0 0; 0 0 aa(8) aa(9) 0 0 0 0;...
                    aa(10) 0 0 aa(11) 0 0 0 0];
            K = 4;
        elseif Setting ==2 %
            B = [0 0 0 0 0 0 0 0; 0.6 0 0 0 0 0 0 0; 0 .7 0 0 0 0 0 0; .7 .4 0 0 0 0 0 0; .5 .4 0 0 0 0 0 0; 0 .4 .7 0 0 0 0 0; 0.6 0 0.5 -.4 0 0 0 0;...
                0.4 -0.4 0.5 0 .6 0 0 0];
            K = 5;
        end
        
        for expId = 1:20
            rng(expId);
            fileName = sprintf('%sN%d_K%d_T%d_expId%d_set%d.txt', resDir, N, K, T, expId, Setting);
            if exist(fileName, 'file')
                continue;
            else
                fid = fopen(fileName,'w');
                fprintf(fid, 'Created on %s', char(datetime));
                fclose(fid);
                Err = normrnd(0,1,N,T);
                
                % non-Gaussian
                for i=1:N
                    Expon = rand * .6 + 1.4;
                    Err(i,:) =  (abs(Err(i,:))).^Expon .* (sign(Err(i,:)));
                end
                A_m =  inv(eye(N) - B);
                X_tilde = A_m * Err;
                
                % noise_std = 0.6 * ones(N,1);
                noise_std = ones(N,1);
                Noise = diag(noise_std) * normrnd(0,1,N,T);
                
                X = X_tilde + Noise;
                % FA for pre-processing
                %     K = 2;
                % [A1,Ph1,LL]=ffa(X',K,100,1E-5);
                %if Setting == 1
                    %[A1,Y1,sigma1,LL] = fa_em(X,K);
                    %sigma1,
                    %[L1, Ph1] = factoran(X', K);
                    %L1 = diag(std(X')) * L1;
                    %Ph1 = (var(X'))' .* Ph1;
                    %Ph1,
                %end
                r = K;
                D = N;
                n = D + r;
                
                A = [A_m(:,1:r) eye(D)];
                mu = repmat([0 0],[n,1]);
                w = repmat([0.8 0.2],[n 1]);
                sigma = repmat([1e-1 1.5],[n 1]);
                w(n-D+1:n,:) = repmat([1 0],[D 1]);
                sigma(n-D+1:n,:) = 2;
                
                % experiments
                parsEM.dimG = D;
                parsEM.nEye = D;
                parsEM.thres = 1e-6;
                parsEM.minIter = 500;
                parsEM.maxIter = 5000;
                initNoise = 0.8;
                parsEM.A = A+initNoise*(rand(D,n)-1/2);
                parsEM.A(:,end-parsEM.nEye+1:end) = eye(parsEM.nEye);
                parsEM.mu = mu;
                parsEM.sigma = sigma;
                parsEM.w = w;
                parsEM.updatePrior = 1;
                parsEM.zeroMean = 0;
                parsEM.noise = 1e-2;
                parsEM.fast = 0;
                tic;
                [A_Hat, sigmaHat, wHat, loglAll, posMeanE] = ifaEM_MN(X, 2, parsEM);
                toc;
                fileName = sprintf('%sN%d_K%d_T%d_expId%d_set%d.mat', resDir, N, K, T, expId, Setting);
                save(fileName, 'Err', 'Noise', 'A_m', 'A_Hat', 'sigmaHat', 'wHat', 'loglAll', 'posMeanE', 'X');
            end
        end
    end
    % Initialization really matters!  Will focus on it for the fMRI data.
    
    %% simul 1: non-Gaussian tilde{E}, gaussian E, two independent causes: B = [0 0 0 0 0 0 0; 0 0 0 0 0 0 0; .5 .6 0 0 0 0 0; .7 .8 0 0 0 0 0 ; .5 .4 0 0 0 0; .8 .4 0 0 0 0 0; 0 .7 0 0 0 0 0];
    %% simul 2: non-Gaussian tilde{E}, gaussian E, two dependent causes: B = [0
    %% 0 0 0 0 0 0; 0.5 0 0 0 0 0 0; .5 .6 0 0 0 0 0; .7 .8 0 0 0 0 0 ; .5 .4 0 0 0 0 0; .8 .4 0 0 0 0 0; 0 .7 0 0 0 0 0];
    %% simul 3: non-Gaussian tilde{E}, gaussian E, three dependent causes: B =
    %% [0 0 0 0 0 0 0; 0.5 0 0 0 0 0 0; .5 .6 0 0 0 0 0; .7 .8 0 0 0 0 0 ; .5 .4 .5 0 0 0 0; .8 .4 0 0 0 0 0; 0 .7 0 0 0 0 0];
else
    % large scale simulations...
    
    N = 20;
    T = 5000;
    B1 = 1.2* (rand(N,N) - .5); for i=1:N B1(i,i:N) = 0; end
    l = 8;
    B1(N-(l-1):N, N-(l-1):N) = 0;
    
    Err = normrnd(0,1,N,T);
    
%     % non-Gaussian
%     for i=1:N
%         Expon = rand * .6 + 1.4;
%         Err(i,:) =  (abs(Err(i,:))).^Expon .* (sign(Err(i,:)));
%     end
    A_m =  inv(eye(N) - B1);
    X_tilde = A_m * Err;
    
    % noise_std = 0.6 * ones(N,1);
    noise_std = .5*ones(N,1);
    Noise = diag(noise_std) * normrnd(0,1,N,T);
    
    X = X_tilde + Noise;
    % FA for pre-processing
    K = N - l;
    % [A1,Ph1,LL]=ffa(X',K,100,1E-5);
    
    % [A1,Y1,sigma1,LL] = fa_em(X,K);
    % figure, plot(sigma1), title('\sigma')
%     [L,Ph,LL]=ffa(X',K,300,0.0000001); figure, plot(sqrt(Ph)), title('var')
    [L, Ph] = factoran(X', K);  % this is only for the correlation matrix... we want to use covrariance matrix
    L = diag(std(X')) * L;
    Ph = (var(X'))' .* Ph;
    figure, plot(sqrt(Ph), '-+'), title('std(E*)')
end


if Real
%% Real data: autism and typical data
Data_index = [272:276 278 282 284 285 291:295 297 298 300:302 304 310 312 314 315  318:321 324 325 327 329 330:333 336:345 347:370 372:377 379 380 381];
% In total there are 79 (30 autistic and 49 typical) subject, for each of
% which there are 295 data points.
Data= [];
for i=1:79
    if i~=47
    if i<= 30
    fname = sprintf('UM_1_autism_0050%3d_rois_aal.1D', Data_index(i));
    else
        fname = sprintf('UM_1_typical_0050%3d_rois_aal.1D', Data_index(i));
    end
    data1 = importdata(fname);
    TT = length(data1.data);
    Data = [Data; data1.data - repmat(mean(data1.data), TT, 1)];
    end
end

[T_total,N] = size(Data);
T_each = T_total/(79-1);
% the two groups
Data1 = Data(1:30*T_each,:);
Data2 = Data(1+30*T_each:end,:);

% FA
% [L1,Ph1,LL1]=ffa(Data1, 116-18,10000,1E-10); figure, plot(sqrt(Ph1)); 
% [L2,Ph2,LL2]=ffa(Data2, 116-18,10000,1E-10); figure, plot(sqrt(Ph2)); 
% [Lo,Pho,LLo]=ffa(Data, 116-18,10000,1E-10); figure, plot(sqrt(Pho)); 

[L1, Ph1] = factoran(Data1, 116-30); L1 = diag(std(Data1)) * L1; Ph1 = (var(Data1))' .* Ph1; figure, plot(sqrt(Ph1), '-+'), title('std(E*)')
[L2, Ph2] = factoran(Data2, 116-30); L2 = diag(std(Data2)) * L2; Ph2 = (var(Data2))' .* Ph2; figure, plot(sqrt(Ph2), '-+'), title('std(E*)')
[Lo, Pho] = factoran(Data, 116-30); Lo = diag(std(Data)) * Lo; Ph2o= (var(Data))' .* Pho; figure, plot(sqrt(Pho), '-+'), title('std(E*)')


Cov1 = Data1'*Data1/length(Data1);
Cov2 = Data2'*Data2/length(Data2);

Cor1 = diag(1./sqrt(diag(Cov1))) * Cov1 * diag(1./sqrt(diag(Cov1)));
Cor2 = diag(1./sqrt(diag(Cov2))) * Cov2 * diag(1./sqrt(diag(Cov2)));

Cov1_s = Cov1 - diag(Ph1);
Cov2_s = Cov2 - diag(Ph2);

Cor1_s = diag(1./sqrt(diag(Cov1_s))) * Cov1_s * diag(1./sqrt(diag(Cov1_s)));
Cor2_s = diag(1./sqrt(diag(Cov2_s))) * Cov2_s * diag(1./sqrt(diag(Cov2_s)));

% or
Cov1_s = Cov1 - diag(Pho);
Cov2_s = Cov2 - diag(Pho);

% save files
i = 1; mm = Cov1_s(i,1:i); dlmwrite('FileCov1_s.txt',  mm,'delimiter','\t','precision','%.4f'); 
for i=2:116 mm = Cov1_s(i,1:i); dlmwrite('FileCov1_s.txt',  mm,'delimiter','\t','precision','%.4f', '-append'); end

i = 1; mm = Cov2_s(i,1:i); dlmwrite('FileCov2_s.txt',  mm,'delimiter','\t','precision','%.4f'); 
for i=2:116 mm = Cov2_s(i,1:i); dlmwrite('FileCov2_s.txt',  mm,'delimiter','\t','precision','%.4f', '-append'); end

Cov_t = Data' * Data/length(Data);
i = 1; mm = Cov_t(i,1:i); dlmwrite('FileCov.txt',  mm,'delimiter','\t','precision','%.4f'); 
for i=2:116 mm = Cov_t(i,1:i); dlmwrite('FileCov.txt',  mm,'delimiter','\t','precision','%.4f', '-append'); end

i = 1; mm = Cov1(i,1:i); dlmwrite('FileCov1.txt',  mm,'delimiter','\t','precision','%.4f'); 
for i=2:116 mm = Cov1(i,1:i); dlmwrite('FileCov1.txt',  mm,'delimiter','\t','precision','%.4f', '-append'); end

i = 1; mm = Cov2(i,1:i); dlmwrite('FileCov2.txt',  mm,'delimiter','\t','precision','%.4f'); 
for i=2:116 mm = Cov2(i,1:i); dlmwrite('FileCov2.txt',  mm,'delimiter','\t','precision','%.4f', '-append'); end


end

% % how can we find non-singular submatrix?
% [Ph1_sort, Ind1] = sort(Ph1, 'ascend');
