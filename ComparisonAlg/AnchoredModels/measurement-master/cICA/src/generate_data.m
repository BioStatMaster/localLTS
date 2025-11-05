function generate_data(pars)
% generate data according to random graph structures

% load graph_structures % : G_record
load(pars.graphFile);

% all variables have dimension 1
rng(1)
N = pars.N;
T = pars.T;
count = 0;
func = [1]; % 1: linear, 2: sinc, 3: cos, 4: tanh
noise = [4]; % 1: Gaussian, 2: uniform, 3: Gamma, 4: MoG
for i = 1:size(G_record,2)
    G = G_record{i};
    B = G;
    for trial = 1:1
        X=zeros(T,N);
        for j = 1:N
            if(j==1)
                noise_id = noise(randi(length(noise),1));
                if(noise_id==1)
                    X(:,j) = randn(T,1);
                end
                if(noise_id==2)
                    X(:,j) = (rand(T,1)-0.5);
                end
                if noise_id==4
                    mu = [0 0];
                    w = [0.8 0.2];
                    sigma = [1e-2 1];
                    MU = mu';
                    SIGMA = reshape(sigma,1,1,2);
                    P = w;
                    obj = gmdistribution(MU,SIGMA,P);
                    E = random(obj, T);
                    X(:,j) = E;
                end
            else
                PA = find(G(:,j)==1);
                nPA = length(PA); % number of parents
                func_id = func(randi(length(func),1));
                noise_id = noise(randi(length(noise),1));
                if(func_id==1)
                    B(PA,j) = ones(nPA,1)*1.7/(nPA+1);
                    X(:,j) = X(:,PA)*B(PA,j);
                end
                if(func_id==2)
                    X(:,j) = sin(X(:,PA)*ones(nPA,1));
                end
                if(func_id==3)
                    X(:,j) = cos(X(:,PA)*ones(nPA,1));
                end
                if(func_id==4)
                    X(:,j) = tanh(X(:,PA)*ones(nPA,1));
                end
                if(noise_id==1)
                    X(:,j) = X(:,j) + 0.4*randn(T,1);
                end
                if(noise_id==2)
                    X(:,j) = X(:,j) + 0.5*(rand(T,1)-0.5);
                end
                if noise_id==4
                    mu = [0 0];
                    w = [0.8 0.2];
                    sigma = [1e-2 1];
                    MU = mu';
                    SIGMA = reshape(sigma,1,1,2);
                    P = w;
                    obj = gmdistribution(MU,SIGMA,P);
                    E = random(obj, T);
                    X(:,j) = X(:,j) + E;
                end
            end
        end
        C = corr(X); % check the faithfulness by estimating the correlation matrix
        sign = 0;
        for ii = 1:N
            for jj = 1:N
                if(G(ii,jj)==1)
                    if(abs(C(ii,jj))<0.3 | abs(C(ii,jj))>0.9)
                        sign = 1;
                        break;
                    end
                end
            end
            if(sign)
                break;
            end
        end
        if(~sign)
            fprintf('%d %d\n',i,trial);  
            count = count+1;
            G_save{count} = G;
            B_save{count} = B;
            Data_save{count} = X;
        end
    end
end
save(pars.dataFile, 'G_save', 'B_save', 'Data_save');
