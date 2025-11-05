% generate the data
T = 200;
entry_should0 = [];
entry_not0 = [];
for epoch = 1:1 % 40
    % WW
    WW = rand(5,5)*.45 + .05;
    WW = tril(WW) .* sign(rand(5,5)-.5);
    AA = inv(WW);
    % sources
    s = [];
    for i=1:5
        ss = normrnd(0,1,1,T);
        exponent = rand/2 + 1.5;
        var_ss = rand*.8 + .2;
        ss = sign(ss) .* abs(ss).^exponent;
        ss = ss - mean(ss);
        ss = ss/std(ss) * sqrt(var_ss);
        s = [s; ss];
    end
    % X
    X = AA * s;
    
    % with adaptive Lasso
    % [y, W, Score] = sparseica_W_adasize_Alasso(log(T)*.75, X);W, %pause,
	[y, W, Score] = sparseica_W_subgra_Alasso(log(T), X);
    
    % with OBD
%     [y, W, W_back, sal_back] = sparseica_path_Sol(X); 

    % compare the new and original W
    cor_sy = corr(s', y');
    %     WA = W*AA;
    new_W = [];
    for i=1:5
        [temp, I] = max(abs(cor_sy(i,:)));
        new_W = [new_W; W(I,:)];
    end

    for i=1:5
        for j=1:5
            if j>i
                entry_should0 = [entry_should0 new_W(i,j)];
            else
                entry_not0 = [entry_not0 new_W(i,j)];
            end
        end
    end
end

False_positive = sum(abs(entry_should0)>1E-3)/length(entry_should0);
False_negative = sum(abs(entry_not0)<1E-5)/length(entry_not0);

% for adaptive Lasso:
sum(abs(entry_should0)>1E-2)/length(entry_should0);
sum(abs(entry_not0)<1E-5)/length(entry_not0);
