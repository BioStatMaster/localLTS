function generate_strcucture(pars)
% generate random graphs with 10 varaibles
% for each variable, have 2 parents at most, and the total degree is at
% most 3

N = pars.N; % number of variables
nPA = pars.nPA; % restrict the maximum number of parents for each variable (we control the degree of parents)
min_PA = pars.min_PA; % restrict the mininum number of parents, except for nodes in the first layer
degree = pars.degree; % restrcit the total number of degree
ordering = 1:N; % fix the ordering: 1,2,3,4,5,6

count = 0;
for i=1:100 % generate 100 instance
    PA = cell(1,N);
    sign = 0;
    for j=N:-1:2 % start from the leaf node
        canPA = ordering(1:j-1); % candidate parents
        combs2 = [];
        for k=min(length(canPA),min_PA):min(length(canPA),nPA) % restrict the number of parents up to 4
            combs{k} = nchoosek(canPA,k); % for each varible, list all the possible parents
            if(size(combs{k},2)<nPA)
                combs{k}=[combs{k},zeros(size(combs{k},1),nPA-size(combs{k},2))];
            end
            combs2 = [combs2;combs{k}];
        end
        %         combs2 = [zeros(1,nPA);combs2]; % it is possible there is no parents
        PA{ordering(j)} = combs2(randi(size(combs2,1)),:); % the parents for each variable
    end
    % generate the graph structure
    G = zeros(N,N);
    for j=1:N
        idPA = find(PA{j});
        if(~isempty(idPA))
            G(PA{j}(idPA),j)=1;
        end
    end
    for j=1:N
        if(length(find(G(j,:)==1)) + length(find(G(:,j)==1)) > degree)
            sign=1;
            break;
        end
    end
    if(~sign)
        % total degree
        tdegree = [];
        for k = 1:N
            tdegree(k) = length(find(G(:,k)==1)) + length(find(G(k,:)==1));
        end
        mean(tdegree)
        if(mean(tdegree)>=2 & mean(tdegree)<3)
            count = count+1;
            G_record{count}=G;
        end
        
        %         % in degree
        %         idegree = [];
        %         for k = 1:N
        %             idegree(k) = length(find(G(:,k)==1));
        %         end
        %
        %         % out degree
        %         odegree = [];
        %         for k = 1:N
        %             odegree(k) = length(find(G(k,:)==1));
        %         end
    end
    clear PA canPA combs combs2 G
    
end

save(pars.graphFile, 'G_record');


