function [IF_causality, Bestcausal,causalperm] = check_causality( Best )


N = length(Best);
n = N;
% Best = (abs(Best)>Tol) .* Best;
causalperm = [];
Bestcausal = Best;
best = Best;
IF_causality = 1;

for i=1:N
    row = 1;
    while sum(abs(best(row,:)))
        row = row + 1;
        if row > N;
            IF_causality = 0;
            break;
        end
    end
    
    if ~IF_causality
        break;
    end    
    
    causalperm = [causalperm row];
    best(:,row) = 0; best(row,row) = 1;
    
    % exchange
%     temp1 = Bestcausal(i,:);
    Bestcausal(i,:) = Best(row,:);
%     Bestcausal(row,:) = temp1;
%     temp2 = Bestcausal(:,i);
%     Bestcausal(:,i) = Bestcausal(:,row);
%     Bestcausal(:,row) = temp2;
end


