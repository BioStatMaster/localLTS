%% evaluation (given B as the ground-truth and W_m)

W_m = W_m - diag(diag(W_m));
% orientations the method fails to discover
pp1 = sum(sum( (abs(B)>0.05) .* ( abs(W_m)<0.002)))/sum(sum( (abs(B)>0.05))),

% extra undirected edges
%% pp2 = sum(sum( ((abs(B)+abs(B')<0.01)) .* ( (abs(W_m)+abs(W_m'))>0.02 )))/sum(sum( ((abs(B)+abs(B')<0.01))))
pp2 = sum(sum( (abs(B + B')<0.03) .* ( abs(W_m)>0.05) ))/sum(sum( (abs(W_m)>0.05))),

% extra directed edges
%% pp3 = sum(sum( (abs(B)<0.05) .* ( abs(W_m)>0.02)))/sum(sum( (abs(B)<0.05)))
pp3 = sum(sum( (abs(B)<0.03) .* ( abs(W_m)>0.05) ))/sum(sum( (abs(W_m)>0.05))),

% missing undirected edges:
pp4 = sum(sum( (abs(B)>0.05) .* ( (abs(W_m)+abs(W_m'))<0.001)))/sum(sum( (abs(B)>0.05)))

% missing directed edges:
pp5 = sum(sum( (abs(B)>0.05) .* ( abs(W_m)<0.001)))/sum(sum( (abs(B)>0.05)))