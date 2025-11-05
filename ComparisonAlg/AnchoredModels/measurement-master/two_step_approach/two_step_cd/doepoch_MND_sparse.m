% Do an epoch

testcost        % Test whether cost function has improved; if not, backtrack

wadapt          % Adapt weights

arrinit         % Initialize arrays

forward_zeromean         % Propagate forward

% since the SCAD penalty is hard to evaluate, it is neglected in evaluating the cost function...
costderiv_MND       % Compute cost function and its derivatives, for backprop

back            % Propagate backwards

% disp('MODIFIED BY ME, with regularization!!!'),
% compgrad_sparse        % Compute gradient relative to weights
compgrad_MND_sparse