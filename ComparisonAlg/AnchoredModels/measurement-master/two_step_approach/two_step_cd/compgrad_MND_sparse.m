% Compute gradient relative to weights

for subnet=1:ninputs
    for column=1:ninputs
        grad12{subnet}(:,column) = sum(jback2i{subnet,column},2);
        grad23{subnet}(:,1:nhidden) = grad23{subnet}(:,1:nhidden) + jback3{subnet,column} * jacob2o{subnet,column}';
        grad13{subnet}(:,column) = sum(jback3{subnet,column},2);
        grad34{subnet}(:,1) = grad34{subnet}(:,1) + jback4i{subnet,column} * jacob3{subnet,column}';
        grad45{subnet}(:,1:nextra) = grad45{subnet}(:,1:nextra) + jback5{column}(subnet,:) * jacob4o{subnet,column}';
    end
    grad12{subnet} = grad12{subnet} + back2{subnet} * output1';
    grad12{subnet}(:,1:ninputs) = grad12{subnet}(:,1:ninputs) + wdecayf12 * weight12{subnet}(:,1:ninputs);
    grad23{subnet} = grad23{subnet} + back3{subnet} * output2{subnet}';
    grad23{subnet}(:,1:nhidden) = grad23{subnet}(:,1:nhidden) + wdecayf23 * weight23{subnet}(:,1:nhidden);
    grad13{subnet} = grad13{subnet} + back3{subnet} * output1';
    grad34{subnet} = grad34{subnet} + back4{subnet} * output3{subnet}';
    grad34{subnet}(:,1) = grad34{subnet}(:,1) + wdecayf34 * weight34{subnet}(:,1);
end   

% with the MND regularization! The regularization coefficient is lambda

% with the SCAD penalty for each direct connection from the input to the
% output. Parameters are lambda_in and a.
% dev is the gradient of the objective function w.r.t. the weights
dev = zeros(ninputs,ninputs); m = 180;

for subnet=1:ninputs
    K{subnet} = cal_K(output3{subnet}(1,:),output1(1:ninputs,:));
    grad23{subnet} = K{subnet} * output2{subnet}' * lambda + grad23{subnet};
    
    % incorporating the SCAD penalty
    for col = 1:ninputs
        if abs(weight13{subnet}(col)) <= lambda_in
            dev(subnet, col) = tanh(m*weight13{subnet}(col));
        else if abs(weight13{subnet}(col)) <= a * lambda_in
                dev(subnet, col) = sign(weight13{subnet}(col))*tanh(m*lambda_in)*(a*lambda_in - abs(weight13{subnet}(col)))/(a-1)/lambda_in;
            end
        end
    end
    grad13{subnet} = K{subnet} * output1' * lambda +grad13{subnet} + lambda_in * [dev(subnet,:) 0] *T;
    for i = 1:ntrain
        grad12{subnet} = diag(g_dev1_2pi_atan(input2{subnet}(:,i))) * weight23{subnet}(:,1:nhidden)'*K{subnet}(i)*output1(:,i)' * lambda + grad12{subnet};
    end    
end