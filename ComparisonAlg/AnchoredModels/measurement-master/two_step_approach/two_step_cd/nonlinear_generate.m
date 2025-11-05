function [x_nonlinear, MSE, A_fit] = nonlinear_generate(ninputs, e)
% generate NONLINEAR mixtures from e

% create and initialize network variables
weightrange = 3;
If_direct = 0;
nhidden = 10;

for subnet=1:ninputs
    weight12{subnet} = (rand(nhidden,ninputs+1) - .5) * (2 * weightrange);
    weight23{subnet} = (rand(1,nhidden+1) - .5) * (2 * weightrange);
    if(If_direct)
        weight13{subnet} = (rand(1,ninputs+1) - .5) * (2 * weightrange);
    else
        weight13{subnet} = (rand(1,ninputs+1) - .5) * 0;
    end
end

% feedforward to calculate the nonlinear outputs
output1(1:ninputs,:) = e;
output1(ninputs+1,:) = 1;

for subnet=1:ninputs
   input2{subnet} = weight12{subnet} * output1;
   %output2{subnet}(1:nhidden,:) = tanh(input2{subnet});           % for tanh sigmoids
   output2{subnet}(1:nhidden,:) = (2/pi) * atan(input2{subnet});   % for arctangent sigmoids
   output2{subnet}(nhidden+1,:) = 1;

   input3{subnet} = weight23{subnet} * output2{subnet} + weight13{subnet} * output1;
   x_nonlinear(subnet,:) = input3{subnet};
end

ee = [e; ones(1,length(e))];
A_fit = (x_nonlinear * ee') * inv(ee*ee');
MSE = var((x_nonlinear-A_fit*ee)')./var(x_nonlinear');