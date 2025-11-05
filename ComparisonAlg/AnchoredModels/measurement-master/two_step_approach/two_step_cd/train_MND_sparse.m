% Train a net

figure(1), clf
a1 = axes('position',[.03,.1,.29,.85]);
a2 = axes('position',[.36,.56,.29,.38]);
a3 = axes('position',[.69,.56,.29,.38]);
a4 = axes('position',[.69,.1,.29,.38]);
a5 = axes('position',[.36,.1,.29,.38]);
cost_back = [];

for epoch=1:nepochs,

    % to avoid small changes
    if ~mod(epoch,50) | mod(epoch,50) == 1
        netinit_during;
        mincost = 1E20;
    end

    doepoch_MND_sparse;
    cost_back = [cost_back cost];
    plotdata_mine
    epochs = epochs + 1;
    reportresults

    % extract the separation result
    for i=1:ninputs temp_y(i,:) = input3{i}(1,:); end;
    % extract the direct connections
    W_n = []; for i=1:ninputs W_n(i,:) = weight13{i}(1:ninputs); end; % diag(1./diag(W_n))*W_n,
    temp=[]; for i=1:ninputs temp(i,:) = weight23{i} * output2{i}; end; 
    % the nonlinear distortion level in the separation system. Plotted.
    fprintf('Nonlinear distortin level: '); disp(var(temp')./var(temp_y')),
    
    %% for causality discovery...
    % row permutation of W_n
    [Wp,rowp] = nzdiaglinprog( W_n );
    fprintf('Done!\n');
    % threshold: to set very small weights to zero
    Tol = 0.02;
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Divide each row of Wp by the diagonal element
    estdisturbancestd = 1./diag(abs(Wp));
    Wp = Wp./(diag(Wp)*ones(1,ninputs));
    Wp = (abs(Wp)>Tol) .* Wp;
    
    % Compute corresponding B
    Best = eye(ninputs)-Wp;
    
    % Check if B is lower-triangular after permutation
    [IF_causality, Bestcausal,causalperm] = check_causality( Best );
    if IF_causality
%         Best, 
        icausal = iperm(causalperm);
        B = Bestcausal(icausal, icausal);
        k = causalperm;
        % plot the result. In the simulation study this is not necessary.
        % plotmodel(B,k);
        yy = W_n * trpattern;
        mixeddata = trpattern; processdata;
        % the scatter plot of the separation results v.s. their linear
        % part.
        figure(9),
        for ii = 1:ninputs
            subplot(ceil(ninputs/2), 2,ii), plot(yy(ii,:), separdata(ii,:),'.'); 
        end
        %% the scatter plot of the separation results v.s. their nonlinear
        %% part
        % figure(10),
        % for ii = 1:ninputs
        %     subplot(ceil(ninputs/2), 2,ii), plot(yy(ii,:), separdata(ii,:)-yy(ii,:),'.'); 
        % end
        break;
    end
end

