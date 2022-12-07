%%
% Compute norm(A) for G2S models
%
% Brian Moore
% brimoor@umich.edu
%

rng(1);

% Knobs
mList = 2:3:50;
nList = 2:3:50;
type  = '1c'; % {'1c','1n','2a','2b','2c'}

% Compute norms
nm = numel(mList);
nn = numel(nList);
normA = nan(nm,nn);
PT = ProgressTimer(nm * nn,10);
for i = 1:nm
    for j = 1:nn
        % Generate model
        m = mList(i);
        n = nList(j);
        M = true(m,n);
        N = randn(m,n,3); 
        
        % Compute norm(A)
        [A, ~, ~]  = G2Smodel(N,M,type);
        normA(i,j) = svds(A,1);
        
        % Update progress
        PT.Update();
    end
end

% Plot results
figure;
surf(nList,mList,normA);
title(sprintf('max(normA(:)) = %.5f',max(normA(:))));
xlabel('n');
ylabel('m');
