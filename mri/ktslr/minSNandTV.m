function [U, cost, opts] = minSNandTV(A,At,D,Dt,U,b,~,opts)
% Syntax:   [U, cost, opts] = minSNandTV(A,At,D,Dt,U,b,~,opts);

% Parse inputs
[m, n, d] = size(U); 

% k-t SLR
fprintf('***** k-t SLR *****\n');
Lam1 = zeros(m,n,d);
Lam2 = zeros(m,n,d);
Lam3 = zeros(m,n,d);
Lam4 = zeros(m,n,d);
Dfwd = D(U);
cost = [];
for out = single(1:opts.outer_iter)
    for in = single(1:opts.inner_iter)
        % Initialize iteration
        itimer = tic;
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % W update (TV)
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        Z1    = Dfwd{1} + Lam1 / opts.beta2;
        Z2    = Dfwd{2} + Lam2 / opts.beta2;
        Z3    = Dfwd{3} + Lam3 / opts.beta2;
        V     = sqrt(abs(Z1).^2 + abs(Z2).^2 + abs(Z3).^2);
        V(~V) = 1;
        V     = max(0,V - 1 / opts.beta2) ./ V;
        Wx    = Z1 .* V;
        Wy    = Z2 .* V;
        Wt    = Z3 .* V;
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Lambda update (matrix ell-p value shrinkage)
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        UL4    = reshape(U + Lam4 / opts.beta1,[],d);
        Lambda = reshape(ellpshrink(UL4,opts.p,opts.beta1),m,n,d);
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % U update (CG)
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%    
        U    = CG_solver(b,A,At,D,Dt,Lambda,Lam1,Lam2,Lam3,Lam4,Wx,Wy,Wt,opts,U,1e-7,10);
        Dfwd = D(U);
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Lagrange multiplier updates
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        Lam1 = Lam1 - 1.618 * opts.beta2 * (Wx - Dfwd{1});
        Lam2 = Lam2 - 1.618 * opts.beta2 * (Wy - Dfwd{2});
        Lam3 = Lam3 - 1.618 * opts.beta2 * (Wt - Dfwd{3});
        Lam4 = Lam4 - 1.618 * opts.beta1 * (Lambda - U);
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Compute cost
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        cost(end + 1) = norm(vec(A(U) - b))^2 + ...
                        opts.mu1 * ellpnorm(reshape(U,[],d),opts.p) + ...
                        opts.mu2 * sum(vec(abs(sqrt(Dfwd{1}.^2 + Dfwd{2}.^2 + Dfwd{3}.^2)))); %#ok
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Display progress
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        fprintf('Iter: %d:%d/%d:%d, cost: %.3e, time: %.2fs\n',out,in,opts.outer_iter,opts.inner_iter,cost(end),toc(itimer));
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Check convergence
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        if (in > 1) && (abs(cost(end) - cost(end - 1)) < 1e-3 * abs(cost(end)))
            disp('Converged!');
            break;
        end
    end 
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Beta updates
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    opts.beta1 = opts.beta1 * opts.beta1rate;
    opts.beta2 = opts.beta2 * opts.beta2rate;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Local functions
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Vectorize data
function x = vec(X)
x = X(:);

% Matrix ell-p norm
function np = ellpnorm(X,p)
[~, Su, ~] = givefastSVD(X);
np         = sum(diag(Su).^p) / p;

% Matrix ell-p shrinkage
function Xp = ellpshrink(X,p,beta)
[U, S, V] = givefastSVD(X);
s         = diag(S);
Xp        = U * diag(max(0,s - (s.^(p - 1)) / beta)) * V';
