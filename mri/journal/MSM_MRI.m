function [X NRMSE its rel] = MSM_MRI(Y,model,opt,GDopt,truth)
%--------------------------------------------------------------------------
% Syntax:       X = MSM_MRI(Y,model);
%               X = MSM_MRI(Y,model,opt);
%               X = MSM_MRI(Y,model,opt,GDopt);
%               X = MSM_MRI(Y,model,opt,GDopt,truth);
%               [X NRMSE its rel] = MSM_MRI(Y,model);
%               [X NRMSE its rel] = MSM_MRI(Y,model,opt);
%               [X NRMSE its rel] = MSM_MRI(Y,model,opt,GDopt);
%               [X NRMSE its rel] = MSM_MRI(Y,model,opt,GDopt,truth);
%               
% Inputs:       Y is an observations matrix                      (Defaults)
%               
%               model is a "model" structure with fields:
%               
%        [OPTIONAL] model.A     = @() System matrix                     (I)
%        [OPTIONAL] model.At    = @() System matrix, transposed         (I)
%        [OPTIONAL] model.mu    = Denoising parameter                   (1)
%                   model.X     = Model structure cell-array
%        [OPTIONAL] model.X0    = Initial iterate                  (A' * Y)
%               
%               where each element of model.X is itself a "regularizer"
%               structure cell-array with fields:
%               
%        [OPTIONAL] model.X{i}{j}.U      = @() L-transform              (I)
%        [OPTIONAL] model.X{i}{j}.Ut     = @() L-transform, transposed  (I)
%        [OPTIONAL] model.X{i}{j}.V      = @() R-transform              (I)
%        [OPTIONAL] model.X{i}{j}.Vt     = @() R-transform, transposed  (I)
%                   model.X{i}{j}.params = Proximal parameters
%                   model.X{i}{j}.prox   = Proximal function string
%                   
%                            f(X) |       prox      | params
%      ---------------------------+-----------------+---------------
%      (element-wise)     \|S\|_1 | 'soft'          | {lambda}
%      (element-wise)   \|S\|_p^p | 'softp'         | {p,lambda}
%      (column-wise)  \|S\|_{2,1} | 'mixed21'       | {lambda}
%      (svd)              \|L\|_* | 'SVT'           | {lambda}
%      (svd)            \|L\|_p^p | 'SVTp'          | {p,lambda}
%                           (RMT) | 'OptShrink'     | {r}
%                       \|X\|_F^2 | 'fro'           | {lambda}
%            I_{\|X-M\|_F <= eps} | 'ball2'         | {M,eps} 
%      (svd)   \sum_i \|P_i L\|_* | 'SVT_LLR'       | {[],nb,lambda}
%                                 |                 | {mask,nb,lambda}
%      (svd) \sum_i \|P_i L\|_p^p | 'SVTp_LLR'      | {[],nb,p,lambda}
%                                 |                 | {mask,nb,p,lambda}
%      (svd)     (block-wise RMT) | 'OptShrink_LLR' | {[],nb,r}
%                                 |                 | {mask,nb,r}
%               
%               opt is an "optimization" structure with fields:
%               
%        [OPTIONAL] opt.rho     = ADMM parameter                      (1e1)
%        [OPTIONAL] opt.epsilon = Stopping tolerance                 (1e-6)
%        [OPTIONAL] opt.Nmin    = Minimum # iterations                  (5)
%        [OPTIONAL] opt.Nmax    = Max # iterations                    (100)
%               
%               GDopt is a "gradient descent" structure with fields:
%               
%        [OPTIONAL] GDopt.k     = # gradient descent iterations         (5)
%        [OPTIONAL] GDopt.normA = Spectral norm of system matrix        (1)
%        [OPTIONAL] GDopt.tau   = Step size - range (0,2) - scaling     (1)
%               
%               truth is a "ground truth" structure with fields:
%               
%        [OPTIONAL] truth.Xtrue     = ground truth X1 + ... + Xm       ([])
%        [OPTIONAL] truth.Xitrue{i} = ground truth of ith component    ([])
%               
% Outputs:      X = {X1,...,Xm} contains the solution to (*)
%               
%               NRMSE.X(j) contains (when truth.Xtrue is specified) NRMSE
%               values of X = X1 + ... + Xm w.r.t. Xtrue at iteration j
%               
%               NRMSE.Xi{i}(j) contains (when truth.Xitrue{j} is specified)
%               NRMSE value of X{i} w.r.t. Xitrue{i} at iteration j
%               
%               its is the number of iterations required for convergence
%               
%               rel(j,:) contains the relative convergence values for X{:}
%               at iteration j
%               
% Description:  This function solves the structured matrix problem
%               
%               minimize 0.5 mu ||Y - A(X1 + ... Xm)||_F^2 +            (*)
%                        sum_ij lambda_ij f_ij(U_ij Xi V_ij')
%               
%               using alternating direction method of multipliers (ADMM)
%               
% Author:       Brian Moore
%               brimoor@umich.edu
%               
% Date:         September 7, 2014
%--------------------------------------------------------------------------

% Parse optional model parameters
if ~isfield(model,'A') || ~isfield(model,'At')
    % Default system matrix
    model.A = @(X) X;
    model.At = @(X) X;
end
if ~isfield(model,'mu')
    % Default denoising regularizer
    model.mu = 1;
end
m = length(model.X);
n = cellfun(@length,model.X);
for i = 1:m
    for j = 1:n(i)
        if (~isfield(model.X{i}{j},'U') || ~isfield(model.X{i}{j},'Ut'))
            % Default (no transformation)
            model.X{i}{j}.U = @(X) X;
            model.X{i}{j}.Ut = @(X) X;
        end
        if (~isfield(model.X{i}{j},'V') || ~isfield(model.X{i}{j},'Vt'))
            % Default (no transformation)
            model.X{i}{j}.V = @(X) X;
            model.X{i}{j}.Vt = @(X) X;
        end
    end
end
if isfield(model,'X0')
    % User-specified initial iterate
    X0 = model.X0 / m;
else
    % Default initial iterate
    X0 = model.At(Y) / m;
end

% Parse optional ADMM parameters
if (~exist('opt','var') || ~isfield(opt,'rho'))
    % Default ADMM step size
    opt.rho = 1e1;
end
if ~isfield(opt,'epsilon')
    % Default stopping tolerance
    opt.epsilon = 1e-6;
end
if ~isfield(opt,'Nmin')
    % Default min # iterations
    opt.Nmin = 5;
end
if ~isfield(opt,'Nmax')
    % Default max # iterations
    opt.Nmax = 100;
end

% Parse optional GD parameters
if (~exist('GDopt','var') || ~isfield(GDopt,'k'))
    % Default # GD iterations
    GDopt.k = 5;
end
if ~isfield(GDopt,'normA')
    % Default system matrix norm
    GDopt.normA = 1;
end
if ~isfield(GDopt,'tau')
    % Default step size scaling
    GDopt.tau = 1;
end

% Parse ground truth parameters
if ~exist('truth','var')
    % Empty ground truth
    truth = struct();
end
computeNRMSE.X = (isfield(truth,'Xtrue') && ~isempty(truth.Xtrue));
for i = 1:m
    computeNRMSE.Xi{i} = (isfield(truth,'Xitrue') && ...
                         (length(truth.Xitrue) >= i) && ...
                         ~isempty(truth.Xitrue{i}));
end

%--------------------------------------------------------------------------
% Alternating direction method of multiplieres (ADMM)
%--------------------------------------------------------------------------
% Parameter initializations
A = model.A;
At = model.At;
mu = model.mu;
rho = opt.rho;
epsilon = opt.epsilon;
Nmin = opt.Nmin;
Nmax = opt.Nmax;
k = GDopt.k;
tau = GDopt.tau / (m * mu * GDopt.normA^2 + max(n) * rho);

% NRMSE function
NRMSEfcn = @(X,Xt) norm(X(:) - Xt(:)) / norm(Xt(:));

% Variable initializations
grad = cell(m,1);
X = cell(m,1);
Xs = cell(m,1);
Lambda = cell(m,1);
for i = 1:m
    X{i} = X0;
    Xs{i} = cell(n(i),1);
    Lambda{i} = cell(n(i),1);
    for j = 1:n(i)
        Xs{i}{j} = X0;
        Lambda{i}{j} = 0;
    end
end

% Iterations
its = 0;
rel = zeros(0,m);
NRMSE = struct('X',[],'Xi',{cell(1,m)});
maxRel = inf;
while ((its < Nmin) || ((its < Nmax) && (maxRel > epsilon)))
    % Start iteration timer
    tic;
    
    % Update iteration count
    its = its + 1;
    
    % X-updates
    Xlast = X;
    for i = 1:m
        % Gradient function handle
        grad{i} = @(X) mu * At(A(sum(cat(3,X{:}),3)) - Y) + ...
                      rho * (n(i) * X{i} - sum(cat(3,Xs{i}{:}),3) + sum(cat(3,Lambda{i}{:}),3));
    end
    %X = GD(grad,X,tau,k);
    X = nGD(grad,X,tau,k);
    
    % Xs-updates
    for i = 1:m
        for j = 1:n(i)
            % Load proximal paramaters
            U = model.X{i}{j}.U;
            Ut = model.X{i}{j}.Ut;
            V = model.X{i}{j}.V;
            Vt = model.X{i}{j}.Vt;
            params = model.X{i}{j}.params;
            str = lower(model.X{i}{j}.prox);
            if ~any(ismember({'optshrink','optshrink_llr','ball2'},str))
                % Scale regularization by ADMM parameter
                params{end} = params{end} / rho;
            end
            
            % Apply proximal operator
            Xs{i}{j} = Ut(V(prox(Vt(U(X{i} + Lambda{i}{j})),str,params)));
        end
    end
    
    % Lambda-updates
    for i = 1:m
        for j = 1:n(i)
            Lambda{i}{j} = Lambda{i}{j} + X{i} - Xs{i}{j};
        end
    end
    
    % Update relative changes
    rel(its,:) = cellfun(@(Xi,Xli) NRMSEfcn(Xi,Xli),X,Xlast);
    maxRel = max(rel(its,:));
    
    % Compute NRMSEs, if necessary
    if (computeNRMSE.X == true)
        NRMSE.X(its) = NRMSEfcn(sum(cat(3,X{:}),3),truth.Xtrue);
    end
    for i = 1:m
        if (computeNRMSE.Xi{i} == true)
            NRMSE.Xi{i}(its) = NRMSEfcn(X{i},truth.Xitrue{i});
        end
    end
    
    % Display progress
    fprintf('Iteration %i - max(dX) = %.2e - Time = %.2fs\n',its,maxRel,toc);
end
%--------------------------------------------------------------------------

% Average over variable splits
for i = 1:m
    X{i} = (1 / (n(i) + 1)) * (X{i} + sum(cat(3,Xs{i}{:}),3));
end

end

%
% Gradient descent
%
function Xk = GD(grad,X0,tau,its)
% Syntax:   Xk = GD(grad,X0,tau,its);
    
% Gradient descent
Xk = X0;
for i = 1:its
    % X-updates
    fcn = @(X,grad) X - tau * grad(Xk);
    Xk = cellfun(fcn,Xk,grad,'UniformOutput',false);
end

end

%
% Nesterov-accelerated gradient descent
%
function Xk = nGD(grad,X0,tau,its)
% Syntax:   Xk = nGD(grad,X0,tau,its);

% Nesterov-accelerated gradient descent
tk = 0;
Xlast = X0;
Xk = X0;
for i = 1:its
    % t-update
    tlast = tk;
    tk = 0.5 * (1 + sqrt(1 + 4 * tk^2));

    % Z-updates
    Zfcn = @(X,Xlast) X + ((tlast - 1) / tk) * (X - Xlast);
    Zk = cellfun(Zfcn,Xk,Xlast,'UniformOutput',false);

    % X-updates
    Xlast = Xk;
    Xfcn = @(Z,grad) Z - tau * grad(Zk);
    Xk = cellfun(Xfcn,Zk,grad,'UniformOutput',false);        
end

end

%
% Apply proximal operator
%
function Yhat = prox(Y,str,params)
% Syntax:   Yhat = prox(Y,str,params);

% Apply proximal operator
switch lower(str)
    case 'soft'
        % Apply soft-thresholding
        Yhat = soft(Y,params{:});
    case 'softp'
        % Apply ell-p shrinkage
        Yhat = pShrink(Y,params{:});
    case 'mixed21'
        % Apply mixed 2-1 norm shrinkage
        Yhat = mixed21(Y,params{:});
    case 'svt'
        % Apply singular value thresholding
        Yhat = SVT(Y,params{:});
    case 'svt_llr'
        % Apply locally low-rank (LLR) singular value thresholding
        mask = params{1};
        params = {params{3:end} params{2}};
        if isempty(mask)
            % No mask
            Yhat = SVT_LLR(Y,params{:});
        else
            % Embed/remask data
            Yhat = masker(SVT_LLR(embed(Y,mask),params{:}),mask);
        end
    case 'svtp'
        % Apply ell-p singular value shrinkage
        Yhat = SVTp(Y,params{:});
    case 'svtp_llr'
        % Apply locally low-rank (LLR) ell-p singular value shrinkage
        mask = params{1};
        params = {params{3:end} params{2}};
        if isempty(mask)
            % No mask
            Yhat = SVTp_LLR(Y,params{:});
        else
            % Embed/remask data
            Yhat = masker(SVTp_LLR(embed(Y,mask),params{:}),mask);
        end
    case 'optshrink'
        % Apply OptShrink
        Yhat = OptShrink(Y,params{:});
    case 'optshrink_llr'
        % Apply locally low-rank (LLR) OptShrink
        mask = params{1};
        params = {params{3:end} params{2}};
        if isempty(mask)
            % No mask
            Yhat = OptShrink_LLR(Y,params{:});
        else
            % Embed/remask data
            Yhat = masker(OptShrink_LLR(embed(Y,mask),params{:}),mask);
        end
    case 'fro'
        % Apply frobenius shrinkage
        Yhat = FroShrink(Y,params{:});
    case 'ball2'
        % Apply Euclidean ball projection
        Yhat = EuclideanBallProject(Y,params{:});
end

end
