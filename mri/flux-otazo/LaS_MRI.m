function [X NRMSE rel] = LaS_MRI(Y,model,opt,GDopt,Xtrue)
%--------------------------------------------------------------------------
% Syntax:       X = LaS_MRI(Y,model);
%               X = LaS_MRI(Y,model,opt);
%               X = LaS_MRI(Y,model,opt,GDopt);
%               X = LaS_MRI(Y,model,opt,GDopt,Xtrue);
%               [X NRMSE its rel] = LaS_MRI(Y,model);
%               [X NRMSE its rel] = LaS_MRI(Y,model,opt);
%               [X NRMSE its rel] = LaS_MRI(Y,model,opt,GDopt);
%               [X NRMSE its rel] = LaS_MRI(Y,model,opt,GDopt,Xtrue);
%               
% Inputs:             Y        : Observations                    (Defaults)
%    [OPTIONAL] model.A        : @() System matrix                      (I)
%    [OPTIONAL] model.At       : @() System matrix, transposed          (I)
%    [OPTIONAL] model.T        : @() Sparse transform                   (I)
%    [OPTIONAL] model.Tt       : @() Sparse transform, transposed       (I)
%               model.lambdaL  : Low-rank regularization parameter
%               model.r        : Sparse regularization parameter
%               model.p        : Ell-p low-rank shrinkage parameter
%               model.lambdaS  : Sparse regularization parameter
%    [OPTIONAL] model.nb       : [nbx nby] block size           (No blocks)
%    [OPTIONAL] model.mask     : Observed entry mask          (true(ny,nx))
%    [OPTIONAL] model.X0       : Initial iterate                   (A' * Y)
%    [OPTIONAL]   opt.beta0    : Initial beta                         (1e6)
%    [OPTIONAL]   opt.betaInc  : Beta increment                        (25)
%    [OPTIONAL]   opt.Ninner   : # inner iterations                    (10)
%    [OPTIONAL]   opt.Nouter   : # outer iterations                     (5)
%    [OPTIONAL]   opt.methodL  : {'OptShrink','PCA', ...
%                                 'SVT','SVTp'}               ('OptShrink')
%    [OPTIONAL]   opt.methodS  : {'soft','mixed21'}                ('soft')
%    [OPTIONAL] GDopt.k        : # gradient descent iterations          (5)
%    [OPTIONAL] GDopt.normA    : Spectral norm of system matrix         (1)
%    [OPTIONAL] GDopt.tau      : Step size - range (0,2) - scaling      (1)
%    [OPTIONAL]       Xtrue    : Ground truth matrix                   ([])
%               
% Outputs:      X is the reconstructed frame sequence
%               
%               NRMSE(i,j) contains (when Xtrue is specified) the NRMSE
%               values for X w.r.t. Xtrue at (inner,outer) iteration (i,j)
%               
%               rel(i,j) contains the relative convergence value for X at
%               (inner,outer) iteration (i,j)
%               
% Description:  This function solves the optimization problem
%               
%     X = argmin 0.5 ||Y - A(X)||^2 + lambdaL ||X||* + lambdaS ||TX||1     
%               
%               using OptShrink(), SVT(), or SVTp() for the low-rank
%               updates
%               
% Author:       Brian Moore
%               brimoor@umich.edu
%               
% Date:         September 1, 2014
%--------------------------------------------------------------------------

% Parse optional model parameters
if ~isfield(model,'A') || ~isfield(model,'At')
    % Default system matrix
    model.A = @(X) X;
    model.At = @(X) X;
end
if (~isfield(model,'T') || ~isfield(model,'Tt'))
    % Default sparse transformation
    model.T = @(X) X;
    model.Tt = @(X) X;
end
if isfield(model,'nb')
    % Use locally low-rank model
    useLLR = true;
    if ~isfield(model,'mask')
        error('mask required for LLR computation');
    end
else
    % No locally low-rank model
    useLLR = false;
end
if isfield(model,'X0')
    % User-specified initial iterate
    X0 = model.X0;
else
    % Default initial iterate
    X0 = model.At(Y);
end

% Parse optional optimization parameters
if (~exist('opt','var') || ~isfield(opt,'beta0'))
    % Default initial beta
    opt.beta0 = 1e6;
end
if ~isfield(opt,'betaInc')
    % Default beta increment
    opt.betaInc = 25;
end
if ~isfield(opt,'Ninner')
    % Default # inner iterations
    opt.Ninner = 10;
end
if ~isfield(opt,'Nouter')
    % Default # outer iterations
    opt.Nouter = 5;
end
if ~isfield(opt,'methodL')
    % Default low-rank update method
    opt.methodL = 'OptShrink';
end
if ~isfield(opt,'methodS')
    % Default sparse update method
    opt.methodS = 'soft';
end

% Parse optional gradient descent parameters
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
if (~exist('Xtrue','var') || isempty(Xtrue))
    % Don't compute NRMSE values
    computeNRMSE = false;
else
    % *Do* compute NRMSE values
    computeNRMSE = true;
end

%--------------------------------------------------------------------------
% Augmented Lagrangian method (ALM)
%--------------------------------------------------------------------------
% Initializations
A = model.A;
At = model.At;
T = model.T;
Tt = model.Tt;
X = X0;

% NRMSE function
NRMSEfcn = @(X,Xt) norm(X(:) - Xt(:)) / norm(Xt(:));

% Outer iterations
rel = [];
NRMSE = [];
beta = opt.beta0;
for j = 1:opt.Nouter
    % Inner iterations
    for i = 1:opt.Ninner
        % Start iteration timer
        tic;
        
        % L-update
        switch lower(opt.methodL)
            case 'svt'
                if (useLLR == false)
                    % Apply SVT
                    lambdaL = model.lambdaL / beta;
                    L = SVT(X,lambdaL);
                else
                    % Apply locally low-rank (LLR) SVT
                    lambdaL = model.lambdaL / beta;
                    L = masker(SVT_LLR(embed(X,model.mask),lambdaL,model.nb),model.mask);
                end
            case 'svtp'
                if (useLLR == false)
                    % Apply ell-p singular value shrinkage
                    lambdaL = model.lambdaL / beta;
                    L = SVTp(X,model.p,lambdaL);
                else
                    % Apply locally low-rank (LLR) ell-p singular value shrinkage
                    lambdaL = model.lambdaL / beta;
                    L = masker(SVTp_LLR(embed(X,model.mask),model.p,lambdaL,model.nb),model.mask);
                end
            case 'optshrink'
                if (useLLR == false)
                    % Apply OptShrink
                    L = OptShrink(X,model.r);
                else
                    % Apply locally low-rank (LLR) OptShrink
                    L = masker(OptShrink_LLR(embed(X,model.mask),model.r,model.nb),model.mask);
                end
            case 'pca'
                if (useLLR == false)
                    % Apply PCA
                    L = PCA(X,model.r);
                else
                    % Apply locally low-rank (LLR) PCA
                    L = masker(PCA_LLR(embed(X,model.mask),model.r,model.nb),model.mask);
                end
            otherwise
                % Method not recognized
                error('method "%s" not supported',opt.methodL);
        end
        
        % S-update
        switch lower(opt.methodS)
            case 'soft'
                % Apply soft thresholding
                S = soft(T(X),model.lambdaS / beta);
            case 'mixed21'
                % Apply mixed 2-1 prox function
                S = mixed21(T(X),model.lambdaS / beta);
            otherwise
                % Method not recognized
                error('method "%s" not supported',opt.methodS);
        end
        
        % X-update
        Xlast = X;
        grad = @(X) At(A(X) - Y) + beta * (2 * X - L - Tt(S));
        tau = GDopt.tau / (GDopt.normA^2 + 2 * beta);
        %X = GD(grad,X,tau,GDopt.k);
        X = nGD(grad,X,tau,GDopt.k);
        
        % Update relative changes
        rel(i,j) = NRMSEfcn(X,Xlast); %#ok
        
        % Compute NRMSEs, if necessary
        if (computeNRMSE == true)
            NRMSE(i,j) = NRMSEfcn(X,Xtrue); %#ok
        end
        
        % Display progress
        fprintf('Iteration (%i,%i)/(%i,%i) - Time = %.2fs\n',i,j,opt.Ninner,opt.Nouter,toc);
    end
    
    % Update beta
    beta = beta * opt.betaInc;
end
%--------------------------------------------------------------------------

end

%
% Gradient descent
%
function Xk = GD(grad,X0,tau,its)
% Syntax:   Xk = GD(grad,X0,tau,its);

% Gradient descent
Xk = X0;
for i = 1:its
    % X-update
    Xk = X - tau * grad(Xk);
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
    
    % Z-update    
    Zk = Xk + ((tlast - 1) / tk) * (Xk - Xlast);
    
    % X-update
    Xlast = Xk;
    Xk = Zk - tau * grad(Zk);
end

end
