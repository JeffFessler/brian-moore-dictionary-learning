function [L S NRMSE its rel] = LpS21_MRI(Y,model,opt,truth)
%--------------------------------------------------------------------------
% Syntax:       [L S] = LpS21_MRI(Y,model);
%               [L S] = LpS21_MRI(Y,model,opt);
%               [L S] = LpS21_MRI(Y,model,opt,truth);
%               [L S NRMSE its rel] = LpS21_MRI(Y,model);
%               [L S NRMSE its rel] = LpS21_MRI(Y,model,opt);
%               [L S NRMSE its rel] = LpS21_MRI(Y,model,opt,truth);
%               
% Inputs:             Y        : Observations                    (Defaults)
%    [OPTIONAL] model.A        : @() System matrix                      (I)
%    [OPTIONAL] model.At       : @() System matrix, transposed          (I)
%    [OPTIONAL] model.T        : @() Sparse transform                   (I)
%    [OPTIONAL] model.Tt       : @() Sparse transform, transposed       (I)
%               model.lambdaL  : Low-rank regularization parameter
%               model.r        : Explicit rank regularization
%               model.p        : Ell-p low-rank shrinkage parameter
%               model.lambdaS  : Sparse regularization parameter
%    [OPTIONAL] model.nb       : [nbx nby] block size           (No blocks)
%    [OPTIONAL] model.mask     : Observed entry mask          (true(ny,nx))
%    [OPTIONAL] model.M0       : Initial iterate                   (A' * Y)
%    [OPTIONAL]   opt.tau      : Proximal gradient step size            (1)
%    [OPTIONAL]   opt.tol      : Stopping tolerance                  (1e-6)
%    [OPTIONAL]   opt.Nmin     : Minimum # iterations                   (5)
%    [OPTIONAL]   opt.Nmax     : Maximum # iterations                 (100)
%    [OPTIONAL]   opt.method   : {'OptShrink','SVT','SVTp'}   ('OptShrink')
%    [OPTIONAL] truth.Xtrue    : Ground truth X = L + S                ([])
%    [OPTIONAL] truth.Ltrue    : Ground truth L                        ([])
%    [OPTIONAL] truth.Strue    : Ground truth S                        ([])
%               
% Outputs:      {L,S} are the reconstructed low-rank and sparse components
%               
%               NRMSE.X(i) contains (when truth.Xtrue is specified) the
%               NRMSE values of X = L + S w.r.t. Xtrue at iteration i
%               
%               NRMSE.L(i) contains (when truth.Ltrue is specified) the
%               NRMSE values of L w.r.t. Ltrue at iteration i
%               
%               NRMSE.S(i) contains (when truth.Strue is specified) the
%               NRMSE values of S w.r.t. Strue at iteration i
%               
%               its is the number of iterations required for convergence
%               
%               rel(i) contains the relative convergence values for M{:}
%               at iteration i
%               
% Description:  This function solves the optimization problem
%               
%{L,S} = argmin 0.5 ||Y - A(L + S)||^2 + lambdaL ||L||* + lambdaS ||TS||2,1
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
if isfield(model,'M0')
    % User-specified initial iterate
    M0 = model.M0;
else
    % Default initial iterate
    M0 = model.At(Y);
end

% Parse optional optimization parameters
if (~exist('opt','var') || ~isfield(opt,'tau'))
    % Default step size
    opt.tau = 1;
end
if ~isfield(opt,'tol')
    % Default stopping tolerance
    opt.tol = 1e-6;
end
if ~isfield(opt,'Nmin')
    % Default minimum # iterations
    opt.Nmin = 5;
end
if ~isfield(opt,'Nmax')
    % Default maximum # iterations
    opt.Nmax = 100;
end
if ~isfield(opt,'method')
    % Default low-rank update method
    opt.method = 'OptShrink';
end

% Parse ground truth parameters
if ~exist('truth','var')
    % Empty ground truth
    truth = struct();
end
computeNRMSE.X = (isfield(truth,'Xtrue') && ~isempty(truth.Xtrue));
computeNRMSE.L = (isfield(truth,'Ltrue') && ~isempty(truth.Ltrue));
computeNRMSE.S = (isfield(truth,'Strue') && ~isempty(truth.Strue));

%--------------------------------------------------------------------------
% Proximal gradient method (PGM)
%--------------------------------------------------------------------------
% Initializations
A = model.A;
At = model.At;
T = model.T;
Tt = model.Tt;
M = M0;
L = zeros(size(M));
S = zeros(size(M));

% NRMSE function
NRMSEfcn = @(X,Xt) norm(X(:) - Xt(:)) / norm(Xt(:));

% Iterations
its = 0;
rel = [];
NRMSE = struct('X',[],'L',[],'S',[]);
done = false;
while ~done
    % Start iteration timer
    tic;
    
    % Update iteration count
	its = its + 1;
    
	% Low-rank update
    Llast = L;
    switch lower(opt.method)
        case 'svt'
            if (useLLR == false)
                % Apply SVT
                lambdaL = opt.tau * model.lambdaL;
                L = SVT(M - S,lambdaL);
            else
                % Apply locally low-rank (LLR) SVT
                lambdaL = opt.tau * model.lambdaL;
                MmS = embed(M - S,model.mask);
                L = masker(SVT_LLR(MmS,lambdaL,model.nb),model.mask);
            end
        case 'svtp'
            if (useLLR == false)
                % Apply ell-p singular value shrinkage
                lambdaL = opt.tau * model.lambdaL;
                L = SVTp(M - S,model.p,lambdaL);
            else
                % Apply locally low-rank (LLR) ell-p singular value shrinkage
                lambdaL = opt.tau * model.lambdaL;
                MmS = embed(M - S,model.mask);
                L = masker(SVTp_LLR(MmS,model.p,lambdaL,model.nb),model.mask);
            end
        case 'optshrink'
            if (useLLR == false)
                % Apply OptShrink
                L = OptShrink(M - S,model.r);
            else
                % Apply locally low-rank (LLR) OptShrink
                MmS = embed(M - S,model.mask);
                L = masker(OptShrink_LLR(MmS,model.r,model.nb),model.mask);
            end
        otherwise
            % Method not recognized
            error('method "%s" not supported',opt.method);
    end
    
	% Sparse update: Mixed 2-1 norm
	%S = Tt(soft(T(M - Llast),opt.tau * model.lambdaS));
	S = Tt(mixed21(T(M - Llast),opt.tau * model.lambdaS));
    
	% Data consistency update
	Mlast = M;
	M = L + S - opt.tau * At(A(L + S) - Y);
    
    % Compute relative convergence
    dM = norm(M(:) - Mlast(:),2) / norm(Mlast(:),2);
    rel(its) = dM; %#ok
    
    % Compute NRMSEs, if necessary
    if (computeNRMSE.X == true)
        NRMSE.X(its) = NRMSEfcn(L + S,truth.Xtrue);
    end
    if (computeNRMSE.L == true)
        NRMSE.L(its) = NRMSEfcn(L,truth.Ltrue);
    end
    if (computeNRMSE.S == true)
        NRMSE.S(its) = NRMSEfcn(S,truth.Strue);
    end
    
    % Display progress
    fprintf('Iteration %i - dM = %.2e - Time = %.2fs\n',its,dM,toc);
    
    % Stopping criteria
    if ((its >= opt.Nmin) && ((its >= opt.Nmax) || (dM <= opt.tol)))
        % Finish updates
        done = true;
    end
end
%--------------------------------------------------------------------------
