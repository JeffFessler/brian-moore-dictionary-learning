function [D, Bt, params, stats] = onlineUnitaryDil(Yt,mu,gamma,varargin)
%
% Syntax:       [D, Bt] = onlineUnitaryDil(Yt,mu,gamma);
%               [D, Bt] = onlineUnitaryDil(Yt,mu,gamma,opts);
%               [D, Bt, stats] = onlineUnitaryDil(Yt,mu,gamma);
%               [D, Bt, stats] = onlineUnitaryDil(Yt,mu,gamma,opts);
%               
% Inputs:       Yt is a d x n data matrix
%               
%               mu >= 0 is the sparsity regularization parameter
%               
%               gamma in [0, 1] is the forgetting factor, where gamma = 0
%               is memoryless, and gamma = 1 is infinite memory
%               
%               [OPTIONAL] opts is a struct containing one or more of the
%               following fields. The default values are in ()
%                   
%                   opts.D0 (dctmtx(d)') is a d x m unitary matrix
%                   containing the initial dictionary for time t (usually
%                   the dictionary from time t - 1). Undercomplete
%                   dictionaries (m < d) are allowed, in which case the
%                   columns of the initial dictionary must be orthonormal
%                   
%                   opts.B0 (zeros(m,n)) is an m x n matrix containing the
%                   initial sparse codes for time t
%                   
%                   opts.params ([]) is a struct containing the cumulative
%                   parameters returned by the previous call to this 
%                   function. At time t = 1, omit this field or pass
%                   params = [] and it will be initialized appropriately
%                   
%                   opts.type ('hard') can be 'hard' or 'soft' and
%                   specifies whether apply hard or soft thresholding to Bt
%                   
%                   opts.fixedD (false) determines whether to fix the
%                   initial dictionary (true) or learn it (false)
%                   
%                   opts.nIters (50) is the number of block coordinate
%                   descent iterations to perform
%                   
%                   opts.flag (1) determines what status updates to print
%                   to the command window. The choices are
%                   
%                       flag = 0: no printing
%                       flag = 1: print iteration status
%               
% Outputs:      D is a d x m unitary matrix containing the dictionary at
%               time t
%               
%               Bt is an m x n matrix of sparse codes for time t
%               
%               params is a struct of cumulative parameters for time t (to
%               be passed back into this function at time t + 1)
%               
%               stats is a struct containing the following fields:
%               
%                   stats.nIters is the number of iterations performed
%                   
%                   stats.cost is a 1 x nIters vector containing the cost
%                   function values at each iteration
%                   
%                   stats.ccost is a 1 x nIters vector containing the full
%                   (cumulative) cost function values at each iteration 
%                   
%                   stats.repError is a 1 x nIters vector containing the
%                   representation error at each iteration, defined as
%                   \|Yt - D_k Bt_k\|_F / \|Yt\|_F
%                   
%                   stats.deltaD is a 1 x nIters vector containing the
%                   relative convergence of D at each iteration, defined as
%                   \|D_k - D_{k - 1}\|_F / \|D_{k - 1}\|_F
%                   
%                   stats.deltaB is a 1 x nIters vector containing the
%                   relative convergence of Bt at each iteration, defined:
%                   \|Bt_k - Bt_{k - 1}\|_F / \|Bt_{k - 1}\|_F
%                   
%                   stats.sparsity is a 1 x nIters vector containing the
%                   sparsity, in percent, of Bt at each iteration
%                   
%                   stats.time is a 1 x nIters vector containing the time,
%                   in seconds, required to perform each outer iteration
%               
% Description:  Performs dictionary learning via the Online Unitary
%               DIctionary Learning (Online Unitary-DIL) problem. That is,
%               when type == 'hard':
%               
%               min_{D,Bt} 0.5 \sum_{j=1}^t gamma^{t-j} \|Yj - D Bj\|_F^2 +
%                          0.5 \mu^2 \|Bt\|_0
%               
%               subject to D' * D == I
%               
%               Or, when type == 'soft':
%               
%               min_{D,Bt} 0.5 \sum_{j=1}^t gamma^{t-j} \|Yj - D Bj\|_F^2 +
%                          \mu \|Bt\|_1
%               
%               subject to D' * D == I
%               
% Author:       Brian Moore
%               brimoor@umich.edu
%               
% Date:         May 15, 2017
%               July 20, 2017
%               August 29, 2017
%

% Parse inputs
[D, Bt, params, type, fixedD, nIters, flag] = parseInputs(Yt,varargin{:});
PRINT_STATS     = (flag > 0);
COMPUTE_STATS   = PRINT_STATS || (nargout == 4);
UPDATE_PARAMS   = (nargout >= 3);
USE_HARD_THRESH = strcmpi(type,'hard');
OPTIMIZE_DICT   = ~fixedD;

% Initialize params/stats
if UPDATE_PARAMS || COMPUTE_STATS
    if strcmpi(type,'hard')
        % Ell-0 regularization
        Btfit = @(Bt) 0.5 * mu^2 * nnz(Bt);
    elseif strcmpi(type,'soft')
        % Ell-1 regularization
        Btfit = @(Bt) mu * norm(vec(Bt),1);
    else
        % Unsupported shrinkage
        error('Unsupported shrinkage type ''%s''',type);
    end
end
if COMPUTE_STATS
    % Cost function
    Ytfit = @(D,Bt) 0.5 * norm(vec(Yt - D * Bt))^2;
    Rtfit = @(D) 0.5 * real(trace(params.Et) + ...
                              trace(D' * (D * params.Gt - 2 * params.Ft)));
    % Cost function (only portion affected by D, Bt)
    Psi = @(D,Bt) Ytfit(D,Bt) + gamma * Rtfit(D) + Btfit(Bt);
    % Cumulative cost function
    CPsi = @(D,Bt) (Psi(D,Bt) + ...
                              gamma * params.bt) / (1 + gamma * params.ct);
    
    % Stats-printing function
    iterFmt = sprintf('%%0%dd',ceil(log10(nIters + 1)));
    out     = printFcn('Iter'    ,iterFmt,'ccost'   ,'%.2f'  , ...
                       'repError','%.3f' ,'deltaD'  ,'%.3e'  , ...
                       'deltaB'  ,'%.3e' ,'sparsity','%.2f%%', ...
                       'time'    ,'%.2fs');
    
    % Initialize stats
    cost     = nan(1,nIters);
    ccost    = nan(1,nIters);
    repError = nan(1,nIters);
    deltaD   = nan(1,nIters);
    deltaB   = nan(1,nIters);
    sparsity = nan(1,nIters);
    time     = nan(1,nIters);
end

% Online Unitary-DIL
if PRINT_STATS
    fprintf('***** Online Unitary-DIL *****\n');
end
for it = 1:nIters
    % Initialize iteration
    itimer = tic;
    if COMPUTE_STATS
        % Save last iterates
        Btlast = Bt;
        Dlast = D;
    end
    
    % Bt-update
    if USE_HARD_THRESH
        % Hard thresholding
        Bt = hard(D' * Yt,mu);
    else
        % Soft thresholding
        Bt = soft(D' * Yt,mu);
    end
    
    % D-update
    if OPTIMIZE_DICT
        D = onlineProcrustes(Yt,Bt,params.Ft,gamma);
    end
    
    % Record stats
    if COMPUTE_STATS
        cost(it)     =  Psi(D,Bt);
        ccost(it)    = CPsi(D,Bt);
        repError(it) = computeNRMSE(D * Bt,Yt);
        deltaD(it)   = computeNRMSE(D,Dlast);
        deltaB(it)   = computeNRMSE(Bt,Btlast);
        sparsity(it) = 100 * (nnz(Bt) / numel(Bt));
        time(it)     = toc(itimer);
        if PRINT_STATS
            out(it,ccost(it),repError(it),deltaD(it),deltaB(it), ...
                sparsity(it),time(it));
        end
    end
end

% Update cumulative parameters
if UPDATE_PARAMS
    params.t  = params.t + 1;
    params.Et = gamma * params.Et + Yt * Yt';
    params.Ft = gamma * params.Ft + Yt * Bt';
    params.Gt = gamma * params.Gt + Bt * Bt';
    params.bt = gamma * params.bt + Btfit(Bt);
    params.ct = gamma * params.ct + 1;
end

% Return stats
if COMPUTE_STATS
    stats.nIters   = nIters;
    stats.cost     = cost;
    stats.ccost    = ccost;
    stats.repError = repError;
    stats.deltaD   = deltaD;
    stats.deltaB   = deltaB;
    stats.sparsity = sparsity;
    stats.time     = time;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Local functions
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Hard thresholding
function Y = hard(X,mu)
Y = X .* (abs(X) > mu);

% Soft thresholding
function Y = soft(X,mu)
Y = sign(X) .* max(abs(X) - mu,0);

% Online Procrustes solution to
% \min_{Q: Q'Q = I} \sum_{j=1}^t \gamma^{t-j}\|Xj - Q Yj\|_F^2
% where
% Ft = \sum_{j=1}^{t-1} \gamma^{t-1-j} Xj Yj'
function Q = onlineProcrustes(Xt,Yt,Ft,gamma)
% NOTE: works when Q is d x m with m < d
[U, ~, V] = svd(gamma * Ft + Xt * Yt','econ');
Q = U * V';

% Print function
function out = printFcn(varargin)
str = [sprintf('%s[%s] ',varargin{:}), '\n'];
out = @(varargin) fprintf(str,varargin{:});

% Compute NRMSE
function err = computeNRMSE(Xhat,X)
denom = norm(X(:));
if isnan(denom)
    err = nan;
elseif denom == 0
    err = 0;
else
    err = norm(Xhat(:) - X(:)) / denom;
end

% Vectorize input
function x = vec(X)
x = X(:);

% Parse inputs
function [D0, B0, params, type, fixedD, nIters, flag] = ...
                                                       parseInputs(Yt,opts)
if ~exist('opts','var') || isempty(opts)
    opts = struct();
end

% Standard args
D0     = parseField(opts,'D0',nan);
B0     = parseField(opts,'B0',nan);
params = parseField(opts,'params',[]);
type   = parseField(opts,'type','hard');
fixedD = parseField(opts,'fixedD',false);
nIters = parseField(opts,'nIters',50);
flag   = parseField(opts,'flag',1);

% Expensive defaults
if isnan(D0), D0 = dctmtx(size(Yt,1))';             end
if isnan(B0), B0 = zeros(size(D0,2),size(Yt,2));    end

% Parse params
if ~isstruct(params)
    params.t  = 0;
    params.Et = zeros(size(D0,1));
    params.Ft = zeros(size(D0));
    params.Gt = zeros(size(D0,2));
    params.bt = 0;
    params.ct = 0; 
end

% Parse struct field
function val = parseField(stats,field,default)
if isfield(stats,field)
    val = stats.(field);
else
    val = default;
end
