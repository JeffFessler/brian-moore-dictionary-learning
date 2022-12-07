function [B, stats] = sparseCoding(P,D,mu,varargin)
%
% Syntax:       B = sparseCoding(P,D,mu);
%               B = sparseCoding(P,D,mu,opts);
%               [B, stats] = sparseCoding(P,D,mu);
%               [B, stats] = sparseCoding(P,D,mu,opts);
%               
% Inputs:       P is an d x n data matrix
%               
%               D is a d x m data matrix
%               
%               mu >= 0 is the patch regularization parameter
%               
%               [OPTIONAL] opts is a struct containing one or more of the
%               following fields. The default values are in ()
%                   
%                   opts.nIters (10) is the number of iterations to
%                   perform
%                   
%                   opts.type ('hard') can be {'hard','soft'} and specifies
%                   whether to use hard or soft thresholding to regularize
%                   the sparse codes B
%                   
%                   opts.B0 (D' * P) is a m x n matrix specifying the
%                   initial B iterate
%                   
%                   opts.Btrue (nan) is the ground truth B matrix to use 
%                   for MSE calculations
%                   
%                   opts.NRMSEfcn (@computeNRMSE) is the function to use
%                   when computing the NRMSE of B after each iteration
%                   
%                   opts.accel (true) specifies whether to use accelerated
%                   proximal gradient steps
%                   
%                   opts.tau ((0.99 + ~accel) / norm(D)^2) is the step size
%                   parameter, and should satisfy
%                   
%                       tau <= 1 / norm(D)^2, when accel = true
%                       tau <  2 / norm(D)^2, when accel = false
%                   
%                   opts.flag (1) determines what status updates to print
%                   to the command window. The choices are
%                   
%                       flag = 0: no printing
%                       flag = 1: print iteration status
%               
% Outputs:      B is the estimated m x n matrix of sparse codes
%               
%               stats is a struct containing the following fields:
%               
%                   stats.nIters is the number of iterations performed
%                   
%                   stats.cost is a 1 x nIters vector containing the cost
%                   function values at each iteration
%                   
%                   stats.nrmse is a 1 x nIters vector containing the NRMSE
%                   of B with respect to Btrue at each iteration
%                   
%                   stats.delta is a 1 x nIters vector containing the
%                   relative convergence of B at each iteration, defined as
%                   \|B_{k + 1} - B_k\|_F / \|B_k\|_F
%                   
%                   stats.time is a 1 x nIters vector containing the time,
%                   in seconds, required to perform each outer iteration
%               
% Description:  Solves one of the following sparse coding problems:
%               
%               When type == 'hard':
%               
%               \min_{B} 0.5 \|P - DB\|_F^2 + 0.5 \mu^2 \|B\|_0
%               
%               When type == 'soft':
%               
%               \min_{B} 0.5 \|P - DB\|_F^2 + \mu \|B\|_1
%               
% Author:       Brian Moore
%               brimoor@umich.edu
%               
% Date:         January 31, 2017
%

% Parse inputs
[nIters, type, B, Btrue, NRMSEfcn, accel, tau, flag] = ...
                                              parseInputs(P,D,varargin{:});
PRINT_STATS    = (flag > 0);
COMPUTE_STATS  = PRINT_STATS || (nargout == 2);

% Parse type
if strcmpi(type,'hard');
    shrink = @(Y,tau) hard(Y,sqrt(tau) * mu);
elseif strcmpi(type,'soft')
    shrink = @(Y,tau) soft(Y,tau * mu);
else
    % Unsupported shrinkage
    error('Unsupported shrinkage type ''%s''',type);
end

% Initialize stats
if COMPUTE_STATS
    % Cost function
    Bfit = @(B) 0.5 * norm(vec(P - D * B))^2;
    if strcmpi(type,'hard')
        % Ell-0 regularization
        Breg = @(B) 0.5 * mu^2 * nnz(B);
    elseif strcmpi(type,'soft')
        % Ell-0 regularization
        Breg = @(B) mu * norm(vec(B),1);
    else
        % Unsupported shrinkage
        error('Unsupported shrinkage type ''%s''',type);
    end
    Psi = @(B) Bfit(B) + Breg(B);
    
    % Stats-printing function
    iterFmt = sprintf('%%0%dd',ceil(log10(nIters + 1)));
    out     = printFcn('Iter' ,iterFmt,'cost','%.2f','nrmse','%.3f', ...
                       'delta','%.3e' ,'time','%.2fs');
    
    % Initialize stats
    cost  = nan(1,nIters);
    nrmse = nan(1,nIters);
    delta = nan(1,nIters);
    time  = nan(1,nIters);
end

% Sparse Coding
if PRINT_STATS
    fprintf('***** Sparse Coding *****\n');
end
if accel
    % Initialize accelerated method
    t     = 0;
    Blast = B;
end
for it = 1:nIters
    % Initilize iteration
    itimer = tic;
    
    % Proximal gradient update
    if accel
        % Accelerated step
        tlast = t;
        t     = 0.5 * (1 + sqrt(1 + 4 * t^2));
        Bbar  = B + ((tlast - 1) / t) * (B - Blast);
        Blast = B;
        Gbar  = D' * (D * Bbar - P);
        B     = shrink(Bbar - tau * Gbar,tau);
    else
        % Standard step
        Blast = B;
        G     = D' * (D * B - P);
        B     = shrink(B - tau * G,tau);
    end
    
    % Record stats
    if COMPUTE_STATS
        cost(it)  = Psi(B);
        nrmse(it) = NRMSEfcn(B,Btrue);
        delta(it) = computeNRMSE(B,Blast);
        time(it)  = toc(itimer);
        if PRINT_STATS
            out(it,cost(it),nrmse(it),delta(it),time(it)); 
        end
    end
end

% Return stats
if COMPUTE_STATS
    stats.nIters = nIters;
    stats.cost   = cost;
    stats.nrmse  = nrmse;
    stats.delta  = delta;
    stats.time   = time;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Local functions
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Hard thresholding
function Y = hard(X,lambda)
Y = X .* (abs(X) > lambda);

% Soft thresholding
function Y = soft(X,lambda)
Y = sign(X) .* max(abs(X) - lambda,0);

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
function [nIters, type, B0, Btrue, NRMSEfcn, accel, tau, flag] = ...
                                                      parseInputs(P,D,opts)
if ~exist('opts','var') || isempty(opts)
    opts = struct();
end

% Standard args
nIters   = parseField(opts,'nIters',10);
type     = parseField(opts,'type','hard');
B0       = parseField(opts,'B0',nan);
Btrue    = parseField(opts,'Btrue',nan);
NRMSEfcn = parseField(opts,'NRMSEfcn',@computeNRMSE);
accel    = parseField(opts,'accel',true);
tau      = parseField(opts,'tau',nan);
flag     = parseField(opts,'flag',1);

% Expensive defaults
if isnan(B0),  B0  = D' * P;                        end
if isnan(tau), tau = (0.99 + ~accel) / norm(D)^2;   end

% Parse struct field
function val = parseField(stats,field,default)
if isfield(stats,field)
    val = stats.(field);
else
    val = default;
end
