function [D, stats] = dictUpdate(P,B,varargin)
%
% Syntax:       D = dictUpdate(P,B);
%               D = dictUpdate(P,B,opts);
%               [D, stats] = dictUpdate(P,B);
%               [D, stats] = dictUpdate(P,B,opts);
%               
% Inputs:       P is an d x n data matrix
%               
%               B is a m x n matrix of sparse codes
%               
%               [OPTIONAL] opts is a struct containing one or more of the
%               following fields. The default values are in ()
%                   
%                   opts.nIters (10) is the number of iterations to
%                   perform
%                   
%                   opts.dr (nan) is a rank constraint on the dictionary
%                   atoms. Note: both dr and ddim must be non-nan to apply
%                   a rank constraint
%                   
%                   opts.ddim ([prod(pdim(1:(end - 1))), pdim(end)]) is a
%                   1 x 2 vector describing how to reshape the dictionary
%                   atoms into a matrix before applying the rank constraint
%                   Note: both dr and ddim must be non-nan to apply a rank
%                   constraint
%                   
%                   opts.D0 (P * B') is a m x n matrix specifying the
%                   initial D iterate
%                   
%                   opts.Dtrue (nan) is the ground truth D matrix to use 
%                   for MSE calculations
%                   
%                   opts.NRMSEfcn (@computeNRMSE) is the function to use
%                   when computing the NRMSE of D after each iteration
%                   
%                   opts.accel (true) specifies whether to use accelerated
%                   proximal gradient steps
%                   
%                   opts.tau ((0.99 + ~accel) / norm(B)^2) is the step size
%                   parameter, and should satisfy
%                   
%                       tau <= 1 / norm(B)^2, when accel = true
%                       tau <  2 / norm(B)^2, when accel = false
%                   
%                   opts.flag (1) determines what status updates to print
%                   to the command window. The choices are
%                   
%                       flag = 0: no printing
%                       flag = 1: print iteration status
%               
% Outputs:      D is a d x m matrix containing the estimated dictionary
%               
%               stats is a struct containing the following fields:
%               
%                   stats.nIters is the number of iterations performed
%                   
%                   stats.cost is a 1 x nIters vector containing the cost
%                   function values at each iteration
%                   
%                   stats.nrmse is a 1 x nIters vector containing the NRMSE
%                   of D with respect to Dtrue at each iteration
%                   
%                   stats.delta is a 1 x nIters vector containing the
%                   relative convergence of D at each iteration, defined as
%                   \|D_{k + 1} - D_k\|_F / \|D_k\|_F
%                   
%                   stats.time is a 1 x nIters vector containing the time,
%                   in seconds, required to perform each outer iteration
%               
% Description:  Solves the dictionary update problem
%               
%                 \min_{D}  \|P - DB\|_F^2
%               
%               subject to  rank(R(D(:,k)) <= dr
%                             \|D(:,k)\|_2 == 1
%               
% Author:       Brian Moore
%               brimoor@umich.edu
%               
% Date:         January 31, 2017
%

% Parse inputs
[nIters, dr, ddim, D, Dtrue, NRMSEfcn, accel, tau, flag] = ...
                                              parseInputs(P,B,varargin{:});
d     = size(P,1);
mOrig = size(B,1);
idx   = any(B,2);
PRINT_STATS   = (flag > 0);
COMPUTE_STATS = PRINT_STATS || (nargout == 2);
UNUSED_ATOMS  = ~all(idx);

% Omit unused atoms, if necessary
if UNUSED_ATOMS
    if isequal(size(Dtrue),size(D))
        Dtrue = Dtrue(:,idx);
    end
    D = D(:,idx);
    B = B(idx,:);
end

% Initialize stats
if COMPUTE_STATS
    % Cost function
    Psi = @(D) 0.5 * norm(vec(P - D * B))^2;
    
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

% Dictionary Update
if PRINT_STATS
    fprintf('***** Dictionary Update *****\n');
end
d0 = [1; zeros(d - 1,1)]; % Zero-reset atom
if nIters <= 0
    % Single update
    itimer = tic;
    
    % Least-squares
    %D = (P * B') / (1e-6 * eye(size(B,1)) + B * B');
    D = P * pinv(B);
    
    % Project atoms
    D = projectAtoms(D,dr,ddim,d0);
    
    % Record stats
    if COMPUTE_STATS
        cost  = Psi(D);
        nrmse = NRMSEfcn(D,Dtrue);
        delta = 0;
        time  = toc(itimer);
        if PRINT_STATS
            out(0,cost,nrmse,delta,time); 
        end
    end
else
    % Proximal gradient
    if accel
        % Initialize accelerated method
        t     = 0;
        Dlast = D;
    end
    for it = 1:nIters
        % Initilize iteration
        itimer = tic;

        % Proximal gradient update
        if accel
            % Accelerated step
            tlast = t;
            t     = 0.5 * (1 + sqrt(1 + 4 * t^2));
            Dbar  = D + ((tlast - 1) / t) * (D - Dlast);
            Dlast = D;
            Gbar  = (Dbar * B - P) * B';
            D     = projectAtoms(Dbar - tau * Gbar,dr,ddim,d0);
        else
            % Standard step
            Dlast = D;
            G     = (D * B - P) * B';
            D     = projectAtoms(D - tau * G,dr,ddim,d0);
        end

        % Record stats
        if COMPUTE_STATS
            cost(it)  = Psi(D);
            nrmse(it) = NRMSEfcn(D,Dtrue);
            delta(it) = computeNRMSE(D,Dlast);
            time(it)  = toc(itimer);
            if PRINT_STATS
                out(it,cost(it),nrmse(it),delta(it),time(it)); 
            end
        end
    end
end

% Reconstitute full dictionary, if necessary
if UNUSED_ATOMS
    tmp = D;
    D   = repmat(d0,1,mOrig);
    D(:,idx) = tmp;
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

% Project atoms
% Note: assumes that d0 is unit-norm with reshaped rank <= dr
function Dp = projectAtoms(D,dr,ddim,d0)
% Parse inputs
LOW_RANK_ATOMS = (dr < min(ddim));
ZERO_TOL = 1e-8;

% Low-rank atoms
if LOW_RANK_ATOMS
    m = size(D,2);
    for k = 1:m
        [Uk, Sk, Vk] = svd(reshape(D(:,k),ddim),'econ');
        D(:,k) = vec(Uk(:,1:dr) * Sk(1:dr,1:dr) * Vk(:,1:dr)');
    end
end

% Normalize atoms
normD = sqrt(sum(abs(D).^2,1));
idx   = (normD > ZERO_TOL);
nZero = nnz(~idx);
if nZero > 0
    Dp = repmat(d0,1,numel(idx));
    Dp(:,idx) = bsxfun(@rdivide,D(:,idx),normD(idx));
else
    Dp = bsxfun(@rdivide,D,normD);
end

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
function [nIters, dr, ddim, D0, Dtrue, NRMSEfcn, accel, tau, flag] = ...
                                                      parseInputs(P,B,opts)
if ~exist('opts','var') || isempty(opts)
    opts = struct();
end

% Standard args
nIters   = parseField(opts,'nIters',10);
dr       = parseField(opts,'dr',nan);
ddim     = parseField(opts,'ddim',nan);
D0       = parseField(opts,'D0',nan);
Dtrue    = parseField(opts,'Dtrue',nan);
NRMSEfcn = parseField(opts,'NRMSEfcn',@computeNRMSE);
accel    = parseField(opts,'accel',true);
tau      = parseField(opts,'tau',nan);
flag     = parseField(opts,'flag',1);

% Expensive defaults
if isnan(D0),  D0  = P * B';                        end
if isnan(tau), tau = (0.99 + ~accel) / norm(B)^2;   end

% Parse struct field
function val = parseField(stats,field,default)
if isfield(stats,field)
    val = stats.(field);
else
    val = default;
end
