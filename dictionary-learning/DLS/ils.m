function [x, stats] = ils(A,b,varargin)
%
% Syntax:       x = ils(A,b);
%               x = ils(A,b,opts);
%               [x, stats] = ils(A,b);
%               [x, stats] = ils(A,b,opts);
%               
% Inputs:       A is an m x n system matrix
%               
%               b is an m x 1 data vector
%               
%               [OPTIONAL] opts is a struct containing one or more of the
%               following fields. The default values are in ()
%                   
%                   opts.nIters (50) is the number of iterations to perform
%                   
%                   opts.x0 (A' * b) is a n x 1 vector specifying the
%                   initial x iterate
%                   
%                   opts.xtrue (nan) is the ground truth x vector to use 
%                   for NRMSE calculations
%                   
%                   opts.NRMSEfcn (@computeNRMSE) is the function to use
%                   when computing the NRMSE of x after each iteration
%                   
%                   opts.accel (true) specifies whether to use accelerated
%                   gradient descent steps
%                   
%                   opts.tau ((0.99 + ~accel) / norm(A)^2) is the step size
%                   parameter, and should satisfy
%                   
%                       tau <= 1 / norm(A)^2, when accel = true
%                       tau <  2 / norm(A)^2, when accel = false
%                   
%                   opts.flag (1) determines what status updates to print
%                   to the command window. The choices are
%                   
%                       flag = 0: no printing
%                       flag = 1: print iteration status
%               
% Outputs:      x is the estimated n x 1 vector
%               
%               stats is a struct containing the following fields:
%               
%                   stats.nIters is the number of iterations performed
%                   
%                   stats.cost is a 1 x nIters vector containing the cost
%                   function values at each iteration
%                   
%                   stats.nrmse is a 1 x nIters vector containing the NRMSE
%                   of x with respect to xtrue at each iteration
%                   
%                   stats.delta is a 1 x nIters vector containing the
%                   relative convergence of x at each iteration, defined as
%                   \|x_{k + 1} - x_k\|_F / \|x_k\|_F
%                   
%                   stats.time is a 1 x nIters vector containing the time,
%                   in seconds, required to perform each outer iteration
%               
% Description:  Solves the least squares problem
%               
%               min_{x} 0.5\|b - Ax\|_2^2
%               
%               via gradient descent
%               
% Author:       Brian Moore
%               brimoor@umich.edu
%               
% Date:         May 10, 2016
%

% Parse inputs
[nIters, x, xtrue, NRMSEfcn, accel, tau, flag] = ...
                                              parseInputs(A,b,varargin{:});
PRINT_STATS   = (flag > 0);
COMPUTE_STATS = PRINT_STATS || (nargout == 2);

% Initialize stats
if COMPUTE_STATS
    % Cost function
    Psi = @(x) 0.5 * norm(b - A * x)^2;
    
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

% Iterative Least Squares
if PRINT_STATS
    if accel
        fprintf('***** ILS *****\n');
    else
        fprintf('***** NILS *****\n');
    end
end
if accel
    % Initialize accelerated method
    t     = 0;
    xlast = x;
end
for it = 1:nIters
    % Initilize iteration
    itimer = tic;
    
    % Gradient descent update
    if accel
        % Accelerated step
        tlast = t;
        t     = 0.5 * (1 + sqrt(1 + 4 * t^2));
        xbar  = x + ((tlast - 1) / t) * (x - xlast);
        xlast = x;
        gbar  = A' * (A * xbar - b);
        x     = xbar - tau * gbar;
    else
        % Standard step
        xlast = x;
        g     = A' * (A * x - b);
        x     = x - tau * g;
    end
    
    % Record stats
    if COMPUTE_STATS
        cost(it)  = Psi(x);
        nrmse(it) = NRMSEfcn(x,xtrue);
        delta(it) = computeNRMSE(x,xlast);
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

% Parse inputs
function [nIters, x0, xtrue, NRMSEfcn, accel, tau, flag] = ...
                                                      parseInputs(A,b,opts)
if ~exist('opts','var') || isempty(opts)
    opts = struct();
end

% Standard args
nIters   = parseField(opts,'nIters',50);
x0       = parseField(opts,'x0',nan);
xtrue    = parseField(opts,'xtrue',nan);
NRMSEfcn = parseField(opts,'NRMSEfcn',@computeNRMSE);
accel    = parseField(opts,'accel',true);
tau      = parseField(opts,'tau',nan);
flag     = parseField(opts,'flag',1);

% Expensive defaults
if isnan(x0),  x0   = A' * b;                      end
if isnan(tau), tau  = (0.99 + ~accel) / norm(A)^2; end

% Parse struct field
function val = parseField(stats,field,default)
if isfield(stats,field)
    val = stats.(field);
else
    val = default;
end
