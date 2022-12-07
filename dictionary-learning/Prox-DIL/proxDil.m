function [D, B, stats] = proxDil(Y,mu,varargin)
%
% Syntax:       [D, B] = proxDil(Y,mu);
%               [D, B] = proxDil(Y,mu,opts);
%               [D, B, stats] = proxDil(Y,mu);
%               [D, B, stats] = proxDil(Y,mu,opts);
%               
% Inputs:       Y is a d x n data matrix
%               
%               mu >= 0 is the sparsity regularization parameter
%               
%               [OPTIONAL] opts is a struct containing one or more of the
%               following fields. The default values are in ()
%                   
%                   opts.D0 (dctmtx(d)') is a d x m matrix containing the
%                   initial dictionary. Overcomplete (m > d) and
%                   undercomplate (m < d) dictionaries are allowed
%                   
%                   opts.B0 (zeros(m,n)) is an m x n matrix containing the
%                   initial sparse codes
%                   
%                   opts.type ('hard') can be 'hard' or 'soft' and
%                   specifies whether apply hard or soft thresholding to B
%                   
%                   opts.dr (nan) is a rank constraint on the dictionary
%                   atoms. Note: both dr and ddim must be non-nan to apply
%                   a rank constraint
%                   
%                   opts.ddim (nan) is a 1 x 2 vector describing how to
%                   reshape the dictionary atoms into a matrix before 
%                   applying the rank constraint. Note: both dr and ddim
%                   must be non-nan to apply a rank constraint
%                   
%                   opts.fixedD (false) determines whether to fix the
%                   initial dictionary (true) or learn it (false)
%                   
%                   opts.nIters (50) is the number of proximal gradient
%                   iterations to perform
%                   
%                   opts.accel (true) specifies whether to use accelerated
%                   proximal gradient steps
%                   
%                   opts.tau (0.99 + ~accel) is the step size parameter to
%                   use, and should satisfy
%                   
%                       tau <= 1, when accel = true
%                       tau <  2, when accel = false
%                   
%                   The actual step sizes used are
%                       
%                       D updates: tau / norm(B)^2
%                       B updates: tau / norm(D)^2
%                   
%                   Alternatively, tau can be a 1 x 2 or nIters x 2 matrix
%                   of the form tau = [tauB, tauD] specifying some fixed
%                   step sizes to use (will not be normalized internally)
%                   
%                   opts.flag (1) determines what status updates to print
%                   to the command window. The choices are
%                   
%                       flag = 0: no printing
%                       flag = 1: print iteration status
%               
% Outputs:      D is a d x m matrix containing the final dictionary
%               
%               B is an m x n matrix containing the final sparse codes
%               
%               stats is a struct containing the following fields:
%               
%                   stats.nIters is the number of iterations performed
%                   
%                   stats.cost is a 1 x nIters vector containing the cost
%                   function values at each iteration
%                   
%                   stats.repError is a 1 x nIters vector containing the
%                   representation error at each iteration, defined as
%                   \|Y - D_k B_k\|_F / \|Y\|_F
%                   
%                   stats.deltaD is a 1 x nIters vector containing the
%                   relative convergence of D at each iteration, defined as
%                   \|D_k - D_{k - 1}\|_F / \|D_{k - 1}\|_F
%                   
%                   stats.deltaB is a 1 x nIters vector containing the
%                   relative convergence of B at each iteration, defined as
%                   \|B_k - B_{k - 1}\|_F / \|B_{k - 1}\|_F
%                   
%                   stats.sparsity is a 1 x nIters vector containing the
%                   sparsity, in percent, of B at each iteration
%                   
%                   stats.time is a 1 x nIters vector containing the time,
%                   in seconds, required to perform each outer iteration
%               
% Description:  Performs proximal gradient dictionary learning (Prox-DIL)
%               by solving one of the following problems:
%               
%               When type == 'hard':
%               
%               min_{D,B} 0.5 \|Y - D B\|_F^2 + 0.5 \mu^2 \|B\|_0
%               
%               Or, when type == 'soft':
%               
%               min_{D,B} 0.5 \|Y - D B\|_F^2 + \mu \|B\|_1
%               
% Author:       Brian Moore
%               brimoor@umich.edu
%               
% Date:         January 31, 2017
%               July 20, 2017
%

% Parse inputs
[D, B, type, dr, ddim, fixedD, nIters, accel, tau, flag] = ...
                                                parseInputs(Y,varargin{:});
PRINT_STATS   = (flag > 0);
COMPUTE_STATS = PRINT_STATS || (nargout == 3);
OPTIMIZE_DICT = ~fixedD;
FIXED_STEPS   = ~isscalar(tau);

% Parse tau
if FIXED_STEPS && (size(tau,1) == 1)
    tau = repmat(tau,nIters,1);
end

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
    Yfit = @(D,B) 0.5 * norm(vec(Y - D * B))^2;
    if strcmpi(type,'hard')
        % Ell-0 regularization
        Bfit = @(B) 0.5 * mu^2 * nnz(B);
    elseif strcmpi(type,'soft')
        % Ell-1 regularization
        Bfit = @(B) mu * norm(vec(B),1);
    else
        % Unsupported shrinkage
        error('Unsupported shrinkage type ''%s''',type);
    end
    Psi = @(D,B) Yfit(D,B) + Bfit(B);
    
    % Stats-printing function
    iterFmt = sprintf('%%0%dd',ceil(log10(nIters + 1)));
    out     = printFcn('Iter'    ,iterFmt,'cost'    ,'%.2f'  , ...
                       'repError','%.3f' ,'deltaD'  ,'%.3e'  , ...
                       'deltaB'  ,'%.3e' ,'sparsity','%.2f%%', ...
                       'time'    ,'%.2fs');
    
    % Initialize stats
    cost     = nan(1,nIters);
    repError = nan(1,nIters);
    deltaD   = nan(1,nIters);
    deltaB   = nan(1,nIters);
    sparsity = nan(1,nIters);
    time     = nan(1,nIters);
end

% Prox-DIL
if PRINT_STATS
    fprintf('***** Prox-DIL *****\n');
end
d  = size(D,1);
d0 = [1; zeros(d - 1,1)]; % Zero-reset atom
if accel
    % Initialize accelerated method
    t     = 0;
    Blast = B;
    Dlast = D;
end
for it = 1:nIters
    % Initialize iteration
    itimer = tic;
    if accel
        % Nesterov parameter
        tlast = t;
        t     = 0.5 * (1 + sqrt(1 + 4 * t^2));
    end
    
    % B step size
    if FIXED_STEPS
        % User-specified step size
        tauB = tau(it,1);
    else
        % Scaled step size
        tauB = tau / norm(D)^2;
    end
    
    % B update
    if accel
        % Accelerated proximal gradient
        Bbar  = B + ((tlast - 1) / t) * (B - Blast);
        Blast = B;
        GBbar = D' * (D * Bbar - Y);
        B     = shrink(Bbar - tauB * GBbar,tauB);
    else
        % Proximal gradient
        Blast = B;
        GB    = D' * (D * B - Y);
        B     = shrink(B - tauB * GB,tauB);
    end
    
    if OPTIMIZE_DICT
        % D step size
        if FIXED_STEPS
            % User-specified step size
            tauD = tau(it,2);
        else
            % Scaled step size
            tauD = tau / norm(B)^2;
        end
        
        % D update
        if accel
            % Accelerated proximal gradient
            Dbar   = D + ((tlast - 1) / t) * (D - Dlast);
            Dlast  = D;
            GDbar  = (Dbar * B - Y) * B';
            D = projectAtoms(Dbar - tauD * GDbar,dr,ddim,d0);
        else
            % Proximal gradient
            Dlast = D;
            GD    = (D * B - Y) * B';
            D     = projectAtoms(D - tauD * GD,dr,ddim,d0);
        end
    end
    
    % Record stats
    if COMPUTE_STATS
        cost(it)     = Psi(D,B);
        repError(it) = computeNRMSE(D * B,Y);
        deltaD(it)   = computeNRMSE(D,Dlast);
        deltaB(it)   = computeNRMSE(B,Blast);
        sparsity(it) = 100 * (nnz(B) / numel(B));
        time(it)     = toc(itimer);
        if PRINT_STATS
            out(it,cost(it),repError(it),deltaD(it),deltaB(it), ...
                sparsity(it),time(it));
        end
    end
end

%{
figure();
subplot(1,2,1); plot(1:nIters,LB,'b-o'); title('LB'); axis tight; padAxis();
subplot(1,2,2); plot(1:nIters,LD,'r-o'); title('LD'); axis tight; padAxis();
%}

% Return stats
if COMPUTE_STATS
    stats.nIters   = nIters;
    stats.cost     = cost;
    stats.repError = repError;
    stats.deltaD   = deltaD;
    stats.deltaB   = deltaB;
    stats.sparsity = sparsity;
    stats.time     = time;
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
function [D0, B0, type, dr, ddim, fixedD, nIters, accel, tau, flag] = ...
                                                        parseInputs(Y,opts)
if ~exist('opts','var') || isempty(opts)
    opts = struct();
end

% Standard args
D0     = parseField(opts,'D0',nan);
B0     = parseField(opts,'B0',nan);
type   = parseField(opts,'type','hard');
dr     = parseField(opts,'dr',nan);
ddim   = parseField(opts,'ddim',nan);
fixedD = parseField(opts,'fixedD',false);
nIters = parseField(opts,'nIters',50);
accel  = parseField(opts,'accel',true);
tau    = parseField(opts,'tau',nan);
flag   = parseField(opts,'flag',1);

% Expensive defaults
if isnan(D0),   D0 = dctmtx(size(Y,1))';            end
if isnan(B0),   B0 = zeros(size(D0,2),size(Y,2));   end
if isnan(tau), tau = (0.99 + ~accel);               end

% Parse struct field
function val = parseField(stats,field,default)
if isfield(stats,field)
    val = stats.(field);
else
    val = default;
end
