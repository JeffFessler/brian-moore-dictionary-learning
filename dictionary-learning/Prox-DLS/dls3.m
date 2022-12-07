function [X, D, B, stats] = dls3(Y,lambda,mu,varargin)
%
% Syntax:       [X, D, B] = dls3(Y,lambda,mu);
%               [X, D, B] = dls3(Y,lambda,mu,opts);
%               [X, D, B, stats] = dls3(Y,lambda,mu);
%               [X, D, B, stats] = dls3(Y,lambda,mu,opts);
%               
% Inputs:       Y is an m x n data matrix
%               
%               lambda >= 0 is the dictionary regularization parameter
%               
%               mu >= 0 is the sparsity regularization parameter.
%               Alternatively, mu can be a 1 x nIters vector specifying a
%               different mu value for each outer iteration
%               
%               [OPTIONAL] opts is a struct containing one or more of the
%               following fields. The default values are in ()
%                   
%                   opts.A (1) is an m x p system matrix
%                   
%                   opts.M (1) is an m x n weight matrix
%                   
%                   opts.xdim (size(X0)) are the dimensions of X to use
%                   when extracting patches. Note that xdim need not be the
%                   same as size(X) as long as the two representations are
%                   are interchangable via reshape(). For example, you
%                   could specify xdim = [p, n1, n2] where n1 * n2 = n
%                   
%                   opts.pdim ([8, 8]) is a vector with numel(xdim)
%                   elements specifying the patch sizes to extract from X
%                   
%                   opts.pgap ([2, 2]) is a vector with numel(xdim) 
%                   elements specfiying the patch strides (i.e., shifts) to
%                   use along each dimension of X
%                   
%                   opts.type ('hard') can be {'hard','soft'} and specifies
%                   whether to use hard or soft thresholding to regularize
%                   the sparse codes B
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
%                   opts.fixedD (false) determines whether to fix the
%                   initial dictionary (true) or learn it (false)
%                   
%                   opts.nIters (50) is the number of outer iterations to
%                   perform
%                   
%                   opts.nItersB (5) is the number of inner B updates
%                   to perform per outer iteration
%                   
%                   opts.nItersD (5) is the number of inner D updates
%                   to perform per outer iteration
%                   
%                   opts.nItersX (5) is the number of inner X updates to
%                   perform per outer iteration. Note that, when A = 1, the
%                   exact X update is always computed and nItersX is
%                   ignored
%                   
%                   opts.X0 (A' * Y) is a p x n matrix specifying the
%                   initial X iterate
%                   
%                   opts.D0 (dctmtx(b)') is a b x c matrix specifying the
%                   initial dictionary. Overcomplete dictionaries (c > b)
%                   are allowed
%                   
%                   opts.B0 (zeros(c,d)) is an c x d matrix specifying 
%                   the initial sparse codes
%                   
%                   opts.Xtrue (nan) is the ground truth X matrix to use 
%                   for NRMSE calculations
%                   
%                   opts.NRMSEfcn (@computeNRMSE) is the function to use
%                   when computing the NRMSE of X after each iteration
%                   
%                   opts.accelX (true) specifies whether to use Nesterov's
%                   acceleration scheme for the X updates
%                   
%                   opts.accelD (true) specifies whether to use Nesterov's
%                   acceleration scheme for the D updates
%                   
%                   opts.accelB (true) specifies whether to use Nesterov's
%                   acceleration scheme for the B updates
%                   
%                   opts.tauX ((0.99 + ~accelX) / norm(A)^2) is the step
%                   size parameter for the X updates, and should satisfy
%                   
%                       tauX <= 1 / norm(A)^2, when accelX = true
%                       tauX <  2 / norm(A)^2, when accelX = false
%                   
%                   opts.tauD (0.99 + ~accelD) is the step size parameter
%                   for the D updates, and should satisfy
%                   
%                       tauB <= 1, when accelD = true
%                       tauB <  2, when accelD = false
%                   
%                   The actual step size used is tauD / norm(B)^2
%                   
%                   opts.tauB (0.99 + ~accelB) is the step size parameter
%                   for the B updates, and should satisfy
%                   
%                       tauB <= 1, when accelB = true
%                       tauB <  2, when accelB = false
%                   
%                   The actual step size used is tauB / norm(D)^2
%                   
%                   opts.flag (1) determines what status updates to print
%                   to the command window. The choices are
%                   
%                       flag = 0: no printing
%                       flag = 1: print outer iteration status
%                       flag = 2: print inner and outer iteration status
%               
% Outputs:      X is the estimated p x n matrix
%               
%               D is a b x c matrix, where b = prod(pdim), c = size(D0,2), 
%               whose columns contain the vectorized dictionary elements
%               
%               B is a c x d matrix, where d = total # patches, whose
%               columns contain the sparse codes for the patches of X
%               
%               stats is a struct containing the following fields:
%               
%                   stats.nIters is the number of outer iterations
%                   performed
%                   
%                   stats.cost is a 1 x nIters vector containing the cost
%                   function values at each iteration
%                   
%                   stats.nrmse is a 1 x nIters vector containing the NRMSE
%                   of X with respect to Xtrue at each iteration
%                   
%                   stats.deltaX is a 1 x nIters vector containing the
%                   relative convergence of X at each iteration, defined as
%                   \|X_{k + 1} - X_k\|_F / \|X_k\|_F
%                   
%                   stats.deltaD is a 1 x nIters vector containing the
%                   relative convergence of D at each iteration, defined as
%                   \|D_{k + 1} - D_k\|_F / \|D_k\|_F
%                   
%                   stats.deltaB is a 1 x nIters vector containing the
%                   relative convergence of B at each iteration, defined as
%                   \|B_{k + 1} - B_k\|_F / \|B_k\|_F
%                   
%                   stats.sparsity is a 1 x nIters vector containing the
%                   sparsity, in percent, of B at each iteration
%                   
%                   stats.time is a 1 x nIters vector containing the time,
%                   in seconds, required to perform each outer iteration
%               
% Description:  Solves one of the following Dictionary Least-Sqaures (DLS) 
%               problems:
%               
%               When type == 'hard':
%               
%               \min_{X,D,B} 0.5 \|\sqrt{M} \odot (Y - A X)\|_F^2 + 
%                   \lambda (0.5 \|P(X) - D B\|_F^2 + 0.5 \mu^2 \|B\|_0)
%               
%               When type == 'soft':
%               
%               \min_{X,D,B} 0.5 \|\sqrt{M} \odot (Y - A X)\|_F^2 + 
%                   \lambda (0.5 \|P(X) - D B\|_F^2 + \mu \|B\|_1)
%               
% Author:       Brian Moore
%               brimoor@umich.edu
%               
% Date:         January 31, 2017
%               June 6, 2017
%

% Parse inputs
[A, M, xdim, pdim, pgap, type, dr, ddim, fixedD, nIters, nItersB, ...
 nItersD, nItersX, X, D, B, Xtrue, NRMSEfcn, accelX, accelD, accelB, ...
 tauX, tauD, tauB, flag] = parseInputs(Y,varargin{:});
if isscalar(mu)
    mu = repmat(mu,1,nIters);
end
PRINT_STATS   = (flag > 0);
COMPUTE_STATS = PRINT_STATS || (nargout == 4);
OPTIMIZE_DICT = ~fixedD;

% Initialize stats
if COMPUTE_STATS
    % Cost function
    if strcmpi(type,'hard')
        % Ell-0 regularization
        Psi = @(X,P,DB,B) 0.5 * norm(vec(sqrt(M) .* (Y - A * X)))^2 + ...
                 lambda * 0.5 * (norm(vec(P - DB))^2 + mu(end)^2 * nnz(B));
    elseif strcmpi(type,'soft')
        % Ell-1 regularization
        Psi = @(X,P,DB,B) 0.5 * norm(vec(sqrt(M) .* (Y - A * X)))^2 + ...
           lambda * (0.5 * norm(vec(P - DB))^2 + mu(end) * norm(vec(B),1));
    else
        % Unsupported shrinkage
        error('Unsupported shrinkage type ''%s''',type);
    end
    
    % Stats-printing function
    iterFmt = sprintf('%%0%dd',ceil(log10(nIters + 1)));
    out     = printFcn('Iter'    ,iterFmt ,'cost'  ,'%.2f', ...
                       'nrmse'   ,'%.3f'  ,'deltaX','%.3e', ...
                       'deltaD'  ,'%.3e'  ,'deltaB','%.3e', ...
                       'sparsity','%.2f%%','time'  ,'%.2fs');
    
    % Initialize stats
    cost     = nan(1,nIters);
    nrmse    = nan(1,nIters);
    deltaX   = nan(1,nIters);
    deltaD   = nan(1,nIters);
    deltaB   = nan(1,nIters);
    sparsity = nan(1,nIters);
    time     = nan(1,nIters);
end

% Initialize patches
[Mp, ~, Np] = patchInds(xdim,pdim,pgap); % Patch indices
Np = reshape(Np,size(X));                % Patch counts
P  = X(Mp);                              % Initial patches

% Initialize sparseCoding opts
optsB.nIters = nItersB;
optsB.type   = type;
optsB.accel  = accelB;
optsB.flag   = flag - 1;

% Initialize updateDict opts
optsD.nIters = nItersD;
optsD.dr     = dr;
optsD.ddim   = ddim;
optsD.accel  = accelD;
optsD.flag   = flag - 1;

% Initialize patchLS opts
optsX.A        = A;
optsX.M        = M;
optsX.nIters   = nItersX;
optsX.Xtrue    = Xtrue;
optsX.NRMSEfcn = NRMSEfcn;
optsX.accel    = accelX;
optsX.tau      = tauX;
optsX.flag     = flag - 1;

% Dictionary Least-Squares
if PRINT_STATS
    fprintf('***** Dictionary Least-Squares *****\n');
end
if isnan(B)
    % Initialize sparse codes
    B = zeros(size(D,2),size(P,2));       
end
for it = 1:nIters
    % Initialize iteration
    itimer = tic;
    if COMPUTE_STATS
        % Save last iterates
        Xlast = X;
        Dlast = D;
        Blast = B;
    end
    
    % B update(s)
    optsB.B0 = B;
    optsB.tau = tauB / norm(D)^2;
    B = sparseCoding(P,D,mu(it),optsB);
    
    % D update
    if OPTIMIZE_DICT
        optsD.D0  = D;
        optsD.tau = tauD / norm(B)^2;
        D = dictUpdate(P,B,optsD);
    end
    DB = D * B;
    
    % X update(s)
    optsX.X0 = X;
    X = patchLS(Y,DB,lambda,Mp,Np,optsX);
    P = X(Mp);
    
    % Record stats
    if COMPUTE_STATS
        cost(it)     = Psi(X,P,DB,B);
        nrmse(it)    = NRMSEfcn(X,Xtrue);
        deltaX(it)   = computeNRMSE(X,Xlast);
        deltaD(it)   = computeNRMSE(D,Dlast);
        deltaB(it)   = computeNRMSE(B,Blast);
        sparsity(it) = 100 * (nnz(B) / numel(B));
        time(it)     = toc(itimer);
        if PRINT_STATS
            out(it,cost(it),nrmse(it),deltaX(it),deltaD(it),deltaB(it), ...
                sparsity(it),time(it)); 
        end
    end
end

% Return stats
if COMPUTE_STATS
    stats.nIters   = nIters;
    stats.cost     = cost;
    stats.nrmse    = nrmse;
    stats.deltaX   = deltaX;
    stats.deltaD   = deltaD;
    stats.deltaB   = deltaB;
    stats.sparsity = sparsity;
    stats.time     = time;
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

% Vectorize input
function x = vec(X)
x = X(:);

% Parse inputs                                                
function [A, M, xdim, pdim, pgap, type, dr, ddim, fixedD, nIters, ...
          nItersB, nItersD, nItersX, X0, D0, B0, Xtrue, NRMSEfcn, ...
          accelX, accelD, accelB, tauX, tauD, tauB, flag] = ...
                                                        parseInputs(Y,opts)
if ~exist('opts','var') || isempty(opts)
    opts = struct();
end

% Standard args
A        = parseField(opts,'A',1);
M        = parseField(opts,'M',1);
%xdim dependent
pdim     = parseField(opts,'pdim',[8, 8]);
pgap     = parseField(opts,'pgap',[2, 2]);
type     = parseField(opts,'type','hard');
dr       = parseField(opts,'dr',nan);
%ddim dependent
fixedD   = parseField(opts,'fixedD',false);
nIters   = parseField(opts,'nIters',50);
nItersB  = parseField(opts,'nItersB',5);
nItersD  = parseField(opts,'nItersD',5);
nItersX  = parseField(opts,'nItersX',5);
X0       = parseField(opts,'X0',nan);
D0       = parseField(opts,'D0',nan);
B0       = parseField(opts,'B0',nan);
Xtrue    = parseField(opts,'Xtrue',nan);
NRMSEfcn = parseField(opts,'NRMSEfcn',@computeNRMSE);
accelX   = parseField(opts,'accelX',true);
accelD   = parseField(opts,'accelD',true);
accelB   = parseField(opts,'accelB',true);
tauX     = parseField(opts,'tauX',nan);
tauD     = parseField(opts,'tauD',nan);
tauB     = parseField(opts,'tauB',nan);
flag     = parseField(opts,'flag',1);

% Expensive defaults
if isnan(X0),     X0  = A' * Y;                         end
if isnan(D0),     D0  = dctmtx(prod(pdim))';            end
if isnan(tauX), tauX = (0.99 + ~accelX) / norm(A)^2;    end
if isnan(tauD), tauD = (0.99 + ~accelD);                end
if isnan(tauB), tauB = (0.99 + ~accelB);                end

% Dependent args
xdim     = parseField(opts,'xdim',size(X0));
ddim     = parseField(opts,'ddim',[prod(pdim(1:(end - 1))), pdim(end)]);

% Parse struct field
function val = parseField(stats,field,default)
if isfield(stats,field)
    val = stats.(field);
else
    val = default;
end
