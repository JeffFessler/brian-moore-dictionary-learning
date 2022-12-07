function zzz_OnAIR(Y,lambda,mu,gamma,varargin)
%
% Syntax:       [Xt, D, Bt] = OnAIR(Y,lambda,mu,gamma);
%               [Xt, D, Bt] = OnAIR(Y,lambda,mu,gamma,opts);
%               [Xt, D, Bt, stats] = OnAIR(Y,lambda,mu,gamma);
%               [Xt, D, Bt, stats] = OnAIR(Y,lambda,mu,gamma,opts);
%               
% Inputs:       Yt is an m x n data matrix
%               
%               lambda >= 0 is the dictionary regularization parameter
%               
%               mu >= 0 is the sparsity regularization parameter.
%               Alternatively, mu can be a 1 x nIters vector specifying a
%               different mu value for each outer iteration
%               
%               gamma in [0, 1] is the forgetting factor, where gamma = 0
%               is memoryless, and gamma = 1 is infinite memory
%               
%               [OPTIONAL] opts is a struct containing one or more of the
%               following fields. The default values are in ()
%                   
%                   opts.A (1) is an m x p system matrix
%                   
%                   opts.M (1) is an m x n weight matrix
%                   
%                   opts.xdim (size(X0)) are the dimensions of Xt to use
%                   when extracting patches. Note that xdim need not be the
%                   same as size(Xt) as long as the two representations are
%                   are interchangable via reshape(). For example, you
%                   could specify xdim = [p, n1, n2] where n1 * n2 = n
%                   
%                   opts.pdim ([8, 8]) is a vector with numel(xdim)
%                   elements specifying the patch sizes to extract from Xt
%                   
%                   opts.pgap ([2, 2]) is a vector with numel(xdim) 
%                   elements specfiying the patch strides (i.e., shifts) to
%                   use along each dimension of Xt
%                   
%                   opts.type ('hard') can be {'hard','soft'} and specifies
%                   whether to use hard or soft thresholding to regularize
%                   the sparse codes Bt
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
%                   opts.unitaryD (false) determines whether to enforce a
%                   unitary constraint on the dictionary (true) or not
%                   (false). Note: if unitaryD == true, the dr and ddim
%                   arguments are ignored
%                   
%                   opts.fixedD (false) determines whether to fix the
%                   initial dictionary (true) or learn it (false)
%                   
%                   opts.nIters (50) is the number of outer iterations to
%                   perform
%                   
%                   opts.nIters0 (0) specifies the number of initial outer
%                   iterations (of the nIters total) for which to hold D
%                   and X fixed (regardless of the value of fixedD)
%                   
%                   opts.nItersDB (1) is the number of inner (D,Bt) updates
%                   to perform per outer iteration
%                   
%                   opts.nItersX (5) is the number of inner Xt updates to
%                   perform per outer iteration. Note that, when A = 1, the
%                   exact Xt update is always computed and nItersX is
%                   ignored
%                   
%                   opts.X0 (A' * Yt) is a p x n matrix specifying the
%                   initial Xt iterate
%                   
%                   opts.D0 (dctmtx(b)') is a b x c matrix specifying the
%                   initial dictionary for time t (usually, the dictionary
%                   from time t - 1). Undercomplete (b > c) and
%                   overcomplete (b < c) dictionaries are allowed, unless
%                   unitaryD == true, in which case we must have b >= c
%                   
%                   opts.B0 (zeros(c,d)) is an c x d matrix specifying 
%                   the initial sparse codes for time t
%                   
%                   opts.params ([]) is a struct containing the cumulative
%                   parameters returned by the previous call to this 
%                   function. At time t = 1, omit this field or pass
%                   params = [] and it will be initialized appropriately
%                   
%                   opts.Xtrue (nan) is the ground truth Xt matrix to use 
%                   for NRMSE calculations
%                   
%                   opts.NRMSEfcn (@computeNRMSE) is the function to use
%                   when computing the NRMSE of Xt after each iteration
%                   
%                   opts.accel (true) specifies whether to use Nesterov's
%                   acceleration scheme for the Xt updates
%                   
%                   opts.tau ((0.99 + ~accel) / norm(A)^2) is the step size
%                   parameter for the Xt updates, and should satisfy
%                   
%                       tau <= 1 / norm(A)^2, when accel = true
%                       tau <  2 / norm(A)^2, when accel = false
%                   
%                   opts.flag (1) determines what status updates to print
%                   to the command window. The choices are
%                   
%                       flag = 0: no printing
%                       flag = 1: print outer iteration status
%                       flag = 2: print inner and outer iteration status
%               
% Outputs:      Xt is the estimated p x n matrix for time t
%               
%               D is a b x c matrix, where b = prod(pdim), c = size(D0,2), 
%               whose columns contain the dictionary atoms for time t
%               
%               Bt is a c x d matrix, where d = total # patches, whose
%               columns are sparse codes for the patches of Xt at time t
%               
%               params is a struct of cumulative parameters for time t (to
%               be passed back into this function at time t + 1)
%               
%               stats is a struct containing the following fields:
%               
%                   stats.nIters is the number of outer iterations
%                   performed
%                   
%                   stats.cost is a 1 x nIters vector containing the cost
%                   function (only the portion affected by Xt, D, and Bt)
%                   values at each iteration 
%                   
%                   stats.ccost is a 1 x nIters vector containing the full
%                   (cumulative) cost function values at each iteration 
%                   
%                   stats.nrmse is a 1 x nIters vector containing the NRMSE
%                   of Xt with respect to Xtrue at each iteration
%                   
%                   stats.deltaX is a 1 x nIters vector containing the
%                   relative convergence of Xt at each iteration, defined:
%                   \|Xt_{k + 1} - Xt_k\|_F / \|Xt_k\|_F
%                   
%                   stats.deltaD is a 1 x nIters vector containing the
%                   relative convergence of D at each iteration, defined as
%                   \|D_{k + 1} - D_k\|_F / \|D_k\|_F
%                   
%                   stats.deltaB is a 1 x nIters vector containing the
%                   relative convergence of Bt at each iteration, defined:
%                   \|Bt_{k + 1} - Bt_k\|_F / \|Bt_k\|_F
%                   
%                   stats.sparsity is a 1 x nIters vector containing the
%                   sparsity, in percent, of Bt at each iteration
%                   
%                   stats.time is a 1 x nIters vector containing the time,
%                   in seconds, required to perform each outer iteration
%               
% Description:  Solves one of the following Online Dictionary Least-Sqaures
%               (Online DLS) problems:
%               
%               When type == 'hard':
%               
%       \min_{Xt,D,Bt} 0.5 \|\sqrt{M} \odot (Yt - A Xt)\|_F^2 + 
%             \lambda (0.5 \sum_{j=1}^t gamma^{t-j} \|Pj(Xj) - D Bj\|_F^2 +
%                      0.5 \mu^2 \|Bt\|_0)
%               
%               When type == 'soft':
%               
%       \min_{Xt,D,Bt} 0.5 \|\sqrt{M} \odot (Yt - A Xt)\|_F^2 + 
%             \lambda (0.5 \sum_{j=1}^t gamma^{t-j} \|Pj(Xj) - D Bj\|_F^2 +
%                          \mu \|Bt\|_1)
%               
%               subject to
%                   \|D(:,k)\|_2 = 1             (unitaryD == false)
%                   rank(R(D(:,k)) <= dr         (~isnan(dr))
%               or
%                   D' D = I                     (unitaryD == true)
%               
% Author:       Brian Moore
%               brimoor@umich.edu
%               
% Date:         August 22, 2017
%

% Parse inputs
[A, M, xdim, pdim, pgap, type, dr, ddim, unitaryD, fixedD, nIters, ...
 nIters0, nItersDB, nItersX, Xt, D, Bt, params, Xtrue, NRMSEfcn, ...
                        accel, tau, flag] = parseInputs(Yt,mu,varargin{:});

T = opts.xdim(3);
M = opts.M;

Xi = opts.X0;


% Data
M1 = M(:,:,1:T);
Xi1 = Xi(:,:,1:T);
Xtrue1 = Xtrue(:,:,1:T);
Y1 = Y(:,:,1:T);

% First minibatch
opts.M         = double(M1);
opts.xdim      = [ny, nx, T];
opts.dr        = dr;
opts.nIters    = vars.nItersi;
opts.X0        = Xi1;
opts.Xtrue     = Xtrue1;
[Xt, D, Bt, params, stats] = onlineDls(Y1,lambda,mu,gamma,opts);

% Initialize reconstruction
Xhat = Xi;
Xhat(:,:,1:T) = Xt;

% Remaining minibatches
tt = 1:dt:(nt + 1 - T);
ni = numel(tt);
for i = 2:ni
    % Data
    t = tt(i);
    Mt = M(:,:,t:(t + T - 1));
    Xit = Xi(:,:,t:(t + T - 1));
    Xtruet = Xtrue(:,:,t:(t + T - 1));
    Yt = Y(:,:,t:(t + T - 1));
    
    % onlineDls (t)
    opts.M      = double(Mt);
    opts.nIters = vars.nIters;
    opts.X0     = interpFill(Xhat,Xit,t,dt);
    opts.D0     = D;
    opts.B0     = Bt;
    opts.params = params;
    opts.Xtrue  = Xtruet;
    [Xt, D, Bt, params, statst] = onlineDls(Yt,lambda,mu,gamma,opts);
    
    % Combine stats
    stats.nIters   = stats.nIters + vars.nIters;
    stats.cost     = [stats.cost, statst.cost];
    stats.nrmse    = [stats.nrmse, statst.nrmse];
    stats.deltaX   = [stats.deltaX, statst.deltaX];
    stats.deltaD   = [stats.deltaD, statst.deltaD];
    stats.deltaB   = [stats.deltaB, statst.deltaB];
    stats.sparsity = [stats.sparsity, statst.sparsity];
    stats.time     = [stats.time, statst.time];
    
    % Update reconstruction
    Xhat = updateRecon(Xhat,Xt,gamma,t,dt);
end
