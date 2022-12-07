function onlineDls_run(parFcn,idx)
%Syntax: onlineDls_run();
%        onlineDls_run(parFcn,idx);

% Get parameters
[vars, opts, vars0, opts0] = parFcn();
IS_OTAZO   = isfield(vars,'ps');
INIT_RPCA  = isfield(vars0,'lambdaL');
INIT_KTSLR = isfield(vars0,'mu1');

% Get parameters for this iteration
if IS_OTAZO
    nP = numel(vars.ps);
else
    nP = numel(vars.nLines);
end
nS = numel(vars.seed);
nL = numel(vars.lambda);
nM = numel(vars.mu);
nR = numel(vars.dr);
nG = numel(vars.gamma);
nT = nP * nS * nL * nM * nR * nG;
if ~exist('idx','var') || isempty(idx)
    % Return job count
    fprintf('Total jobs: %d\n',nT);
    return;
elseif (idx < 1) || (idx > nT)
    % Invalid index
    fprintf('Index %d is out of range [1,%d]\n',idx,nT);
    return;
else
    % Valid index
    fprintf('Running index %d/%d\n',idx,nT);
end
[ii, jj, kk, ll, mm, nn] = ind2sub([nP, nS, nL, nM, nR, nG],idx);
if IS_OTAZO
    ps = vars.ps(ii);
else
    nLines = vars.nLines(ii);
end
seed   = vars.seed(jj);
lambda = vars.lambda(kk);
mu     = vars.mu(ll);
mu2    = vars.mu2(ll); % coupled
dr     = vars.dr(mm);
gamma  = vars.gamma(nn);
gamma2 = vars.gamma2(nn); % coupled

% Parse path
[path, name, ext] = fileparts(vars.rawpath);
if ~exist(path,'dir')
    % Make results directory
    mkdir(path);
end
out = sprintf('%s/%s%d%s',path,name,idx,ext);
if exist(out,'file')
    % Output data already exists
    fprintf('Output file ''%s'' already exists\n',out);
    return;
end

% Set random seed
rng(seed);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Load data
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if IS_OTAZO
    % Otazo data
    [Y, Xtrue, Xfft, b1] = generateCardiacPerfData2(ps,vars.SNR,vars.inpath);
    M = (abs(Y(:,:,:,1)) ~= 0);
    A = Emat_xyt(M,b1,size(M));
else
    % Invivo/PINCAT data
    [Y, Xtrue, Xfft] = generateInvivoData2(nLines,vars.SNR,vars.inpath);
    M = (abs(Y) ~= 0);
    A = Afft(M,size(M));
end
[ny, nx, nt] = size(Xtrue);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Initialization
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if INIT_RPCA
    % Add dependencies to path
    addpath('./deps_rpca');
    
    % Run RPCA
    opts0.A      = A;
    opts0.T      = TempFFT(3,[ny, nx, nt]);
    opts0.r      = vars0.r;
    opts0.nIters = vars0.nIters;
    opts0.L0     = reshape(Xfft,[],nt);
    opts0.S0     = zeros(ny * nx,nt);
    opts0.Xtrue  = reshape(Xtrue,[],nt);
    [L0, S0] = robustPCA(Y,vars0.lambdaL,vars0.lambdaS,opts0);
    X0 = reshape(L0 + S0,[ny, nx, nt]);
elseif INIT_KTSLR
    % Add dependencies to path
    addpath('./deps_ktslr');
    
    % Run k-t SLR
    opts0.nItersO  = vars0.nItersO;
    opts0.nItersI  = vars0.nItersI;
    opts0.nItersCG = vars0.nItersCG;
    X0 = ktSLR(Y,A,vars0.p,vars0.mu1,vars0.mu2,Xfft,opts0);
    X0 = reshape(X0,[ny, nx, nt]);
else
    % Baseline
    X0 = Xfft;
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% onlineDls
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Add dependencies to path
addpath('./deps_dls');

% Relevant variables
T = vars.T;
dt = vars.dt;
pdim = opts.pdim;
pgap = opts.pgap;
D0 = opts.D0;

% Modify initial dictionary, if requested
if vars.np > 0
    % Append atoms (= random patches from first T frames of X0)
    Mp = patchInds([ny, nx, T],pdim,pgap);
    m = size(Mp,2);
    idx = randperm(m);
    inds = Mp(:,idx(1:min(vars.np,m)));
    Dnew = X0(inds);
    Dnew = bsxfun(@rdivide,Dnew,sqrt(sum(abs(Dnew).^2,1)));
    D0 = [D0, Dnew];
elseif vars.np < 0
    % Remove random atoms
    m = size(D0,2);
    idx = randperm(m);
    D0 = D0(:,idx(1:max(1,m + vars.np)));
end

% Initialization
Xhat = X0;
stats = struct();
stats.nIters   = 0;
stats.cost     = [];
stats.nrmse    = [];
stats.deltaX   = [];
stats.deltaD   = [];
stats.deltaB   = [];
stats.sparsity = [];
stats.time     = [];

% Multiple online passes
for kk = 1:vars.nReps
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % t = 1
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Data
    M1 = M(:,:,1:T);
    if IS_OTAZO
        Y1 = Y(:,:,1:T,:);
        A1 = Emat_xyt(M1,b1,[ny, nx, T]);
    else
        Y1 = Y(:,:,1:T);
        A1 = Afft(M1,[ny, nx, T]);
    end
    X01 = Xhat(:,:,1:T);
    Xtrue1 = Xtrue(:,:,1:T);
    if mu2 == 0
        % Initialize with zeros
        B01 = nan;
    else
        % Initialize with fixed sparsity
        sparsity = abs(mu2);
        B01 = sparseCodesInit(D0,X01,pdim,pgap,sparsity);
    end
    
    % onlineDls
    opts.A      = A1;
    opts.xdim   = [ny, nx, T];
    opts.dr     = dr;
    opts.nIters = vars.nItersi;
    opts.X0     = reshape(X01,[],T);
    opts.D0     = D0;
    opts.B0     = B01;
    opts.Xtrue  = reshape(Xtrue1,[],T);
    [Xt, D, Bt, params, stats1] = onlineDls(Y1,lambda,mu,gamma,opts);
    Xt = reshape(Xt,[ny, nx, T]);
    
    % Combine stats
    stats.nIters   = stats.nIters + vars.nItersi;
    stats.cost     = [stats.cost, stats1.cost];
    stats.nrmse    = [stats.nrmse, stats1.nrmse];
    stats.deltaX   = [stats.deltaX, stats1.deltaX];
    stats.deltaD   = [stats.deltaD, stats1.deltaD];
    stats.deltaB   = [stats.deltaB, stats1.deltaB];
    stats.sparsity = [stats.sparsity, stats1.sparsity];
    stats.time     = [stats.time, stats1.time];
    
    % Update reconstruction
    Xhat(:,:,1:T) = Xt;
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % t = 2,3,...
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    tt = 1:dt:(nt + 1 - T);
    ni = numel(tt);
    for i = 2:ni
        % Data
        t = tt(i);
        Mt = M(:,:,t:(t + T - 1));
        if IS_OTAZO
            Yt = Y(:,:,t:(t + T - 1),:);
            At = Emat_xyt(Mt,b1,[ny, nx, T]);
        else
            Yt = Y(:,:,t:(t + T - 1));
            At = Afft(Mt,[ny, nx, T]);
        end
        X0t = Xhat(:,:,t:(t + T - 1));
        Xtruet = Xtrue(:,:,t:(t + T - 1));
        if mu2 < 0
            % Initialize with last sparse codes
            B0t = Bt;
        else
            % Initialize with fixed sparsity
            %sparsity = 100 * nnz(Bt) / numel(Bt);
            sparsity = mu2;
            B0t = sparseCodesInit(D,X0t,pdim,pgap,sparsity);
        end
        
        % onlineDls
        opts.A      = At;
        opts.nIters = vars.nIters;
        opts.X0     = reshape(X0t,[],T);
        opts.D0     = D;
        opts.B0     = B0t;
        opts.params = params;
        opts.Xtrue  = reshape(Xtruet,[],T);
        [Xt, D, Bt, params, statst] = onlineDls(Yt,lambda,mu,gamma,opts);
        Xt = reshape(Xt,[ny, nx, T]);
        
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
        Xhat = updateRecon(Xhat,Xt,gamma2,t,dt);
    end
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Compute error metrics
err = computeErrorMetrics(Xhat,Xtrue); %#ok

% Aggregate params
params = struct();
if IS_OTAZO
    params.ps = ps;
else
    params.nLines = nLines;
end
params.seed   = seed;
params.lambda = lambda;
params.mu     = mu;
params.mu2    = mu2;
params.dr     = dr;
params.gamma  = gamma;
params.gamma2 = gamma2; %#ok

% Save results
save(out,'Xhat','D','params','err','stats');
fprintf('DONE\n');
