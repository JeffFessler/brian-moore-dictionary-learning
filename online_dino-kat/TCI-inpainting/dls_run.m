function dls_run(parFcn,idx)
% Syntax: dls_run();
%         dls_run(parFcn,idx);

% Get parameters
[vars, opts, vars0] = parFcn();

% Get parameters for this iteration
nP = numel(vars.p);
nS = numel(vars.seed);
nL = numel(vars.lambda);
nM = numel(vars.mu);
nR = numel(vars.dr);
nT = nP * nS * nL * nM * nR;
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
[ii, jj, kk, ll, mm] = ind2sub([nP, nS, nL, nM, nR],idx);
p        = vars.p(ii);
seed     = vars.seed(jj);
lambda   = vars.lambda(kk);
mu       = vars.mu(ll);
dr       = vars.dr(mm);

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

% Add dependencies to path
addpath('./deps');

% Set random seed
rng(seed);

% Load data
Xtrue = loadVideo(vars.inpath,vars.idim);
[ny, nx, nt] = size(Xtrue);

% Add noise, if requested
if isfield(vars,'SNR') && isfinite(vars.SNR)
    sigma = 10^(-vars.SNR / 20) * norm(Xtrue(:)) / sqrt(numel(Xtrue));
    Y = Xtrue + sigma * randn(ny,nx,nt);
else
    Y = Xtrue;
end

% Missing data
M = (rand(ny,nx,nt) > p);
Y(~M) = 0;

% Interpolated video
Xi = interpVideo(Y,M,vars0.algo,vars0.method);

% Modify initial dictionary, if requested
if isfield(vars,'np')
    if vars.np > 0
        % Append atoms (= random patches from Xi)
        Mp = patchInds([ny, nx, nt],opts.pdim,opts.pgap);
        m = size(Mp,2);
        idx = randperm(m);
        inds = Mp(:,idx(1:min(vars.np,m)));
        Dnew = Xi(inds);
        Dnew = bsxfun(@rdivide,Dnew,sqrt(sum(abs(Dnew).^2,1)));
        opts.D0 = [opts.D0, Dnew];
    elseif vars.np < 0
        % Remove random atoms
        m = size(opts.D0,2);
        idx = randperm(m);
        opts.D0 = opts.D0(:,idx(1:max(1,m + vars.np)));
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% dls
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
opts.M      = double(M);
opts.xdim   = [ny, nx, nt];
opts.dr     = dr;
opts.nIters = vars.nIters;
opts.X0     = Xi;
opts.Xtrue  = Xtrue;
[Xhat, D, ~, stats] = dls(Y,lambda,mu,opts); %#ok
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Compute error metrics
err = computeErrorMetrics(Xhat,Xtrue); %#ok

% Aggregate params
params.p      = p;
params.seed   = seed;
params.lambda = lambda;
params.mu     = mu;
params.dr     = dr; %#ok

% Save results
%save(out,'Xhat','D','params','err','stats');
save(out,'params','err','stats');
fprintf('DONE\n');
