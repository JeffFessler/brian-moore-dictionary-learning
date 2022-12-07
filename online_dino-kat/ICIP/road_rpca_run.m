function road_rpca_run(parFcn,idx)
% Syntax: road_rpca_run();
%         road_rpca_run(parFcn,idx);

% Get parameters
[vars, opts] = parFcn();

% Get parameters for this iteration
nP = numel(vars.p);
nS = numel(vars.seed);
nL = numel(vars.lambda);
nM = numel(vars.mu);
nT = nP * nS * nL * nM;
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
[ii, jj, kk, ll] = ind2sub([nP, nS, nL, nM],idx);
p        = vars.p(ii);
seed     = vars.seed(jj);
lambda   = vars.lambda(kk);
mu       = vars.mu(ll);

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
addpath('./deps_rpca');

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
if isfield(vars,'algo') && isfield(vars,'method')
    % Specified method
    Xi = interpVideo(Y,M,vars.algo,vars.method);
else
    % Default method
    Xi = interpVideo(Y,M,'paint2',1);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% robustPCA
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Reshape data
Y2 = reshape(Y,ny * nx,nt);
M2 = reshape(M,ny * nx,nt);
Xi2 = reshape(Xi,ny * nx,nt);
Xtrue2 = reshape(Xtrue,ny * nx,nt);

% robustPCA
lambdaL = lambda;
lambdaS = lambda * mu;
opts.M        = double(M2);
opts.nIters   = vars.nIters;
opts.L0       = Xi2;
opts.S0       = zeros(ny * nx,nt);
opts.Xtrue    = Xtrue2;
[Lhat, Shat, stats] = robustPCA(Y2,lambdaL,lambdaS,opts); %#ok
Lhat = reshape(Lhat,[ny, nx, nt]);
Shat = reshape(Shat,[ny, nx, nt]);
Xhat = Lhat + Shat;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Compute error metrics
err = computeErrorMetrics(Xhat,Xtrue); %#ok

% Aggregate params
params.p      = p;
params.seed   = seed;
params.lambda = lambda;
params.mu     = mu; %#ok

% Save results
save(out,'Lhat','Shat','params','err','stats');
fprintf('DONE\n');
