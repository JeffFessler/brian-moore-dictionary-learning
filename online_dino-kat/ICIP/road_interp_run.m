function road_interp_run(parFcn,idx)
% Syntax: road_interp_run();
%         road_interp_run(parFcn,idx);

% Get parameters
[vars, opts] = parFcn();

% Get parameters for this iteration
nP = numel(vars.p);
nS = numel(vars.seed);
nT = nP * nS;
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
[ii, jj] = ind2sub([nP, nS],idx);
p        = vars.p(ii);
seed     = vars.seed(jj);

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
addpath('./deps_dls');

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

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% interpVideo
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
itimer = tic();
Xhat = interpVideo(Y,M,opts.algo,opts.method);
stats.time = toc(itimer);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Compute error metrics
err = computeErrorMetrics(Xhat,Xtrue); %#ok

% Aggregate params
params.p    = p;
params.seed = seed; %#ok

% Save results
save(out,'Xhat','params','err','stats');
fprintf('DONE\n');
