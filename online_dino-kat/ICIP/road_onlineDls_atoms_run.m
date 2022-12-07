function road_onlineDls_atoms_run(parFcn,idx)
% Syntax: road_onlineDls_atoms_run();
%         road_onlineDls_atoms_run(parFcn,idx);

% Get parameters
[vars, opts] = parFcn();

% Get parameters for this iteration
nN = numel(vars.np);
nT = nN;
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
[ii] = ind2sub([nN],idx);
p        = vars.p;
seed     = vars.seed;
lambda   = vars.lambda;
mu       = vars.mu;
dr       = vars.dr;
gamma    = vars.gamma;
gamma2   = vars.gamma2;
np       = vars.np(ii);

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
dt = vars.dt;
T = vars.T;
Xtrue = loadVideo(vars.inpath,vars.idim,T,dt);
[ny, nx, nt] = size(Xtrue);

% Missing data
M = (rand(ny,nx,nt) > p);
Y = Xtrue;
Y(~M) = 0;

% Interpolated video
Xi = interpVideo(Y,M,'paint2',1);

% Modify initial dictionary, if requested
if np > 0
    % Append atoms (= random patches from first T frames of Xi)
    Mp = patchInds([ny, nx, T],opts.pdim,opts.pgap);
    m = size(Mp,2);
    idx  = randperm(m);
    inds = Mp(:,idx(1:min(np,m)));
    Dnew = Xi(inds);
    Dnew = bsxfun(@rdivide,Dnew,sqrt(sum(abs(Dnew).^2,1)));
    opts.D0 = [opts.D0, Dnew];
elseif np < 0
    % Remove random atoms
    m  = size(opts.D0,2);
    idx = randperm(m);
    opts.D0 = opts.D0(:,idx(1:max(1,m + np)));
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% onlineDls
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Data
M1 = M(:,:,1:T);
Xi1 = Xi(:,:,1:T);
Xtrue1 = Xtrue(:,:,1:T);
Y1 = Y(:,:,1:T);

% onlineDls (t = 1)
opts.M         = M1;
opts.xdim      = [ny, nx, T];
opts.dr        = dr;
opts.nIters    = vars.nItersi;
opts.X0        = Xi1;
opts.Xtrue     = Xtrue1;
[Xt, D, Bt, params, stats] = onlineDls(Y1,lambda,mu,gamma,opts);

% Initialize reconstruction
Xhat = zeros(ny,nx,nt);
Xhat(:,:,1:T) = Xt;

% t = 2,3,...
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
    opts.M      = Mt;
    opts.nIters = vars.nIters;
    %opts.X0    = zeroFill(Xhat,Yt,Mt,t,dt);
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
    Xhat = updateRecon(Xhat,Xt,gamma2,t,dt);
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Compute error metrics
err = computeErrorMetrics(Xhat,Xtrue); %#ok

% Aggregate params
params = struct();
params.np     = np; %#ok

% Save results
save(out,'Xhat','D','params','err','stats');
fprintf('DONE\n');
