function pincat_RPCA_full_run(parFcn,idx)
% Syntax: pincat_RPCA_full_run();
%         pincat_RPCA_full_run(parFcn,idx);

% Get parameters
[vars, opts] = parFcn();

% Get parameters for this iteration
nL = numel(vars.nLines);
nS = numel(vars.seed);
nT = nL * nS;
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
[ii, jj] = ind2sub([nL, nS],idx);
nLines   = vars.nLines(ii);
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
addpath('./deps_rpca');

% Generate undersampled data
if ~isfield(vars,'SNR')
    % No noise, by default
    vars.SNR = inf;
end
[Y, A, ~, Xtrue, Xfft] = generateInvivoData(nLines,vars.SNR,seed,vars.inpath);
[ny, nx, nt] = size(Xtrue);
T = TempFFT(3,[ny, nx, nt]);

% Run RPCA
opts.A      = A;
opts.T      = T;
opts.r      = vars.r;
opts.nIters = vars.nIters;
opts.L0     = reshape(Xfft,[],nt);
opts.S0     = zeros(ny * nx,nt);
opts.Xtrue  = reshape(Xtrue,[],nt);
[Lhat, Shat, stats] = robustPCA(Y,vars.lambdaL,vars.lambdaS,opts); %#ok

% Aggregate params
params.nLines = nLines;
params.seed   = seed; %#ok

% Save results
save(out,'Lhat','Shat','params','stats');
fprintf('DONE\n');
