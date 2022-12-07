function RPCA_R8_run(parFcn,idx)
% Syntax: RPCA_R8_run();
%         RPCA_R8_run(parFcn,idx);

% Get parameters
[vars, opts] = parFcn();

% Get parameters for this iteration
nR = numel(vars.r);
nL = numel(vars.lambdaL);
nS = numel(vars.lambdaS);
nT = nR * nL * nS;
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
[ii, jj, kk] = ind2sub([nR, nL, nS],idx);
r        = vars.r(ii);
lambdaL  = vars.lambdaL(jj);
lambdaS  = vars.lambdaS(kk);

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

% Load undersampled data (Y, mask, samp, Xfft, Xtrue)
data = load(vars.inpath);
[ny, nx, nt] = size(data.Xtrue);
A = Emat_xyt(data.mask,data.samp,[ny, nx, nt]);
T = TempFFT(3,[ny, nx, nt]);

% Run RPCA
opts.A      = A;
opts.T      = T;
opts.r      = r;
opts.nIters = vars.nIters;
opts.L0     = reshape(data.Xfft,[],nt);
opts.S0     = zeros(ny * nx,nt);
opts.Xtrue  = reshape(data.Xtrue,[],nt);
[Lhat, Shat, stats] = robustPCA(data.Y,lambdaL,lambdaS,opts); %#ok

% Aggregate params
params.r        = r;
params.lambdaL  = lambdaL;
params.lambdaS  = lambdaS; %#ok

% Save results
save(out,'Lhat','Shat','params','stats');
fprintf('DONE\n');
