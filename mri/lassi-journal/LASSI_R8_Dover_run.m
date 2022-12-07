function LASSI_R8_Dover_run(parFcn,idx)
% Syntax: LASSI_R8_Dover_run();
%         LASSI_R8_Dover_run(parFcn,idx);

% Get parameters
[vars, opts] = parFcn();

% Get parameters for this iteration
nP = numel(vars.np);
nT = nP;
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
[ii] = ind2sub([nP],idx);
r        = vars.r;
lambdaL  = vars.lambdaL;
lambdaS  = vars.lambdaS;
lambdaB  = vars.lambdaB(min(ii,numel(vars.lambdaB)));
lambdaB0 = vars.lambdaB0(min(ii,numel(vars.lambdaB0)));
dr       = vars.dr;
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
addpath('./deps_lassi');

% Load undersampled data (Y, mask, samp, Xfft, Xtrue)
data = load(vars.inpath);
[ny, nx, nt] = size(data.Xtrue);
A = Emat_xyt(data.mask,data.samp,[ny, nx, nt]);

% Load L + S initialization (Lhat, Shat, lambdaL, lambdaS, nIters)
init = load(vars.initpath);

% Run LASSI
opts.A       = A;
opts.r       = r;
opts.sdim    = [ny, nx, nt];
opts.np      = np; % Custom dictionary size
opts.dr      = dr;
opts.nIters  = vars.nIters;
opts.L0      = reshape(init.Lhat,[],nt);
opts.S0      = reshape(init.Shat,[],nt);
opts.Xtrue   = reshape(data.Xtrue,[],nt);
lB = logspace(log10(lambdaB0),log10(lambdaB),vars.nIters);
[Lhat, Shat, ~, ~, stats] = drpca(data.Y,lambdaL,lambdaS,lB,opts); %#ok

% Aggregate params
params.r        = r;
params.lambdaL  = lambdaL;
params.lambdaS  = lambdaS;
params.lambdaB  = lambdaB;
params.lambdaB0 = lambdaB0;
params.dr       = dr;
params.np       = np; %#ok

% Save results
save(out,'Lhat','Shat','params','stats');
fprintf('DONE\n');
