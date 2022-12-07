function invivo_LASSI_R8_Diter_run(parFcn,idx)
% Syntax: invivo_LASSI_R8_Diter_run();
%         invivo_LASSI_R8_Diter_run(parFcn,idx);

% Get parameters
[vars, opts] = parFcn();

% Get parameters for this iteration
nI = numel(vars.nIters);
nT = nI;
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
[ii] = ind2sub([nI],idx);
r        = vars.r;
lambdaL  = vars.lambdaL;
lambdaS  = vars.lambdaS;
lambdaB  = vars.lambdaB;
lambdaB0 = vars.lambdaB0;
dr       = vars.dr;
nIters   = vars.nIters(ii);
nItersDB = vars.nItersDB(ii);
nItersLS = vars.nItersLS(ii);

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

% Load undersampled data (Y, mask, Xfft, Xtrue)
data = load(vars.inpath);
[ny, nx, nt] = size(data.Xtrue);
A = Afft(data.mask,[ny, nx, nt]);

% Load L + S initialization (Lhat, Shat, lambdaL, lambdaS, nIters)
init = load(vars.initpath);

% Run LASSI
opts.A        = A;
opts.r        = r;
opts.sdim     = [ny, nx, nt];
opts.dr       = dr;
opts.nIters   = nIters;     % Custom # outer iterations
opts.nItersDB = nItersDB;   % Custom # inner (D,B) iterations
opts.nItersLS = nItersLS;   % Custom # inner (L,S) iterations
opts.L0       = reshape(init.Lhat,[],nt);
opts.S0       = reshape(init.Shat,[],nt);
opts.Xtrue    = reshape(data.Xtrue,[],nt);
lB = logspace(log10(lambdaB0),log10(lambdaB),nIters);
[Lhat, Shat, ~, ~, stats] = drpca(data.Y,lambdaL,lambdaS,lB,opts); %#ok

% Aggregate params
params.r        = r;
params.lambdaL  = lambdaL;
params.lambdaS  = lambdaS;
params.lambdaB  = lambdaB;
params.lambdaB0 = lambdaB0;
params.dr       = dr;
params.nIters   = nIters; 
params.nItersDB = nItersDB;
params.nItersLS = nItersLS; %#ok

% Save results
save(out,'Lhat','Shat','params','stats');
fprintf('DONE\n');
