function SOUP_true_run(parFcn,idx)
% Syntax: SOUP_true_run();
%         SOUP_true_run(parFcn,idx);

% Get parameters
[vars, opts] = parFcn();

% Get parameters for this iteration
nB = numel(vars.lambdaB);
nD = numel(vars.dr);
nT = nB * nD;
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
[ii, jj] = ind2sub([nB, nD],idx);
lambdaB  = vars.lambdaB(ii);
dr       = vars.dr(jj);

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

% Load data (Xtrue)
data = load(vars.inpath);

% Extract patches
sdim = size(data.Xtrue);
Mp   = patchInds(sdim,vars.pdim,vars.pgap);
Y    = data.Xtrue(Mp);

% Run SOUP-DIL
opts.dr     = dr;
opts.nIters = vars.nIters;
opts.B0     = zeros(size(opts.D0,2),size(Y,2));
[Dhat, ~, stats] = soupDil(Y,lambdaB,opts); %#ok

% Aggregate params
params.lambdaB = lambdaB;
params.dr      = dr; %#ok

% Save results
save(out,'Dhat','params','stats');
fprintf('DONE\n');
