function invivo_KTSLR_full_run(parFcn,idx)
% Syntax: invivo_KTSLR_full_run();
%         invivo_KTSLR_full_run(parFcn,idx);

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
nLines = vars.nLines(ii);
seed   = vars.seed(jj);

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
addpath('./deps_ktslr');

% Generate undersampled data
[Y, A, ~, Xtrue, Xfft] = generateInvivoData(nLines,inf,seed,vars.inpath);

% Run k-t SLR
opts.nItersO  = vars.nItersO;
opts.nItersI  = vars.nItersI;
opts.nItersCG = vars.nItersCG;
opts.Xtrue    = Xtrue;
[Xhat, stats] = ktSLR(Y,A,vars.p,vars.mu1,vars.mu2,Xfft,opts); %#ok

% Aggregate params
params.nLines = nLines;
params.seed   = seed; %#ok

% Save results
save(out,'Xhat','params','stats');
fprintf('DONE\n');
