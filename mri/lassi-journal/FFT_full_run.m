function FFT_full_run(parFcn,idx)
% Syntax: FFT_full_run();
%         FFT_full_run(parFcn,idx);

% Get parameters
vars = parFcn();

% Get parameters for this iteration
nP = numel(vars.ps);
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
ps   = vars.ps(ii);
seed = vars.seed(jj);

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

% Generete undersampled data
ftimer = tic();
[~, ~, ~, Xtrue, Xfft] = generateCardiacPerfData(ps,inf,seed,vars.inpath);

% Compute FFT stats
stats.nrmse = computeNRMSE(Xfft,Xtrue);
stats.time  = toc(ftimer);

% Aggregate params
params.ps   = ps;
params.seed = seed; %#ok

% Save results
save(out,'Xfft','params','stats');
fprintf('DONE\n');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Local functions
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Compute NRMSE
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function err = computeNRMSE(Xhat,Xtrue)
denom = norm(Xtrue(:));
if isnan(denom)
    err = nan;
elseif denom == 0
    err = 0;
else
    err = norm(Xhat(:) - Xtrue(:)) / denom;
end
