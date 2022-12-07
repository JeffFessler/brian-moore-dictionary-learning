function pincat_KTSLR_R8_run(parFcn,idx)
% Syntax: pincat_KTSLR_R8_run();
%         pincat_KTSLR_R8_run(parFcn,idx);

% Get parameters
[vars, opts] = parFcn();

% Get parameters for this iteration
nP = numel(vars.p);
n1 = numel(vars.mu1);
n2 = numel(vars.mu2);
nT = nP * n1 * n2;
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
[ii, jj, kk] = ind2sub([nP, n1, n2],idx);
p   = vars.p(ii);
mu1 = vars.mu1(jj);
mu2 = vars.mu2(kk);

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

% Load undersampled data (Y, mask, Xfft, Xtrue)
data = load(vars.inpath);
[ny, nx, nt] = size(data.Xtrue);
A = Afft(data.mask,[ny, nx, nt]);

% Run k-t SLR
opts.nItersO  = vars.nItersO;
opts.nItersI  = vars.nItersI;
opts.nItersCG = vars.nItersCG;
opts.Xtrue    = data.Xtrue;
[Xhat, stats] = ktSLR(data.Y,A,p,mu1,mu2,data.Xfft,opts); %#ok

% Aggregate params
params.mu1      = mu1;
params.mu2      = mu2;
params.p        = p; %#ok

% Save results
save(out,'Xhat','params','stats');
fprintf('DONE\n');
