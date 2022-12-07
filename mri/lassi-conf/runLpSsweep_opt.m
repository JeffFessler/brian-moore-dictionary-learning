function runLpSsweep_opt(idx)
% Syntax:   runLpSsweep_opt(idx);
%           runLpSsweep_opt();

% Knobs
inpath  = 'otazo_R8.mat';
outpath = 'otazo_R8_lps_opt/data.mat';
r       = 1:10;
lambdaS = 0.01 * logspace(-1,1,25);
nIters  = 250;

% Get parameters for this iteration
nS = numel(lambdaS);
nT = nS;
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
[jj] = ind2sub([nS],idx);
lambdaSidx = lambdaS(jj);

% Parse path
[path,name,ext] = fileparts(outpath);
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

% Load undersampled data
load(inpath);
A = Emat_xyt(mask,samp);

% Reconstructions
nR     = numel(r);
time   = nan(nR,nIters);
cost   = nan(nR,nIters);
mse    = nan(nR,nIters);
delta  = nan(nR,nIters);
mmse   = inf;
Lhat   = nan(size(Xtrue));
Shat   = nan(size(Xtrue));
count  = 0;
for ii = 1:nR
    % L + S
    count = count + 1;
    fprintf('\nSimulation %i/%i\n\n',count,numel(mse) / nIters);
    [Lh,Sh,it,ti,co,ms,de] = runLpS(A,Y,Xtrue,Xfft,zeros(size(Xfft)),{'opt',r(ii)},lambdaSidx,nIters);
    
    % Save results
    time(ii,1:it)  = ti;
    cost(ii,1:it)  = co;
    mse(ii,1:it)   = ms;
    delta(ii,1:it) = de;
    if ms(it) < mmse
        % Update optimal solution
        mmse = ms(it);
        Lhat = Lh;
        Shat = Sh;
    end
end

% Save results
save(out,'Lhat','Shat','time','cost','mse','delta');
fprintf('DONE\n');
