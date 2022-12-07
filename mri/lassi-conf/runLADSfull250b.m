function runLADSfull250b(idx)
% Syntax: runLADSfull250b(idx);
%         runLADSfull250b();

%{
% Knobs
% L + S init: 250 iteraions w/ lL = 1.1955 and lS = 0.01
inpath    = 'otazo_full.mat';
initpath  = 'otazo_R8_lps_visual.mat';
outpath   = 'otazo_full_lads250/data.mat';
ps        = [1/64, 1/32, 1/16, 1/8, 1/4, 1/2];
seed      = 1;
lambdaL   = 0.5;
lambdaS   = 0.01;
lambdaB   = 0.03;
lambdaB0  = 0.03;
p         = 3;
nIters    = 50;
%}

% Knobs
% L + S init: 250 iteraions w/ lL = 0.525 and lS = 0.01
inpath    = 'otazo_full.mat';
initpath  = 'otazo_R8_lps_mse.mat';
outpath   = 'otazo_full_lads250b/data.mat';
ps        = [1/64, 1/32, 1/16, 1/8, 1/4, 1/2];
seed      = 1;
lambdaL   = 0.5;
lambdaS   = 0.01;
lambdaB   = 0.04;
lambdaB0  = 0.04;
p         = 3;
nIters    = 50;

%{
% Knobs
% L + S init: 20 iteraions w/ lL = 1.1955 and lS = 0.01
inpath    = 'otazo_full.mat';
initpath  = 'otazo_R8_lps_visual20.mat';
outpath   = 'otazo_full_lads20/data.mat';
ps        = [1/64, 1/32, 1/16, 1/8, 1/4, 1/2];
seed      = 1;
lambdaL   = 0.5;
lambdaS   = 0.01;
lambdaB   = 0.03;
lambdaB0  = 0.06;
p         = 3;
nIters    = 50;
%}

% Get parameters for this iteration
nP = numel(ps);
nS = numel(seed);
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
[ii,jj] = ind2sub([nP nS],idx);
psidx   = ps(ii);
seedidx = seed(jj);

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

% Load full data
data = load(inpath);

% Load L + S initialization (Lhat, Shat, lambdaL, lambdaS, nIters)
init = load(initpath);
L0   = init.Lhat;
S0   = init.Shat;

% Generate data
[Y, A, ~, Xtrue, ~] = generateCardiacPerfData(psidx,inf,seedidx,data);

% LADS
lB = logspace(log10(lambdaB0),log10(lambdaB),nIters);
[Lhat,Shat,~,~,ite,time,cost,mse,delta,sparsity] = runLADS(A,Y,Xtrue,L0,S0,lambdaL,lambdaS,lB,p,nIters); %#ok

% Save results
save(out,'Lhat','Shat','ite','time','cost','mse','delta','sparsity');
fprintf('DONE\n');
