function runLADSsweep250b(idx)
% Syntax: runLADSsweep250b();
%         runLADSsweep250b(idx);

%{
% Knobs
% L + S init: 250 iteraions w/ lL = 1.1955 and lS = 0.01
inpath   = 'otazo_R8.mat';
initpath = 'otazo_R8_lps_visual.mat';
outpath  = 'otazo_R8_lads250/data.mat';
lambdaL  = [0.5, 0.7, 0.9, 1.2, 1.5];
lambdaS  = [0.005, 0.01, 0.02, 0.03, 0.05];
lambdaB  = [0.02, 0.03, 0.04, 0.07];
lambdaB0 = [0.02, 0.03, 0.04, 0.07];
p        = 3;
nIters   = 50;
%}

%{
% Knobs
% L + S init: 250 iteraions w/ lL = 1.1955 and lS = 0.01
inpath   = 'otazo_R8.mat';
initpath = 'otazo_R8_lps_visual.mat';
outpath  = 'otazo_R8_lads250a/data.mat';
lambdaL  = [0.05, 0.1, 0.3, 0.5, 0.7, 0.9];
lambdaS  = [0.005, 0.01, 0.02, 0.03];
lambdaB  = [0.02, 0.03, 0.04];
lambdaB0 = [0.02, 0.03, 0.04];
p        = 3;
nIters   = 50;
%}

% Knobs
% L + S init: 250 iteraions w/ lL = 0.525 and lS = 0.01
inpath   = 'otazo_R8.mat';
initpath = 'otazo_R8_lps_mse.mat';
outpath  = 'otazo_R8_lads250b/data.mat';
lambdaL  = [0.05, 0.1, 0.3, 0.5, 0.7, 0.9];
lambdaS  = [0.005, 0.01, 0.02, 0.03];
lambdaB  = [0.02, 0.03, 0.04];
lambdaB0 = [0.02, 0.03, 0.04];
p        = 3;
nIters   = 50;

%{
% Knobs
% L + S init: 20 iteraions w/ lL = 1.1955 and lS = 0.01
inpath   = 'otazo_R8.mat';
initpath = 'otazo_R8_lps_visual20.mat';
outpath  = 'otazo_R8_lads20/data.mat';
lambdaL  = [0.05, 0.1, 0.3, 0.5, 0.7, 0.9];
lambdaS  = [0.005, 0.01, 0.02, 0.03];
lambdaB  = [0.02, 0.03, 0.04];
lambdaB0 = [0.04, 0.06, 0.08];
p        = 3;
nIters   = 50;
%}

% Get parameters for this iteration
nL = numel(lambdaL);
nS = numel(lambdaS);
nB = numel(lambdaB);
nT = nL * nS * nB;
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
[ii,jj,kk] = ind2sub([nL nS nB],idx);
lambdaLidx = lambdaL(ii);
lambdaSidx = lambdaS(jj);
lambdaBidx = logspace(log10(lambdaB0(kk)),log10(lambdaB(kk)),nIters);

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

% Load undersampled data (Y, mask, samp, Xfft, Xtrue)
load(inpath);
A = Emat_xyt(mask,samp);

% Load L + S initialization (Lhat, Shat, lambdaL, lambdaS, nIters)
init = load(initpath);
L0   = init.Lhat;
S0   = init.Shat;

% Run LADS
[Lhat,Shat,Dhat,Bhat,ite,time,cost,mse,delta,sparsity] = runLADS(A,Y,Xtrue,L0,S0,lambdaLidx,lambdaSidx,lambdaBidx,p,nIters); %#ok

% Save results
save(out,'Lhat','Shat','ite','time','cost','mse','delta','sparsity');
fprintf('DONE\n');
