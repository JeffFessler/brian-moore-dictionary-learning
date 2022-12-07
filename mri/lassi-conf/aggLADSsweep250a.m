function aggLADSsweep250a()
% Syntax: aggLADSsweep250a();

%{
% Knobs
% L + S init: 250 iteraions w/ lL = 1.1955 and lS = 0.01
inpath   = 'otazo_R8_lads250/data.mat';
outpath  = 'otazo_R8_lads_mse250.mat';
lambdaL  = [0.5, 0.7, 0.9, 1.2, 1.5];
lambdaS  = [0.005, 0.01, 0.02, 0.03, 0.05];
lambdaB  = [0.02, 0.03, 0.04, 0.07];
lambdaB0 = [0.02, 0.03, 0.04, 0.07];
p        = 3;
nIters   = 50;
%}

% Knobs
% L + S init: 250 iteraions w/ lL = 1.1955 and lS = 0.01
inpath   = 'otazo_R8_lads250a/data.mat';
outpath  = 'otazo_R8_lads_mse250a.mat';
lambdaL  = [0.05, 0.1, 0.3, 0.5, 0.7, 0.9];
lambdaS  = [0.005, 0.01, 0.02, 0.03];
lambdaB  = [0.02, 0.03, 0.04];
lambdaB0 = [0.02, 0.03, 0.04];
p        = 3;
nIters   = 50;

%{
% Knobs
% L + S init: 250 iteraions w/ lL = 0.525 and lS = 0.01
inpath   = 'otazo_R8_lads250b/data.mat';
outpath  = 'otazo_R8_lads_mse250b.mat';
lambdaL  = [0.05, 0.1, 0.3, 0.5, 0.7, 0.9];
lambdaS  = [0.005, 0.01, 0.02, 0.03];
lambdaB  = [0.02, 0.03, 0.04];
lambdaB0 = [0.02, 0.03, 0.04];
p        = 3;
nIters   = 50;
%}

%{
% Knobs
% L + S init: 20 iteraions w/ lL = 1.1955 and lS = 0.01
inpath   = 'otazo_R8_lads20/data.mat';
outpath  = 'otazo_R8_lads_mse20.mat';
lambdaL  = [0.05, 0.1, 0.3, 0.5, 0.7, 0.9];
lambdaS  = [0.005, 0.01, 0.02, 0.03];
lambdaB  = [0.02, 0.03, 0.04];
lambdaB0 = [0.04, 0.06, 0.08];
p        = 3;
nIters   = 50;
%}

% Initialize storage
nL = numel(lambdaL);
nS = numel(lambdaS);
nB = numel(lambdaB);
time     = nan(nL,nS,nB,nIters);
cost     = nan(nL,nS,nB,nIters);
mse      = nan(nL,nS,nB,nIters);
delta    = nan(nL,nS,nB,nIters);
sparsity = nan(nL,nS,nB,nIters);

% Load data
files  = findMatchingFiles(inpath);
nFiles = numel(files);
for i = 1:nFiles
    % Extract index
    [~,namei,~] = fileparts(files{i});
    idx         = str2double(regexp(namei,'\d+','match'));
    [ii,jj,kk]  = ind2sub([nL nS nB],idx);
    
    % Save data
    datai = load(files{i},'ite','time','cost','mse','delta','sparsity');
    time(ii,jj,kk,:)     = datai.time;
    cost(ii,jj,kk,:)     = datai.cost;
    mse(ii,jj,kk,:)      = datai.mse;
    delta(ii,jj,kk,:)    = datai.delta;
    sparsity(ii,jj,kk,:) = datai.sparsity;
end

% Load optimal reconstruction
mseFinal   = mse(:,:,:,nIters);
[~,idxopt] = nanmin(mseFinal(:));
[pathopt,fileopt,extopt] = fileparts(inpath);
inpathopt = sprintf('%s/%s%d%s',pathopt,fileopt,idxopt,extopt);
dataopt   = load(inpathopt);
Lhat = dataopt.Lhat; %#ok
Shat = dataopt.Shat; %#ok

% Save results
save(outpath,'Lhat','Shat','time','cost','mse','delta', ...
             'sparsity','lambdaL','lambdaS','lambdaB','lambdaB0','p','nIters');
fprintf('DONE\n');
