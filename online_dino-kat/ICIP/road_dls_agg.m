function road_dls_agg(parFcn)
% Syntax: road_dls_agg(parFcn);

% Get parameters
[vars, ~] = parFcn();

% Initialize storage
nP = numel(vars.p);
nS = numel(vars.seed);
nL = numel(vars.lambda);
nM = numel(vars.mu);
nR = numel(vars.dr);
nIters = vars.nIters;
nt     = vars.idim(3);
stats.cost     = nan(nP,nS,nL,nM,nR,nIters);
stats.nrmse    = nan(nP,nS,nL,nM,nR,nIters);
stats.deltaX   = nan(nP,nS,nL,nM,nR,nIters);
stats.deltaD   = nan(nP,nS,nL,nM,nR,nIters);
stats.deltaB   = nan(nP,nS,nL,nM,nR,nIters);
stats.sparsity = nan(nP,nS,nL,nM,nR,nIters);
stats.time     = nan(nP,nS,nL,nM,nR,nIters);
err.NRMSE  = nan(nP,nS,nL,nM,nR,1);
err.pSNR   = nan(nP,nS,nL,nM,nR,1);
err.pNRMSE = nan(nP,nS,nL,nM,nR,nt);
err.ppSNR  = nan(nP,nS,nL,nM,nR,nt);

% Load data
files  = findMatchingFiles(vars.rawpath);
nFiles = numel(files);
for i = 1:nFiles
    itimer = tic;
    
    % Extract index
    [~, namei, ~] = fileparts(files{i});
    idx = str2double(regexp(namei,'\d+','match'));
    [ii, jj, kk, ll, mm] = ind2sub([nP, nS, nL, nM, nR],idx);
    
    % Save data
    datai = load(files{i},'stats','err');
    stats.cost(ii,jj,kk,ll,mm,:)     = datai.stats.cost;
    stats.nrmse(ii,jj,kk,ll,mm,:)    = datai.stats.nrmse;
    stats.deltaX(ii,jj,kk,ll,mm,:)   = datai.stats.deltaX;
    stats.deltaD(ii,jj,kk,ll,mm,:)   = datai.stats.deltaD;
    stats.deltaB(ii,jj,kk,ll,mm,:)   = datai.stats.deltaB;
    stats.sparsity(ii,jj,kk,ll,mm,:) = datai.stats.sparsity;
    stats.time(ii,jj,kk,ll,mm,:)     = datai.stats.time;
    err.NRMSE(ii,jj,kk,ll,mm,:)  = datai.err.NRMSE;
    err.pSNR(ii,jj,kk,ll,mm,:)   = datai.err.pSNR;
    err.pNRMSE(ii,jj,kk,ll,mm,:) = datai.err.pNRMSE;
    err.ppSNR(ii,jj,kk,ll,mm,:)  = datai.err.ppSNR;
    
    % Display progress
    fprintf('File %i/%i processed [Time = %.2fs]\n',i,nFiles,toc(itimer));
end

% Load min-NRMSE reconstruction
[~, im] = nanmin(err.NRMSE(:));
datam   = load(files{im},'Xhat','D');
Xhat    = datam.Xhat; %#ok
D       = datam.D; %#ok

% Save results
save(vars.outpath,'vars','stats','err','Xhat','D');
fprintf('DONE\n');
