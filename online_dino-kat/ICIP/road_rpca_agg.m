function road_rpca_agg(parFcn)
% Syntax: road_rpca_agg(parFcn);

% Get parameters
[vars, ~] = parFcn();

% Initialize storage
nP = numel(vars.p);
nS = numel(vars.seed);
nL = numel(vars.lambda);
nM = numel(vars.mu);
nIters = vars.nIters;
nt     = vars.idim(3);
stats.cost  = nan(nP,nS,nL,nM,nIters);
stats.nrmse = nan(nP,nS,nL,nM,nIters);
stats.delta = nan(nP,nS,nL,nM,nIters);
stats.time  = nan(nP,nS,nL,nM,nIters);
err.NRMSE  = nan(nP,nS,nL,nM,1);
err.pSNR   = nan(nP,nS,nL,nM,1);
err.pNRMSE = nan(nP,nS,nL,nM,nt);
err.ppSNR  = nan(nP,nS,nL,nM,nt);

% Load data
files  = findMatchingFiles(vars.rawpath);
nFiles = numel(files);
for i = 1:nFiles
    itimer = tic;
    
    % Extract index
    [~, namei, ~] = fileparts(files{i});
    idx = str2double(regexp(namei,'\d+','match'));
    [ii, jj, kk, ll] = ind2sub([nP, nS, nL, nM],idx);
    
    % Save data
    datai = load(files{i},'stats','err');
    stats.cost(ii,jj,kk,ll,:)  = datai.stats.cost;
    stats.nrmse(ii,jj,kk,ll,:) = datai.stats.nrmse;
    stats.delta(ii,jj,kk,ll,:) = datai.stats.delta;
    stats.time(ii,jj,kk,ll,:)  = datai.stats.time;
    err.NRMSE(ii,jj,kk,ll,:)  = datai.err.NRMSE;
    err.pSNR(ii,jj,kk,ll,:)   = datai.err.pSNR;
    err.pNRMSE(ii,jj,kk,ll,:) = datai.err.pNRMSE;
    err.ppSNR(ii,jj,kk,ll,:)  = datai.err.ppSNR;
    
    % Display progress
    fprintf('File %i/%i processed [Time = %.2fs]\n',i,nFiles,toc(itimer));
end

% Load min-NRMSE reconstruction
[~, im] = nanmin(err.NRMSE(:));
datam   = load(files{im},'Lhat','Shat');
Lhat    = datam.Lhat; %#ok
Shat    = datam.Shat; %#ok

% Save results
save(vars.outpath,'vars','stats','err','Lhat','Shat');
fprintf('DONE\n');
