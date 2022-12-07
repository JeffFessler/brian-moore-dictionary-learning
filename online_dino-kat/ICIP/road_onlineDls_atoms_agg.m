function road_onlineDls_atoms_agg(parFcn)
% Syntax: road_onlineDls_atoms_agg(parFcn);

% Get parameters
[vars, ~] = parFcn();

% Initialize storage
nN = numel(vars.np);
nt = vars.idim(3);
ni = numel(1:vars.dt:(nt - vars.T + 1));
nIters = vars.nItersi + vars.nIters * (ni - 1);
stats.cost     = nan(nN,nIters);
stats.nrmse    = nan(nN,nIters);
stats.deltaX   = nan(nN,nIters);
stats.deltaD   = nan(nN,nIters);
stats.deltaB   = nan(nN,nIters);
stats.sparsity = nan(nN,nIters);
stats.time     = nan(nN,nIters);
err.NRMSE  = nan(nN,1);
err.pSNR   = nan(nN,1);
err.pNRMSE = nan(nN,nt);
err.ppSNR  = nan(nN,nt);

% Load data
files  = findMatchingFiles(vars.rawpath);
nFiles = numel(files);
for i = 1:nFiles
    itimer = tic;
    
    % Extract index
    [~, namei, ~] = fileparts(files{i});
    idx = str2double(regexp(namei,'\d+','match'));
    [ii] = ind2sub([nN],idx);
    
    % Save data
    datai = load(files{i},'stats','err');
    stats.cost(ii,:)     = datai.stats.cost;
    stats.nrmse(ii,:)    = datai.stats.nrmse;
    stats.deltaX(ii,:)   = datai.stats.deltaX;
    stats.deltaD(ii,:)   = datai.stats.deltaD;
    stats.deltaB(ii,:)   = datai.stats.deltaB;
    stats.sparsity(ii,:) = datai.stats.sparsity;
    stats.time(ii,:)     = datai.stats.time;
    err.NRMSE(ii,:)  = datai.err.NRMSE;
    err.pSNR(ii,:)   = datai.err.pSNR;
    err.pNRMSE(ii,:) = datai.err.pNRMSE;
    err.ppSNR(ii,:)  = datai.err.ppSNR;
    
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
