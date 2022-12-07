function onlineDls_agg(parFcn)
% Syntax: onlineDls_agg(parFcn);

% Get parameters
[vars, ~] = parFcn();

% Initialize storage
nP = numel(vars.p);
nS = numel(vars.seed);
nL = numel(vars.lambda);
nM = numel(vars.mu);
nR = numel(vars.dr);
nG = numel(vars.gamma);
nt = vars.idim(3);
ni = numel(1:vars.dt:(nt - vars.T + 1));
nIters = vars.nItersi + vars.nIters * (ni - 1);
stats.cost     = nan(nP,nS,nL,nM,nR,nG,nIters);
stats.nrmse    = nan(nP,nS,nL,nM,nR,nG,nIters);
stats.deltaX   = nan(nP,nS,nL,nM,nR,nG,nIters);
stats.deltaD   = nan(nP,nS,nL,nM,nR,nG,nIters);
stats.deltaB   = nan(nP,nS,nL,nM,nR,nG,nIters);
stats.sparsity = nan(nP,nS,nL,nM,nR,nG,nIters);
stats.time     = nan(nP,nS,nL,nM,nR,nG,nIters);
err.NRMSE  = nan(nP,nS,nL,nM,nR,nG,1);
err.pSNR   = nan(nP,nS,nL,nM,nR,nG,1);
err.pNRMSE = nan(nP,nS,nL,nM,nR,nG,nt);
err.ppSNR  = nan(nP,nS,nL,nM,nR,nG,nt);

% Load data
files  = findMatchingFiles(vars.rawpath);
nFiles = numel(files);
for i = 1:nFiles
    itimer = tic;
    
    % Extract index
    [~, namei, ~] = fileparts(files{i});
    idx = str2double(regexp(namei,'\d+','match'));
    [ii, jj, kk, ll, mm, nn] = ind2sub([nP, nS, nL, nM, nR, nG],idx);
    
    % Save data
    datai = load(files{i},'stats','err');
    stats.cost(ii,jj,kk,ll,mm,nn,:)     = datai.stats.cost;
    stats.nrmse(ii,jj,kk,ll,mm,nn,:)    = datai.stats.nrmse;
    stats.deltaX(ii,jj,kk,ll,mm,nn,:)   = datai.stats.deltaX;
    stats.deltaD(ii,jj,kk,ll,mm,nn,:)   = datai.stats.deltaD;
    stats.deltaB(ii,jj,kk,ll,mm,nn,:)   = datai.stats.deltaB;
    stats.sparsity(ii,jj,kk,ll,mm,nn,:) = datai.stats.sparsity;
    stats.time(ii,jj,kk,ll,mm,nn,:)     = datai.stats.time;
    err.NRMSE(ii,jj,kk,ll,mm,nn,:)  = datai.err.NRMSE;
    err.pSNR(ii,jj,kk,ll,mm,nn,:)   = datai.err.pSNR;
    err.pNRMSE(ii,jj,kk,ll,mm,nn,:) = datai.err.pNRMSE;
    err.ppSNR(ii,jj,kk,ll,mm,nn,:)  = datai.err.ppSNR;
    
    % Display progress
    fprintf('File %i/%i processed [Time = %.2fs]\n',i,nFiles,toc(itimer));
end

%{
% Load min-NRMSE reconstruction
[~, im] = nanmin(err.NRMSE(:));
datam   = load(files{im},'Xhat','D');
Xhat    = datam.Xhat; %#ok
D       = datam.D; %#ok

% Save results
save(vars.outpath,'vars','stats','err','Xhat','D');
%}
save(vars.outpath,'vars','stats','err');
fprintf('DONE\n');
