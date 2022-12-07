function pincat_LASSI_R8_agg(parFcn)
% Syntax: pincat_LASSI_R8_agg(parFcn);

% Get parameters
[vars, ~] = parFcn();

% Initialize storage
nR     = numel(vars.r);
nL     = numel(vars.lambdaL);
nS     = numel(vars.lambdaS);
nB     = numel(vars.lambdaB);
nD     = numel(vars.dr);
nIters = vars.nIters;
stats.cost     = nan(nR,nL,nS,nB,nD,nIters);
stats.nrmse    = nan(nR,nL,nS,nB,nD,nIters);
stats.delta    = nan(nR,nL,nS,nB,nD,nIters);
stats.sparsity = nan(nR,nL,nS,nB,nD,nIters);
stats.time     = nan(nR,nL,nS,nB,nD,nIters);

% Load data
files  = findMatchingFiles(vars.rawpath);
nFiles = numel(files);
for i = 1:nFiles
    itimer = tic;
    
    % Extract index
    [~,namei,~] = fileparts(files{i});
    idx         = str2double(regexp(namei,'\d+','match'));
    [ii, jj, kk, ll, mm] = ind2sub([nR, nL, nS, nB, nD],idx);
    
    % Save data
    datai = load(files{i},'stats');
    stats.cost(ii,jj,kk,ll,mm,:)     = datai.stats.cost;
    stats.nrmse(ii,jj,kk,ll,mm,:)    = datai.stats.nrmse;
    stats.delta(ii,jj,kk,ll,mm,:)    = datai.stats.delta;
    stats.sparsity(ii,jj,kk,ll,mm,:) = datai.stats.sparsity;
    stats.time(ii,jj,kk,ll,mm,:)     = datai.stats.time;
    
    % Display progress
    fprintf('File %i/%i processed [Time = %.2fs]\n',i,nFiles,toc(itimer));
end

% Load min-nrmse reconstruction
nrmsef  = stats.nrmse(:,:,:,:,:,end);
[~, im] = min(nrmsef(:));
datam   = load(files{im},'Lhat','Shat');
Lhat    = datam.Lhat; %#ok
Shat    = datam.Shat; %#ok

% Save results
save(vars.outpath,'vars','stats','Lhat','Shat');
fprintf('DONE\n');
