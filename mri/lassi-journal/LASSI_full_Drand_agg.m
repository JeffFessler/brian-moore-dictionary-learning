function LASSI_full_Drand_agg(parFcn)
% Syntax: LASSI_full_Drand_agg(parFcn);

% Get parameters
[vars, ~] = parFcn();

% Initialize storage
nP     = numel(vars.ps);
nS     = numel(vars.seed);
nD     = numel(vars.D0seed);
nIters = vars.nIters;
stats.cost     = nan(nP,nS,nD,nIters);
stats.nrmse    = nan(nP,nS,nD,nIters);
stats.delta    = nan(nP,nS,nD,nIters);
stats.sparsity = nan(nP,nS,nD,nIters);
stats.time     = nan(nP,nS,nD,nIters);

% Load data
files  = findMatchingFiles(vars.rawpath);
nFiles = numel(files);
for i = 1:nFiles
    itimer = tic;
    
    % Extract index
    [~,namei,~] = fileparts(files{i});
    idx         = str2double(regexp(namei,'\d+','match'));
    [ii, jj, kk] = ind2sub([nP, nS, nD],idx);
    
    % Save data
    datai = load(files{i},'stats');
    stats.cost(ii,jj,kk,:)     = datai.stats.cost;
    stats.nrmse(ii,jj,kk,:)    = datai.stats.nrmse;
    stats.delta(ii,jj,kk,:)    = datai.stats.delta;
    stats.sparsity(ii,jj,kk,:) = datai.stats.sparsity;
    stats.time(ii,jj,kk,:)     = datai.stats.time;
    
    % Display progress
    fprintf('File %i/%i processed [Time = %.2fs]\n',i,nFiles,toc(itimer));
end

% Save results
save(vars.outpath,'vars','stats');
fprintf('DONE\n');
