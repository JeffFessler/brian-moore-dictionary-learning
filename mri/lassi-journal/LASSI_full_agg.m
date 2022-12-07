function LASSI_full_agg(parFcn)
% Syntax: LASSI_full_agg(parFcn);

% Get parameters
[vars, ~] = parFcn();

% Initialize storage
nP     = numel(vars.ps);
nS     = numel(vars.seed);
nIters = vars.nIters;
stats.cost     = nan(nP,nS,nIters);
stats.nrmse    = nan(nP,nS,nIters);
stats.delta    = nan(nP,nS,nIters);
stats.sparsity = nan(nP,nS,nIters);
stats.time     = nan(nP,nS,nIters);

% Load data
files  = findMatchingFiles(vars.rawpath);
nFiles = numel(files);
for i = 1:nFiles
    itimer = tic;
    
    % Extract index
    [~,namei,~] = fileparts(files{i});
    idx         = str2double(regexp(namei,'\d+','match'));
    [ii, jj] = ind2sub([nP, nS],idx);
    
    % Save data
    datai = load(files{i},'stats');
    stats.cost(ii,jj,:)     = datai.stats.cost;
    stats.nrmse(ii,jj,:)    = datai.stats.nrmse;
    stats.delta(ii,jj,:)    = datai.stats.delta;
    stats.sparsity(ii,jj,:) = datai.stats.sparsity;
    stats.time(ii,jj,:)     = datai.stats.time;
    
    % Display progress
    fprintf('File %i/%i processed [Time = %.2fs]\n',i,nFiles,toc(itimer));
end

% Save results
save(vars.outpath,'vars','stats');
fprintf('DONE\n');
