function invivo_LASSI_R8_Diter_agg(parFcn)
% Syntax: invivo_LASSI_R8_Diter_agg(parFcn);

% Get parameters
[vars, ~] = parFcn();

% Initialize storage
nI     = numel(vars.nIters);
nIters = max(vars.nIters);
stats.cost     = nan(nI,nIters);
stats.nrmse    = nan(nI,nIters);
stats.delta    = nan(nI,nIters);
stats.sparsity = nan(nI,nIters);
stats.time     = nan(nI,nIters);

% Load data
files  = findMatchingFiles(vars.rawpath);
nFiles = numel(files);
for i = 1:nFiles
    itimer = tic;
    
    % Extract index
    [~,namei,~] = fileparts(files{i});
    idx         = str2double(regexp(namei,'\d+','match'));
    [ii]        = ind2sub([nI],idx);
    
    % Save data
    gap   = nIters / vars.nIters(ii);
    datai = load(files{i},'stats');
    stats.cost(ii,gap:gap:end)     = datai.stats.cost;
    stats.nrmse(ii,gap:gap:end)    = datai.stats.nrmse;
    stats.delta(ii,gap:gap:end)    = datai.stats.delta;
    stats.sparsity(ii,gap:gap:end) = datai.stats.sparsity;
    stats.time(ii,gap:gap:end)     = datai.stats.time;
    
    % Display progress
    fprintf('File %i/%i processed [Time = %.2fs]\n',i,nFiles,toc(itimer));
end

% Load min-nrmse reconstruction
nrmsef  = stats.nrmse(:,end);
[~, im] = min(nrmsef(:));
datam   = load(files{im},'Lhat','Shat');
Lhat    = datam.Lhat; %#ok
Shat    = datam.Shat; %#ok

% Save results
save(vars.outpath,'vars','stats','Lhat','Shat');
fprintf('DONE\n');
