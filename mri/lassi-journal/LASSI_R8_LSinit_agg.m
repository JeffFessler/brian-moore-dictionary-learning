function LASSI_R8_LSinit_agg(parFcn)
% Syntax: LASSI_R8_LSinit_agg(parFcn);

% Get parameters
[vars, ~] = parFcn();

% Initialize storage
nI     = numel(vars.initpath);
nIters = vars.nIters;
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
    datai = load(files{i},'stats');
    stats.cost(ii,:)     = datai.stats.cost;
    stats.nrmse(ii,:)    = datai.stats.nrmse;
    stats.delta(ii,:)    = datai.stats.delta;
    stats.sparsity(ii,:) = datai.stats.sparsity;
    stats.time(ii,:)     = datai.stats.time;
    
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
