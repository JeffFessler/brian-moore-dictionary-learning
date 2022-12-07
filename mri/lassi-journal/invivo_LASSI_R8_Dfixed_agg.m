function invivo_LASSI_R8_Dfixed_agg(parFcn)
% Syntax: invivo_LASSI_R8_Dfixed_agg(parFcn);

% Get parameters
[vars, ~] = parFcn();

% Initialize storage
nD     = numel(vars.fixedD);
nIters = vars.nIters;
stats.cost     = nan(nD,nIters);
stats.nrmse    = nan(nD,nIters);
stats.delta    = nan(nD,nIters);
stats.sparsity = nan(nD,nIters);
stats.time     = nan(nD,nIters);

% Load data
files  = findMatchingFiles(vars.rawpath);
nFiles = numel(files);
for i = 1:nFiles
    itimer = tic;
    
    % Extract index
    [~,namei,~] = fileparts(files{i});
    idx         = str2double(regexp(namei,'\d+','match'));
    [ii]        = ind2sub([nD],idx);
    
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
