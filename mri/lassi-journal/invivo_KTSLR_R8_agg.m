function invivo_KTSLR_R8_agg(parFcn)
% Syntax: invivo_KTSLR_R8_agg(parFcn);

% Get parameters
[vars, ~] = parFcn();

% Initialize storage
nP     = numel(vars.p);
n1     = numel(vars.mu1);
n2     = numel(vars.mu2);
nIters = vars.nItersO * vars.nItersI;
stats.cost  = nan(nP,n1,n2,nIters);
stats.nrmse = nan(nP,n1,n2,nIters);
stats.delta = nan(nP,n1,n2,nIters);
stats.time  = nan(nP,n1,n2,nIters);

% Load data
files  = findMatchingFiles(vars.rawpath);
nFiles = numel(files);
for i = 1:nFiles
    itimer = tic;
    
    % Extract index
    [~,namei,~] = fileparts(files{i});
    idx         = str2double(regexp(namei,'\d+','match'));
    [ii, jj, kk] = ind2sub([nP, n1, n2],idx);
    
    % Save data
    datai = load(files{i},'stats');
    stats.cost(ii,jj,kk,:)  = datai.stats.cost;
    stats.nrmse(ii,jj,kk,:) = datai.stats.nrmse;
    stats.delta(ii,jj,kk,:) = datai.stats.delta;
    stats.time(ii,jj,kk,:)  = datai.stats.time;
    
    % Display progress
    fprintf('File %i/%i processed [Time = %.2fs]\n',i,nFiles,toc(itimer));
end

% Load min-nrmse reconstruction
nrmsef  = stats.nrmse(:,:,:,end);
[~, im] = min(nrmsef(:));
datam   = load(files{im},'Xhat');
Xhat    = datam.Xhat; %#ok

% Save results
save(vars.outpath,'vars','stats','Xhat');
fprintf('DONE\n');
