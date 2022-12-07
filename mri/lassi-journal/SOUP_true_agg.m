function SOUP_true_agg(parFcn)
% Syntax: SOUP_true_agg(parFcn);

% Get parameters
[vars, ~] = parFcn();

% Initialize storage
nB     = numel(vars.lambdaB);
nD     = numel(vars.dr);
nIters = vars.nIters;
stats.cost     = nan(nB,nD,nIters);
stats.repError = nan(nB,nD,nIters);
stats.deltaD   = nan(nB,nD,nIters);
stats.deltaB   = nan(nB,nD,nIters);
stats.sparsity = nan(nB,nD,nIters);
stats.time     = nan(nB,nD,nIters);

% Load data
files  = findMatchingFiles(vars.rawpath);
nFiles = numel(files);
for i = 1:nFiles
    itimer = tic;
    
    % Extract index
    [~,namei,~] = fileparts(files{i});
    idx         = str2double(regexp(namei,'\d+','match'));
    [ii, jj]    = ind2sub([nB, nD],idx);
    
    % Save data
    datai = load(files{i},'stats');
    stats.cost(ii,jj,:)     = datai.stats.cost;
    stats.repError(ii,jj,:) = datai.stats.repError;
    stats.deltaD(ii,jj,:)   = datai.stats.deltaD;
    stats.deltaB(ii,jj,:)   = datai.stats.deltaB;
    stats.sparsity(ii,jj,:) = datai.stats.sparsity;
    stats.time(ii,jj,:)     = datai.stats.time;
    
    % Display progress
    fprintf('File %i/%i processed [Time = %.2fs]\n',i,nFiles,toc(itimer));
end

% Load min-repError reconstruction
repErrors = stats.repError(:,:,end);
[~, im]   = min(repErrors(:));
datam     = load(files{im},'Dhat');
Dhat      = datam.Dhat; %#ok

% Save results
save(vars.outpath,'vars','stats','Dhat');
fprintf('DONE\n');
