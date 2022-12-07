function FFT_full_agg(parFcn)
% Syntax: FFT_full_agg(parFcn);

% Get parameters
vars = parFcn();

% Initialize storage
nP     = numel(vars.ps);
nS     = numel(vars.seed);
stats.nrmse = nan(nP,nS);
stats.time  = nan(nP,nS);

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
    stats.nrmse(ii,jj) = datai.stats.nrmse;
    stats.time(ii,jj)  = datai.stats.time;
    
    % Display progress
    fprintf('File %i/%i processed [Time = %.2fs]\n',i,nFiles,toc(itimer));
end

% Save results
save(vars.outpath,'vars','stats');
fprintf('DONE\n');
