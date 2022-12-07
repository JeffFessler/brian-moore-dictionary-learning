function aggLpSsweep_opt()
% Syntax: aggLpSsweep_opt();

% Knobs
inpath  = 'otazo_R8_lps_opt/data.mat';
outpath = 'otazo_R8_lps_opt.mat';
r       = 1:10;
lambdaS = 0.01 * logspace(-1,1,25);
nIters  = 250;

% Initialize storage
nR = numel(r);
nS = numel(lambdaS);
time  = nan(nR,nS,nIters);
cost  = nan(nR,nS,nIters);
mse   = nan(nR,nS,nIters);
delta = nan(nR,nS,nIters);

% Load data
files  = findMatchingFiles(inpath);
nFiles = numel(files);
for i = 1:nFiles
    % Extract index
    [~,namei,~] = fileparts(files{i});
    idx         = str2double(regexp(namei,'\d+','match'));
    [ii]  = ind2sub([nS],idx);
    
    % Save data
    datai = load(files{i},'time','cost','mse','delta');
    time(:,ii,:)     = datai.time;
    cost(:,ii,:)     = datai.cost;
    mse(:,ii,:)      = datai.mse;
    delta(:,ii,:)    = datai.delta;
end

% Load optimal reconstruction
mseFinal   = mse(:,:,nIters);
[~,idxopt] = nanmin(nanmin(mseFinal,[],1));
[pathopt,fileopt,extopt] = fileparts(inpath);
inpathopt = sprintf('%s/%s%d%s',pathopt,fileopt,idxopt,extopt);
dataopt   = load(inpathopt);
Lhat = dataopt.Lhat; %#ok
Shat = dataopt.Shat; %#ok

% Save results
save(outpath,'Lhat','Shat','time','cost','mse','delta', ...
             'r','lambdaS','nIters');
fprintf('DONE\n');
