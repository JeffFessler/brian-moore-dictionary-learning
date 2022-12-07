function aggLADSfull_p1_250()
% Syntax: aggLADSfull_p1_250();

% Knobs
% L + S init: 250 iteraions w/ lL = 1.1955 and lS = 0.01
inpath    = 'otazo_full_lads_p1_250/data.mat';
outpath   = 'otazo_full_lads_p1_250.mat';
ps        = [1/20, 1/16, 1/12, 1/8, 1/4];
seed      = 1;
lambdaL   = 0.5;
lambdaS   = 0.01;
lambdaB   = 0.03;
lambdaB0  = 0.03;
p         = 1;      % !!!!!!!!!!!!!!!!!!
nIters    = 50;

% Initialize storage
nP = numel(ps);
nS = numel(seed);
time     = nan(nP,nIters,nS);
cost     = nan(nP,nIters,nS);
mse      = nan(nP,nIters,nS);
delta    = nan(nP,nIters,nS);
sparsity = nan(nP,nIters,nS);

% Load data
files  = findMatchingFiles(inpath);
nFiles = numel(files);
for i = 1:nFiles
    % Extract index
    [~,namei,~] = fileparts(files{i});
    idx         = str2double(regexp(namei,'\d+','match'));
    [ii,jj]     = ind2sub([nP nS],idx);
    
    % Save data
    datai = load(files{i},'ite','time','cost','mse','delta','sparsity');
    time(ii,:,jj)     = datai.time;
    cost(ii,:,jj)     = datai.cost;
    mse(ii,:,jj)      = datai.mse;
    delta(ii,:,jj)    = datai.delta;
    sparsity(ii,:,jj) = datai.sparsity;
end

% Save results
save(outpath,'time','cost','mse','delta','sparsity','ps','seed','lambdaL','lambdaS','lambdaB','lambdaB0','p','nIters');
fprintf('DONE\n');
