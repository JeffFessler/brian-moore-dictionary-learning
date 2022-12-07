function aggLADSp250b()
% Syntax: aggLADSp250b();

%{
% Knobs
% L + S init: 250 iteraions w/ lL = 1.1955 and lS = 0.01
inpath    = 'otazo_p_lads250/data.mat';
outpath   = 'otazo_p_lads250.mat';
lambdaL   = 0.5;
lambdaS   = 0.01;
lambdaB   = 0.03;
lambdaB0  = 0.03;
p         = [1, 2, 3, 4, 5];
nIters    = 50;
%}

% Knobs
% L + S init: 250 iteraions w/ lL = 0.525 and lS = 0.01
inpath    = 'otazo_p_lads250b/data.mat';
outpath   = 'otazo_p_lads250b.mat';
lambdaL   = 0.5;
lambdaS   = 0.01;
lambdaB   = 0.04;
lambdaB0  = 0.04;
p         = [1, 2, 3, 4, 5];
nIters    = 50;

%{
% Knobs
% L + S init: 20 iteraions w/ lL = 1.1955 and lS = 0.01
inpath    = 'otazo_p_lads20/data.mat';
outpath   = 'otazo_p_lads20.mat';
lambdaL   = 0.5;
lambdaS   = 0.01;
lambdaB   = 0.03;
lambdaB0  = 0.06;
p         = [1, 2, 3, 4, 5];
nIters    = 50;
%}

% Initialize storage
nP = numel(p);
time     = nan(nP,nIters);
cost     = nan(nP,nIters);
mse      = nan(nP,nIters);
delta    = nan(nP,nIters);
sparsity = nan(nP,nIters);

% Load data
files  = findMatchingFiles(inpath);
nFiles = numel(files);
for i = 1:nFiles
    % Extract index
    [~,namei,~] = fileparts(files{i});
    idx         = str2double(regexp(namei,'\d+','match'));
    [ii]        = ind2sub([nP],idx);
    
    % Save data
    datai = load(files{i},'ite','time','cost','mse','delta','sparsity');
    time(ii,:)     = datai.time;
    cost(ii,:)     = datai.cost;
    mse(ii,:)      = datai.mse;
    delta(ii,:)    = datai.delta;
    sparsity(ii,:) = datai.sparsity;
end

% Save results
save(outpath,'time','cost','mse','delta','sparsity','lambdaL','lambdaS','lambdaB','lambdaB0','p','nIters');
fprintf('DONE\n');
