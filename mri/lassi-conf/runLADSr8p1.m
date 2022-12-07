function runLADSr8p1(idx)
% Syntax: runLADSr8p1(idx);

% Knobs
inpath    = 'otazo_R8.mat';
initpath  = 'otazo_R8_lps_visual.mat';
outpath   = 'otazo_R8_lads_p1_mse250.mat';
lambdaL   = 0.5;
lambdaS   = 0.01;
lambdaB   = 0.03;
lambdaB0  = 0.03;
p         = 1;
nIters    = 50;

% Check index
if ~exist('idx','var') || isempty(idx)
    % Return job count
    fprintf('Total jobs: %d\n',1);
    return;
elseif idx ~= 1
    fprintf('Only idx = 1 supported\n');
    return;
end

% Set random seed
rng(1);

% Load undersampled data
load(inpath);
A = Emat_xyt(mask,samp);

% Load L + S initialization
init = load(initpath);
L0   = init.Lhat;
S0   = init.Shat;

% LADS
lB = logspace(log10(lambdaB0),log10(lambdaB),nIters);
[Lhat,Shat,~,~,ite,time,cost,mse,delta,sparsity] = runLADS(A,Y,Xtrue,L0,S0,lambdaL,lambdaS,lB,p,nIters); %#ok

% Save results
save(outpath,'Lhat','Shat','time','cost','mse','delta','sparsity', ...
             'lambdaL','lambdaS','lambdaB','lambdaB0','p','nIters');
fprintf('DONE\n');
