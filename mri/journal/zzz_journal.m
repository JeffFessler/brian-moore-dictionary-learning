%% Load Otazo's perfusion data

clear;
rng(1);

% Load data
data = load('../cardiac_perf_R8.mat');
[ny nx nt nc] = size(data.kdata);
nd = ny * nx;
A = Emat_xyt(data.kdata(:,:,:,1) ~= 0,data.b1); % Aquisition operator
T = TempFFT(3); % Sparsifying transformation : 3D temporal DFT
Y = data.kdata; % Sub-sampled k-space data
truth = struct('Xtrue',[],'Ltrue',[],'Strue',[]);

% Clear extra variables
clearvars -except A T Y nx ny nd nt nc truth

%% Load PINCAT phantom data

clear;
rng(1);

% Knobs
SNR = 50; % SNR, in dB

% SNR to variance conversion function
SNR2sigma = @(SNR,Ytrue) exp(-SNR / 20) * (norm(Ytrue(:)) / sqrt(numel(Ytrue))) / sqrt(2);

% Load true data (and scale to [0 1])
data = load('../ktslr/aperiodic_pincat.mat');
Xtrue = data.new;
[ny nx nt] = size(Xtrue);
nd = ny * nx;
nc = 12; % From Otazo data
Xtrue = Xtrue - min(Xtrue(:));
Xtrue = Xtrue / max(Xtrue(:));
truth = struct('Xtrue',Xtrue,'Ltrue',[],'Strue',[]);

% Load system matrix
A = GenerateSystemMatrix(ny,nx,nt);

% Sparsifying transformation
T = TempFFT(3); % 3D temporal DFT

% Generate observations
Y = A * Xtrue;
M = (Y ~= 0);
nnzM = nnz(M);
sigma = SNR2sigma(SNR,Y(M));
if isreal(Xtrue)
    % Real noise
    Y(M) = Y(M) + sigma * randn(nnzM,1);
else
    % Complex noise
    Y(M) = Y(M) + (sigma / sqrt(2)) * (randn(nnzM,1) + 1i * randn(nnzM,1));
end

% Clear extra variables
clearvars -except A T Y nx ny nd nt nc truth sigma

%% Load in vivo perfusion data

clear;
rng(1);

% Knobs
SNR = 50; % SNR, in dB

% SNR to variance conversion function
SNR2sigma = @(SNR,Ytrue) exp(-SNR / 20) * (norm(Ytrue(:)) / sqrt(numel(Ytrue))) / sqrt(2);

% Load true data (and scale to [0 1])
data = load('../ktslr/invivo_perfusion.mat');
Xtrue = data.x;
[ny nx nt] = size(Xtrue);
nd = ny * nx;
nc = 12; % From Otazo data
Xtrue = Xtrue - min(Xtrue(:));
Xtrue = Xtrue / max(Xtrue(:));
truth = struct('Xtrue',Xtrue,'Ltrue',[],'Strue',[]);

% Load system matrix
A = GenerateSystemMatrix(ny,nx,nt);

% Sparsifying transformation
T = TempFFT(3); % 3D temporal DFT

% Generate observations
Y = A * Xtrue;
M = (Y ~= 0);
nnzM = nnz(M);
sigma = SNR2sigma(SNR,Y(M));
if isreal(Xtrue)
    % Real noise
    Y(M) = Y(M) + sigma * randn(nnzM,1);
else
    % Complex noise
    Y(M) = Y(M) + (sigma / sqrt(2)) * (randn(nnzM,1) + 1i * randn(nnzM,1));
end

% Clear extra variables
clearvars -except A T Y nx ny nd nt nc truth sigma

%% Run SVT algos

% Knobs
outpath = 'C:\Users\brimoor\Documents\Brian''s Stuff\OptShrink4MRI\LpS_SVT.mat';
%outpath = 'C:\Users\brimoor\Documents\Brian''s Stuff\OptShrink4MRI\LpS2_SVT.mat';
%outpath = 'C:\Users\brimoor\Documents\Brian''s Stuff\OptShrink4MRI\LaS_SVT.mat';
%outpath = 'C:\Users\brimoor\Documents\Brian''s Stuff\OptShrink4MRI\LaS2_SVT.mat';
algoStr = 'L+S PGM'; % {'L+S PGM','L+S ADMM','L&S ALM','L&S ADMM'}
%algoStr = 'L+S ADMM'; % {'L+S PGM','L+S ADMM','L&S ALM','L&S ADMM'}
%algoStr = 'L&S ALM'; % {'L+S PGM','L+S ADMM','L&S ALM','L&S ADMM'}
%algoStr = 'L&S ADMM'; % {'L+S PGM','L+S ADMM','L&S ALM','L&S ADMM'}
C = max(svd(reshape(A' * Y,[nd nt])));
minLambdaL = 0.001;
maxLambdaL = 0.001;
minLambdaS = 0.01;
maxLambdaS = 0.01;
NlambdaL = 1;
NlambdaS = 1;
%nb = [8 8]; % {8,24} for Fessler phantom

% Compute spectral norm of system matrix (normA = 20.275 when nt = 40)
%AAtfcn = @(x) reshape(A' * (A * reshape(x,[ny nx nt])),[],1);
%opts = struct('issym',true,'isreal',false,'disp',1); % 'maxit',5
%normA = sqrt(abs(eigs(AAtfcn,ny * nx * nt,1,'LM',opts)));
normA = 20.275;

% Compute parameter scaling
alpha = abs(normA)^2;
%{
if exist('truth','var') && isfield(truth,'Xtrue') && ~isempty(truth.Xtrue)
    % "Makes" max(abs(Xtrue(:))) == 1
    alpha = alpha / max(abs(truth.Xtrue(:)));
end
if exist('sigma','var')
    % Standardize for noise variance
    alpha = alpha * sigma^2;
end
%}

% Initializations
lambdaL = alpha * logspace(log10(minLambdaL),log10(maxLambdaL),NlambdaL);
lambdaS = alpha * logspace(log10(minLambdaS),log10(maxLambdaS),NlambdaS);
NRMSE_SVT = repmat(struct('X',[],'L',[],'S',[]),[NlambdaL NlambdaS]);
recon_SVT = repmat(struct('X',[],'L',[],'S',[]),[NlambdaL NlambdaS]);

% Update output path
if exist('nb','var')
    % Append locally low-rank (LLR) block size
    outpath = sprintf('%s_nb_%i%s',outpath(1:(end - 4)),nb,outpath((end - 3):end));
end

% Perform reconstructions
for j = 1:NlambdaS
    % SVT formulations
    for i = 1:NlambdaL
        timer = tic;
        method = struct();
        method.normA = normA;
        method.str = 'SVT';
        method.lambdaL = C * lambdaL(i);
        method.lambdaS = lambdaS(j);
        if exist('nb','var')
            % Do locally low-rank (LLR)
            method.nb = nb;
        end
        [NRMSE recon] = RunLSAlgo(algoStr,Y,A,T,nx,ny,nt,method,truth);
        NRMSE_SVT(i,j) = NRMSE;
        recon_SVT(i,j) = recon;
        fprintf('*** SVT (%i,%i)/(%i,%i) complete - Time = %.2fs ***\n',j,i,NlambdaS,NlambdaL,toc(timer));
    end
end

% Save data
tic;
save(outpath,'NRMSE_SVT','recon_SVT','lambdaL','lambdaS','algoStr','nx','ny','nt','-v7.3');
D = dir(outpath);
fprintf('Data saved - Size = %.2fGB - Time = %.2fs\n',D.bytes / 1e9,toc);

%% Run OptShrink algos

% Knobs
outpath = 'C:\Users\brimoor\Documents\Brian''s Stuff\OptShrink4MRI\LpS_OptShrink.mat';
%outpath = 'C:\Users\brimoor\Documents\Brian''s Stuff\OptShrink4MRI\LpS2_OptShrink.mat';
%outpath = 'C:\Users\brimoor\Documents\Brian''s Stuff\OptShrink4MRI\LaS_OptShrink.mat';
%outpath = 'C:\Users\brimoor\Documents\Brian''s Stuff\OptShrink4MRI\LaS2_OptShrink.mat';
algoStr = 'L+S PGM'; % {'L+S PGM','L+S ADMM','L&S ALM','L&S ADMM'}
%algoStr = 'L+S ADMM'; % {'L+S PGM','L+S ADMM','L&S ALM','L&S ADMM'}
%algoStr = 'L&S ALM'; % {'L+S PGM','L+S ADMM','L&S ALM','L&S ADMM'}
%algoStr = 'L&S ADMM'; % {'L+S PGM','L+S ADMM','L&S ALM','L&S ADMM'}
minLambdaS = 1e-4;
maxLambdaS = 10;
NlambdaS = 20;
Nr = 5;
%nb = [8 8]; % {8,24} for Fessler phantom

% Compute spectral norm of system matrix (normA = 20.275 when nt = 40)
%AAtfcn = @(x) reshape(A' * (A * reshape(x,[ny nx nt])),[],1);
%opts = struct('issym',true,'isreal',false,'disp',1); % 'maxit',5
%normA = sqrt(abs(eigs(AAtfcn,ny * nx * nt,1,'LM',opts)));
normA = 20.275;

% Compute parameter scaling
alpha = abs(normA)^2;
if exist('truth','var') && isfield(truth,'Xtrue') && ~isempty(truth.Xtrue)
    % "Makes" max(abs(Xtrue(:))) == 1
    alpha = alpha / max(abs(truth.Xtrue(:)));
end
if exist('sigma','var')
    % Standardize for noise variance
    alpha = alpha * sigma^2;
end

% Initializations
r = 1:Nr;
lambdaS = alpha * logspace(log10(minLambdaS),log10(maxLambdaS),NlambdaS);
NRMSE_OptShrink = repmat(struct('X',[],'L',[],'S',[]),[Nr NlambdaS]);
recon_OptShrink = repmat(struct('X',[],'L',[],'S',[]),[Nr NlambdaS]);

% Update output path
if exist('nb','var')
    % Append locally low-rank (LLR) block size
    outpath = sprintf('%s_nb_%i%s',outpath(1:(end - 4)),nb,outpath((end - 3):end));
end

% Perform reconstructions
for j = 1:NlambdaS
    % OptShrink formulations
    for i = 1:Nr
        timer = tic;
        method = struct();
        method.normA = normA;
        method.str = 'OptShrink';
        method.r = r(i);
        if exist('nb','var')
            % Do locally low-rank (LLR)
            method.nb = nb;
        end
        method.lambdaS = lambdaS(j);
        [NRMSE recon] = RunLSAlgo(algoStr,Y,A,T,nx,ny,nt,method,truth);
        NRMSE_OptShrink(i,j) = NRMSE;
        recon_OptShrink(i,j) = recon;
        fprintf('*** OptShrink (%i,%i)/(%i,%i) complete - Time = %.2fs ***\n',j,i,NlambdaS,Nr,toc(timer));
    end
end

% Save data
tic;
save(outpath,'NRMSE_OptShrink','recon_OptShrink','lambdaS','r','algoStr','nx','ny','nt','-v7.3');
D = dir(outpath);
fprintf('Data saved - Size = %.2fGB - Time = %.2fs\n',D.bytes / 1e9,toc);
