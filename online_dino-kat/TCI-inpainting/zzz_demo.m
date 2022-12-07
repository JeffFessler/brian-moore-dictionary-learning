%%
% Online DINO-KAT demo
% coastguard dataset, noisy
%
% Brian Moore
% brimoor@umich.edu
%
% September 9, 2018
%

rng(1);

% Corruption
SNR = 19;
p = 0.7;

% Load data
data = load('data/coastguard2.mat');
Xtrue = im2double(data.M(:,:,1:150));
[ny, nx, nt] = size(Xtrue);

% Add noise
sigma = 10^(-SNR / 20) * norm(Xtrue(:)) / sqrt(numel(Xtrue));
Y = Xtrue + sigma * randn(ny,nx,nt);

% Missing data
M = (rand(ny,nx,nt) > p);
Y(~M) = 0;

% Cubic interpolation
X0 = interpVideo(Y,M,'grid2','cubic');


%% Online [r = 5]

T = 5;
dt = 1;

lambda = 0.05;
mu     = 0.15;
gamma  = 0.9;

opts.A        = 1;
opts.M        = M;
opts.xdim     = [ny, nx, T];
opts.pdim     = [8, 8, 5];
opts.pgap     = [2, 2, 1];
opts.type     = 'hard';
opts.dr       = 5;
opts.ddim     = [prod(opts.pdim(1:(end - 1))), opts.pdim(end)];
opts.fixedD   = false;
opts.unitaryD = false;
opts.nItersi  = 50;
opts.nIters   = 10;
opts.nItersDB = 1;
opts.nItersX  = -1;
opts.X0       = X0;
opts.D0       = dctmtx(prod(opts.pdim))';
opts.Xtrue    = Xtrue;
opts.flag     = 1;

%% Online [unitary]

lambda   = 0.05;
mu       = 0.075;
dr       = nan;
gamma    = [0.9];
gamma2   = [0.9];

opts.A        = 1;
%opts.M
%opts.xdim
opts.pdim     = [8, 8, 5];
opts.pgap     = [2, 2, 1];
opts.type     = 'hard';
%opts.dr
opts.ddim     = [prod(opts.pdim(1:(end - 1))), opts.pdim(end)];
opts.fixedD   = false;
opts.unitaryD = true;
%opts.nIters
opts.nItersDB = 1;
opts.nItersX  = -1;
%opts.X0
opts.D0       = dctmtx(prod(opts.pdim))';
%opts.B0
%opts.params
%opts.Xtrue
opts.flag     = 1;

%% Batch [r = 5]

vars.nIters   = 20;
vars.lambda   = 0.05;
vars.mu       = 0.15;
vars.dr       = 5;

opts.A        = 1;
%opts.M
%opts.xdim
opts.pdim     = [8, 8, 5];
opts.pgap     = [4, 4, 1];
opts.type     = 'hard';
%opts.dr
opts.ddim     = [prod(opts.pdim(1:(end - 1))), opts.pdim(end)];
opts.fixedD   = false;
%opts.nIters
opts.nItersDB = 1;
opts.nItersX  = -1;
opts.D0       = dctmtx(prod(opts.pdim))';
%opts.X0
%opts.Xtrue
opts.flag     = 1;
