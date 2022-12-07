function [vars, opts] = road_onlineDls_par13c2()
% Syntax:   [vars, opts] = road_onlineDls_par13c2();

% Knobs
vars.inpath   = 'data/road100.mat';
vars.rawpath  = 'road_onlineDls_par13c2/data.mat';
vars.outpath  = 'road_onlineDls_par13c2.mat';
vars.idim     = [128, nan, 100];
vars.T        = 5;
vars.dt       = 1;
vars.nItersi  = 50;
vars.nIters   = 10;
vars.np       = 0;
vars.algo     = 'paint2';
vars.method   = 1;
vars.SNR      = 19; % yields pSNR = 25dB
% sweep
vars.p        = [0.5, 0.6, 0.7, 0.8, 0.9];
vars.seed     = [1];
vars.lambda   = [0.01];
vars.mu       = [0.13];
vars.dr       = [1, 3, 5];
vars.gamma    = [0.95];
vars.gamma2   = [0.95]; % coupled

% onlineDls opts
opts.A        = 1;
%opts.M
%opts.xdim
opts.pdim     = [8, 8, 5];
opts.pgap     = [1, 1, 1];
opts.type     = 'hard';
%opts.dr
opts.ddim     = [prod(opts.pdim(1:(end - 1))), opts.pdim(end)];
opts.fixedD   = false;
%opts.nIters
opts.nItersDB = 1;
opts.nItersX  = -1;
%opts.X0
opts.D0       = dctmtx(prod(opts.pdim))';
%opts.B0
%opts.params
%opts.Xtrue
opts.flag     = 1;
