function [vars, opts] = road_onlineDls_par10b()
% Syntax:   [vars, opts] = road_onlineDls_par10b();

% Knobs
vars.inpath   = 'data/road100.mat';
vars.rawpath  = 'road_onlineDls_par10b/data.mat';
vars.outpath  = 'road_onlineDls_par10b.mat';
vars.idim     = [128, nan, 100];
vars.T        = 5;
vars.dt       = 1;
vars.nItersi  = 50;
vars.nIters   = 10;
vars.np       = 0;
vars.algo     = 'paint2';
vars.method   = 1;
vars.SNR      = inf;
% sweep
vars.p        = [0.7];
vars.seed     = [1];
vars.lambda   = [1e-6, 1e-5];
vars.mu       = [0.12, 0.13, 0.14, 0.15, 0.16, 0.17, 0.18];
vars.dr       = [1];
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
