function [vars, opts] = road_onlineDls_par1()
% Syntax:   [vars, opts] = road_onlineDls_par1();

% Knobs
vars.inpath   = 'data/road396.mat';
vars.rawpath  = 'road_onlineDls_par1/data.mat';
vars.outpath  = 'road_onlineDls_par1.mat';
vars.idim     = [128, nan, 200];
vars.T        = 5;
vars.dt       = 1;
vars.nItersi  = 20;
vars.nIters   = 10;
% sweep
vars.p        = [0.7];
vars.seed     = [1];
vars.lambda   = 0.01 * logspace(-1,1,5);
vars.mu       = 0.05 * logspace(-1,1,5);
vars.dr       = [1, 5];
vars.gamma    = [0.9];
vars.gamma2   = [0.9]; % coupled

% onlineDls opts
opts.A        = 1;
%opts.M
%opts.xdim
opts.pdim     = [8, 8, 5];
opts.pgap     = [2, 2, 1];
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
