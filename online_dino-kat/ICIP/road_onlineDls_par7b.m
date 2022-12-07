function [vars, opts] = road_onlineDls_par7b()
% Syntax:   [vars, opts] = road_onlineDls_par7b();

% Knobs
vars.inpath   = 'data/road396.mat';
vars.rawpath  = 'road_onlineDls_par7b/data.mat';
vars.outpath  = 'road_onlineDls_par7b.mat';
vars.idim     = [128, nan, 200];
vars.T        = 8;
vars.dt       = 1;
vars.nItersi  = 50;
vars.nIters   = 10;
vars.np       = 0;
vars.algo     = 'paint2';
vars.method   = 1;
% sweep
vars.p        = [0.5, 0.6, 0.7, 0.8, 0.9];
vars.seed     = [1];
vars.lambda   = [0.001];
vars.mu       = [0.08];
vars.dr       = [nan]; % FIXED_D ******************************************
vars.gamma    = [0.8];
vars.gamma2   = [0.8]; % coupled

% onlineDls opts
opts.A        = 1;
%opts.M
%opts.xdim
opts.pdim     = [8, 8, 8];
opts.pgap     = [1, 1, 1];
opts.type     = 'hard';
%opts.dr
opts.ddim     = [prod(opts.pdim(1:(end - 1))), opts.pdim(end)];
opts.fixedD   = true; % FIXED_D *******************************************
%opts.nIters
opts.nItersDB = 1;
opts.nItersX  = -1;
%opts.X0
opts.D0       = dctmtx(prod(opts.pdim))';
%opts.B0
%opts.params
%opts.Xtrue
opts.flag     = 1;
