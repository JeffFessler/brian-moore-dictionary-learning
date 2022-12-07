function [vars, opts] = road_dls_par3a()
% Syntax:   [vars, opts] = road_dls_par3a();

% Knobs
vars.inpath   = 'data/road100.mat';
vars.rawpath  = 'road_dls_par3a/data.mat';
vars.outpath  = 'road_dls_par3a.mat';
vars.idim     = [128, nan, 100];
vars.nIters   = 20;
vars.np       = 0;
vars.algo     = 'paint2';
vars.method   = 1;
vars.SNR      = inf;
% sweep
vars.p        = [0.5, 0.6, 0.7, 0.8, 0.9];
vars.seed     = [1];
vars.lambda   = [1e-5];
vars.mu       = [0.15];
vars.dr       = [1];

% dls opts
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
opts.D0       = dctmtx(prod(opts.pdim))';
%opts.X0
%opts.Xtrue
opts.flag     = 1;
