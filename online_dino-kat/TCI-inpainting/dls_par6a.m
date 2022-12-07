function [vars, opts, vars0] = dls_par6a()
% Syntax:   [vars, opts, vars0] = dls_par6a();

% Knobs
vars.inpath   = 'data/gbus.mat';
vars.rawpath  = 'dls_par6a/data.mat';
vars.outpath  = 'dls_par6a.mat';
vars.idim     = [nan, nan, 150]; % 89,150
vars.nIters   = 20;
vars.np       = 0;
vars.SNR      = inf;
% sweep
vars.p        = [0.5, 0.6, 0.7, 0.8, 0.9];
vars.seed     = [1];
vars.lambda   = 0.0005;
vars.mu       = 0.34;
vars.dr       = 1;

% dls opts
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

% interpVideo opts
vars0.algo     = 'grid2';
vars0.method   = 'cubic';
