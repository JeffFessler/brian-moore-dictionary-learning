function [vars, opts, vars0] = onlineDls_par8e()
% Syntax:   [vars, opts, vars0] = onlineDls_par8e();

% Knobs
vars.inpath   = 'data/gflower.mat';
vars.rawpath  = 'onlineDls_par8e/data.mat';
vars.outpath  = 'onlineDls_par8e.mat';
vars.idim     = [nan, nan, 150]; % 89,150
vars.T        = 5;
vars.dt       = 1;
vars.nItersi  = 50;
vars.nIters   = 15;
vars.np       = 0;
vars.SNR      = inf;
% sweep
vars.p        = [0.5, 0.6, 0.7, 0.8, 0.9];
vars.seed     = [1];
vars.lambda   = 0.0005;
vars.mu       = 0.18;
vars.dr       = 5;
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
opts.unitaryD = false;
%opts.nIters
opts.nIters0  = 5; %% UPDATE B FIRST %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
opts.nItersDB = 1;
opts.nItersX  = -1;
%opts.X0
opts.D0       = dctmtx(prod(opts.pdim))';
%opts.B0
%opts.params
%opts.Xtrue
opts.flag     = 1;

% interpVideo opts
vars0.algo     = 'grid2';
vars0.method   = 'cubic';
