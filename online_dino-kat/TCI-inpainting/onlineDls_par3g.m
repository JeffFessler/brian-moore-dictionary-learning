function [vars, opts, vars0] = onlineDls_par3g()
% Syntax:   [vars, opts, vars0] = onlineDls_par3g();

% Knobs
vars.inpath   = 'data/coastguard2.mat';
vars.rawpath  = 'onlineDls_par3g/data.mat';
vars.outpath  = 'onlineDls_par3g.mat';
vars.idim     = [nan, nan, 150]; % 89,150
vars.T        = 5;
vars.dt       = 1;
vars.nItersi  = 50;
vars.nIters   = 10;
vars.np       = 0;
vars.SNR      = inf;
% sweep
vars.p        = [0.7];
vars.seed     = [1];
vars.lambda   = 0.0005;
vars.mu       = 0.08 * logspace(-0.7,0.7,9);
vars.dr       = nan; %% N/A %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
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
opts.unitaryD = true; %% UNITARY %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%opts.nIters
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
