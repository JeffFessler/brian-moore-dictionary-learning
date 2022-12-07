function [vars, opts, vars0, opts0] = pincat_onlineDls_par3b()
% Syntax:   [vars, opts, vars0, opts0] = pincat_onlineDls_par3b();

% Knobs
vars.inpath   = 'data/pincat_full.mat';
vars.rawpath  = 'pincat_onlineDls_par3b/data.mat';
vars.outpath  = 'pincat_onlineDls_par3b.mat';
vars.T        = 5;
vars.dt       = 1;
vars.nItersi  = 50;
vars.nIters   = 10;
vars.np       = 0;
vars.SNR      = 45; % add noise
vars.nReps    = 1; % optimized
% sweep
vars.nLines   = [5, 10, 15, 20, 25, 30]; % sweep
vars.seed     = [1];
vars.lambda   = [0.20]; % optimized
vars.mu       = [0.08]; % optimized
% mu2 < 0: B01 = abs(mu2)% sparsity, B0t = Bt
% mu2 > 0: B01 = abs(mu2)% sparsity, B0t = abs(mu2)% sparsity
vars.mu2      = [-10]; % coupled
vars.dr       = [1];
vars.gamma    = [0.9]; % optimized
vars.gamma2   = [0.9]; % coupled

% onlineDls opts
%opts.A;
opts.M        = 1;
%opts.xdim
opts.pdim     = [8, 8, 5];
opts.pgap     = [2, 2, 1];
opts.type     = 'hard';
%opts.dr
opts.ddim     = [prod(opts.pdim(1:(end - 1))), opts.pdim(end)];
opts.fixedD   = false;
%opts.nIters
opts.nItersDB = 1;
opts.nItersX  = 10;
%opts.X0
opts.D0       = dctmtx(prod(opts.pdim))';
%opts.B0
%opts.params
%opts.Xtrue
opts.accel    = true;
opts.tau      = 1;
opts.flag     = 1;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Baseline initialization
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
vars0 = struct();
opts0 = struct();
