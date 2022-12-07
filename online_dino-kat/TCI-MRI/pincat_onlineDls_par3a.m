function [vars, opts, vars0, opts0] = pincat_onlineDls_par3a()
% Syntax:   [vars, opts, vars0, opts0] = pincat_onlineDls_par3a();

% Knobs
vars.inpath   = 'data/pincat_full.mat';
vars.rawpath  = 'pincat_onlineDls_par3a/data.mat';
vars.outpath  = 'pincat_onlineDls_par3a.mat';
vars.T        = 5;
vars.dt       = 1;
vars.nItersi  = 50;
vars.nIters   = 10;
vars.nIters0  = 0;
vars.np       = 0;
vars.SNR      = 45; % add noise
vars.nReps    = 1;
% sweep
vars.nLines   = [5, 10, 15, 20, 25, 30];
vars.seed     = [1, 2, 3];
vars.lambda   = 0.35;
vars.mu       = 0.045;
% mu2 < 0: B01 = abs(mu2)% sparsity, B0t = Bt
% mu2 > 0: B01 = abs(mu2)% sparsity, B0t = abs(mu2)% sparsity
vars.mu2      = 0;
vars.dr       = [nan];
vars.gamma    = [0.9];
vars.gamma2   = [0.9]; % coupled
% Xmode = 0: use previous recons for new frames
% Xmode = 1: synthesized new frame recons
vars.Xmode    = 1;
% Dmode = 0: use previous dictionary at each minibatch
% Dmode = 1: use initial dictionary at each minibatch
vars.Dmode    = -1;

% onlineDls opts
%opts.A;
opts.M        = 1;
%opts.xdim
opts.pdim     = [8, 8, 5];
opts.pgap     = [2, 2, 1];
opts.type     = 'hard';
%opts.dr
opts.ddim     = [prod(opts.pdim(1:(end - 1))), opts.pdim(end)];
opts.unitaryD = false;
opts.fixedD   = true;
%opts.nIters
%opts.nIters0
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
