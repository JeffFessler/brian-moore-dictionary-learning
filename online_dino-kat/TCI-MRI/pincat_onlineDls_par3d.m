function [vars, opts, vars0, opts0] = pincat_onlineDls_par3d()
% Syntax:   [vars, opts, vars0, opts0] = pincat_onlineDls_par3d();

% Knobs
vars.inpath   = 'data/pincat_full.mat';
vars.rawpath  = 'pincat_onlineDls_par3d/data.mat';
vars.outpath  = 'pincat_onlineDls_par3d.mat';
vars.T        = 5;
vars.dt       = 1;
%vars.nItersi  = 100; % MORE
%vars.nIters   = 20; % MORE
vars.nItersi  = 15; % RUNTIMES ONLY
vars.nIters   = 15; % RUNTIMES ONLY
vars.nIters0  = 5;
vars.np       = 0;
vars.SNR      = 45; % add noise
vars.nReps    = 2; % MORE
% sweep
vars.nLines   = [5, 10, 15, 20, 25, 30];
vars.seed     = [1, 2, 3];
vars.lambda   = 0.20;
vars.mu       = 0.08;
% mu2 < 0: B01 = abs(mu2)% sparsity, B0t = Bt
% mu2 > 0: B01 = abs(mu2)% sparsity, B0t = abs(mu2)% sparsity
vars.mu2      = 0;
vars.dr       = [1];
vars.gamma    = [0.9];
vars.gamma2   = [0.9]; % coupled
% Xmode = 0: use previous recons for new frames
% Xmode = 1: synthesized new frame recons
vars.Xmode    = 1;
% Dmode = -1: use previous dictionary (including t = 0 when nReps > 1)
% Dmode = 0: use previous dictionary at each minibatch
% Dmode = 1: use initial dictionary at each minibatch
vars.Dmode    = -1;
vars.resetParams = false;

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
opts.fixedD   = false;
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
