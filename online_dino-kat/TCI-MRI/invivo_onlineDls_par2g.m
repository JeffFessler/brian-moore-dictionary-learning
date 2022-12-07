function [vars, opts, vars0, opts0] = invivo_onlineDls_par2g()
% Syntax:   [vars, opts, vars0, opts0] = invivo_onlineDls_par2g();

% Knobs
vars.inpath   = 'data/invivo_full.mat';
vars.rawpath  = 'invivo_onlineDls_par2g/data.mat';
vars.outpath  = 'invivo_onlineDls_par2g.mat';
vars.T        = 5;
vars.dt       = 1;
vars.nItersi  = 100; % MORE
vars.nIters   = 20; % MORE
vars.np       = 0;
vars.SNR      = inf;
vars.nReps    = 1;
% sweep
vars.nLines   = [15]; % [5, 10, 15, 20, 25, 30];
vars.seed     = [1];
vars.lambda   = 0.01 * logspace(-0.5,0.5,5);
vars.mu       = 0.04 * logspace(-0.5,0.5,5);
% mu2 < 0: B01 = abs(mu2)% sparsity, B0t = Bt
% mu2 > 0: B01 = abs(mu2)% sparsity, B0t = abs(mu2)% sparsity
vars.mu2      = -10 * ones(size(vars.mu)); % coupled
vars.dr       = [1];
vars.gamma    = [0.9];
vars.gamma2   = [0.9]; % coupled
% Xmode = 0: use previous recons for new frames
% Xmode = 1: synthesized new frame recons
vars.Xmode    = 1;
% Dmode = 0: use previous dictionary at each minibatch
% Dmode = 1: use initial dictionary at each minibatch
vars.Dmode    = 0;

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
% nIters0 > 0: hold D fixed for nIters0 outer iterations
opts.nIters0  = 5;
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
