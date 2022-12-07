function [vars, opts, vars0, opts0] = otazo_onlineUnitaryDls_par3a()
% Syntax:   [vars, opts, vars0, opts0] = otazo_onlineUnitaryDls_par3a();

% Knobs
vars.inpath   = 'data/otazo_full.mat';
vars.rawpath  = 'otazo_onlineUnitaryDls_par3a/data.mat';
vars.outpath  = 'otazo_onlineUnitaryDls_par3a.mat';
vars.T        = 5;
vars.dt       = 1;
vars.nItersi  = 100;
vars.nIters   = 20;
vars.np       = 0;
vars.SNR      = inf;
vars.nReps    = 2;
% sweep
vars.ps       = 1 ./ [8]; % 1 ./ [4, 8, 12, 16, 20, 24];
vars.seed     = [1];
vars.lambda   = 0.045 * logspace(-0.5,0.5,5);
vars.mu       = 0.040 * logspace(-0.5,0.5,5);
% mu2 < 0: B01 = abs(mu2)% sparsity, B0t = Bt
% mu2 > 0: B01 = abs(mu2)% sparsity, B0t = abs(mu2)% sparsity
vars.mu2      = zeros(size(vars.mu)); % coupled
vars.dr       = [nan]; %% N/A %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
vars.gamma    = [0.9];
vars.gamma2   = [0.9]; % coupled
% Xmode = 0: use previous recons for new frames
% Xmode = 1: synthesized new frame recons
vars.Xmode    = 1;
% Dmode = -1: use previous dictionary (including t = 0 when nReps > 1)
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
opts.ddim     = [prod(opts.pdim(1:(end - 1))), opts.pdim(end)]; %% N/A %%%%
opts.unitaryD = true; %% UNITARY-DIL %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
opts.fixedD   = false;
%opts.nIters
% nIters0 > 0: hold D fixed for nIters0 outer iterations
opts.nIters0  = 0;
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
