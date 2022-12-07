function [vars, opts, vars0, opts0] = otazo_onlineUnitaryDls_par5b()
% Syntax:   [vars, opts, vars0, opts0] = otazo_onlineUnitaryDls_par5b();

% Knobs
vars.inpath   = 'data/otazo_full.mat';
vars.rawpath  = 'otazo_onlineUnitaryDls_par5b/data.mat';
vars.outpath  = 'otazo_onlineUnitaryDls_par5b.mat';
vars.T        = 5;
vars.dt       = 1;
vars.nItersi  = 75;
vars.nIters   = 15;
vars.nIters0  = 0;
vars.np       = 0;
vars.SNR      = inf;
vars.nReps    = 1;
% sweep
vars.ps       = 1 ./ [8]; % 1 ./ [4, 8, 12, 16, 20, 24];
vars.seed     = [1];
vars.lambda   = 0.035 * logspace(-1,1,5);
vars.mu       = 0.040 * logspace(-1,1,5);
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
opts.pgap     = [1, 1, 1]; % MORE PATCHES
opts.type     = 'hard';
%opts.dr
opts.ddim     = [prod(opts.pdim(1:(end - 1))), opts.pdim(end)]; %% N/A %%%%
opts.unitaryD = true; %% UNITARY-DIL %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
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
