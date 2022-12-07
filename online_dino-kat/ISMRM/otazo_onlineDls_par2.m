function [vars, opts, vars0, opts0] = otazo_onlineDls_par2()
% Syntax:   [vars, opts, vars0, opts0] = otazo_onlineDls_par2();

% Knobs
vars.inpath   = 'data/otazo_full.mat';
vars.rawpath  = 'otazo_onlineDls_par2/data.mat';
vars.outpath  = 'otazo_onlineDls_par2.mat';
vars.T        = 5;
vars.dt       = 1;
vars.nItersi  = 50;
vars.nIters   = 10;
vars.np       = 0;
vars.SNR      = inf;
vars.nReps    = 3;
% sweep
vars.ps       = 1 / 8; % % 1 ./ [4, 8, 12, 16, 20, 24]
vars.seed     = [1];
vars.lambda   = [0.02];
vars.mu       = [0.10];
% mu2 < 0: B01 = abs(mu2)% sparsity, B0t = Bt
% mu2 > 0: B01 = abs(mu2)% sparsity, B0t = abs(mu2)% sparsity
vars.mu2      = -10 * ones(size(vars.mu)); % coupled
vars.dr       = [1];
vars.gamma    = linspace(0.7,1.0,7);
vars.gamma2   = linspace(0.7,1.0,7); % coupled

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