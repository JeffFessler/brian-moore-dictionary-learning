function [vars, opts, vars0] = onlineDls_npar3c()
% Syntax:   [vars, opts, vars0] = onlineDls_npar3c();

% Knobs
vars.inpath   = 'data/coastguard2.mat';
vars.rawpath  = 'onlineDls_npar3c/data.mat';
vars.outpath  = 'onlineDls_npar3c.mat';
vars.idim     = [nan, nan, 150]; % 89,150
vars.T        = 5;
vars.dt       = 1;
vars.nItersi  = 50;
vars.nIters   = 10;
vars.np       = 0;
vars.SNR      = 19; %% NOISY %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% sweep
vars.p        = [0.7];
vars.seed     = [1];
vars.lambda   = 0.01 * logspace(-0.7,0.7,5);
vars.mu       = 0.15 * logspace(-0.7,0.7,5);
vars.dr       = nan; %%%%% N/A %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
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
opts.fixedD   = true; %%%% FIXED DICTIONARY %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
opts.unitaryD = false;
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
