function [vars, opts, vars0, opts0] = invivo_onlineDls_par5()
% Syntax:   [vars, opts, vars0, opts0] = invivo_onlineDls_par5();

% Knobs
vars.inpath   = 'data/invivo_full.mat';
vars.rawpath  = 'invivo_onlineDls_par5/data.mat';
vars.outpath  = 'invivo_onlineDls_par5.mat';
vars.T        = 5;
vars.dt       = 1;
vars.nItersi  = 50;
vars.nIters   = 10;
vars.np       = 0;
vars.SNR      = inf;
vars.nReps    = 1; % optimized
% sweep
vars.nLines   = [5, 10, 15, 20, 25, 30]; % sweep
vars.seed     = [1];
vars.lambda   = [0.001]; % optimized
vars.mu       = [0.01]; % optimized
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
% RPCA initialization
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

vars0.r        = nan;
vars0.lambdaL  = 0.158;
vars0.lambdaS  = 0.0025;
vars0.nIters   = 250;

% RPCA opts
%opts0.A      = A;
%opts0.T      = T;
%opts0.r      = r;
%opts0.nIters = vars.nIters;
%opts0.L0     = reshape(data.Xfft,[],nt);
%opts0.S0     = zeros(ny * nx,nt);
%opts0.Xtrue  = reshape(data.Xtrue,[],nt);
opts0.accel   = false;
opts0.tau     = 1;
opts0.flag    = 1;
