function [vars, opts, vars0, opts0] = otazo_onlineUnitaryDls_par4d()
% Syntax:   [vars, opts, vars0, opts0] = otazo_onlineUnitaryDls_par4d();

% Knobs
vars.inpath   = 'data/otazo_full.mat';
vars.rawpath  = 'otazo_onlineUnitaryDls_par4d/data.mat';
vars.outpath  = 'otazo_onlineUnitaryDls_par4d.mat';
vars.T        = 5;
vars.dt       = 1;
vars.nItersi  = 50;
vars.nIters   = 15;
vars.nIters0  = 0;
vars.np       = 0;
vars.SNR      = inf;
vars.nReps    = 1;
% sweep
vars.ps       = 1 ./ [4, 8, 12, 16, 20, 24];
vars.seed     = [1, 2, 3];
vars.lambda   = 0.03;
vars.mu       = 0.04;
% mu2 < 0: B01 = abs(mu2)% sparsity, B0t = Bt
% mu2 > 0: B01 = abs(mu2)% sparsity, B0t = abs(mu2)% sparsity
vars.mu2      = 0;
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
% RPCA initialization
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Knobs
vars0.r       = nan;
vars0.lambdaL = 0.53;
vars0.lambdaS = 0.01;
vars0.nIters  = 250;

% RPCA opts
%opts0.A      = A;
%opts0.T      = T;
%opts0.r      = vars0.r;
%opts0.nIters = vars0.nIters;
%opts0.L0     = reshape(data.Xfft,[],nt);
%opts0.S0     = zeros(ny * nx,nt);
%opts0.Xtrue  = reshape(data.Xtrue,[],nt);
opts0.accel   = false;
opts0.tau     = 1;
opts0.flag    = 1;
