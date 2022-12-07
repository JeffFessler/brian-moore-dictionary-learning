function [vars, opts, vars0, opts0] = pincat_onlineUnitaryDls_par6f()
% Syntax:   [vars, opts, vars0, opts0] = pincat_onlineUnitaryDls_par6f();

% Knobs
vars.inpath   = 'data/pincat_full.mat';
vars.rawpath  = 'pincat_onlineUnitaryDls_par6f/data.mat';
vars.outpath  = 'pincat_onlineUnitaryDls_par6f.mat';
vars.T        = 5;
vars.dt       = 1;
vars.nItersi  = [50, 50];
vars.nIters   = [10, 15];
vars.nIters0  = 0;
vars.np       = 0;
vars.SNR      = 45; % add noise
vars.nReps    = 2;
% sweep
vars.nLines   = [5, 10, 15, 20, 25, 30];
vars.seed     = [1, 2, 3];
vars.lambda   = 0.125;
vars.mu       = 0.05;
vars.mu2      = 0;
vars.dr       = [nan]; %% N/A %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
vars.gamma    = [0.9];
vars.gamma2   = [0.9]; % coupled
% Xmode = 0: use previous recons for new frames
% Xmode = 1: synthesized new frame recons
vars.Xmode    = 1;
% Dmode = -1: use previous dictionary (including t = 0 when nReps > 1)
% Dmode = 0: use previous dictionary at each minibatch
% Dmode = 1: use initial dictionary at each minibatch
vars.Dmode    = -1;
vars.resetParams = true; % reset params for 2nd rep

vars.fixedD   = [true, false]; % specify for each rep individually
vars.lambda0  = 0.35;          % override param for 1st rep
vars.mu0      = 0.045;         % override param for 1st rep
vars.pgap     = {[2, 2, 1], [1, 1, 1]}; % specify for each rep individually

% onlineDls opts
%opts.A;
opts.M        = 1;
%opts.xdim
opts.pdim     = [8, 8, 5];
%opts.pgap
opts.type     = 'hard';
%opts.dr
opts.ddim     = [prod(opts.pdim(1:(end - 1))), opts.pdim(end)]; %% N/A %%%%
opts.unitaryD = true; %% UNITARY-DIL %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%opts.fixedD
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
