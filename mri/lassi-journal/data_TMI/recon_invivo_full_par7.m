function [vars, opts, vars0, opts0] = recon_invivo_full_par7()
% Syntax:   [vars, opts, vars0, opts0] = recon_invivo_full_par7();

% Knobs
vars.inpath   = 'invivo_full.mat';
vars.outpath  = 'recon_invivo_full_par7.mat';
vars.nLines   = 25; % [5, 10, 15, 20, 25, 30];
vars.seed     = 5;
vars.r        = nan;
vars.lambdaL  = 1e100;
vars.lambdaS  = 0.001;
vars.lambdaB  = 0.01;
vars.lambdaB0 = 0.01;
vars.dr       = 1;
vars.nIters   = 50;

% DRPCA opts
%opts.A       = A;
%opts.sdim    = [ny, nx, nt];
opts.pdim     = [8, 8, 5];
opts.pgap     = [2, 2, 2];
opts.type     = 'hard';
%opts.dr      = dr;
opts.ddim     = [prod(opts.pdim(1:(end - 1))), opts.pdim(end)];
%opts.nIters  = vars.nIters;
opts.nItersDB = 1;
opts.nItersLS = 5;
%opts.L0      = reshape(init.Lhat,[],nt);
%opts.S0      = reshape(init.Shat,[],nt);
opts.D0       = dctmtx(prod(opts.pdim));
%opts.Xtrue   = reshape(data.Xtrue,[],nt);
opts.accel    = false;
opts.tau      = 2; % Since lambdaL = inf
opts.flag     = 1;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% ktSLR initialization
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Knobs
vars0.p        = 1;
vars0.mu1      = 0;
vars0.mu2      = 0.001;
vars0.nItersO  = 5;
vars0.nItersI  = 10;
vars0.nItersCG = 5;

% ktSLR opts
%opts0.nItersO  = vars0.nItersO;
%opts0.nItersI  = vars0.nItersI;
%opts0.nItersCG = vars0.nItersCG;
opts0.beta1     = 1;
opts0.beta2     = 1;
opts0.beta1rate = 50;
opts0.beta2rate = 25;
opts0.stepSize  = [1, 1, 0.303];
