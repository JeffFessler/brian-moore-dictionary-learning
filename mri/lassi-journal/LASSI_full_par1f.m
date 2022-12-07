function [vars, opts, vars0, opts0] = LASSI_full_par1f()
% Syntax:   [vars, opts, vars0, opts0] = LASSI_full_par1f();

% Knobs
vars.inpath    = 'otazo_full.mat';
vars.rawpath   = 'LASSI_full_par1f/data.mat';
vars.outpath   = 'LASSI_full_par1f.mat';
vars.ps        = 1 ./ [4, 8, 12, 16, 20, 24];
vars.seed      = 1:5;
vars.r         = nan;
vars.lambdaL   = 0.50;
vars.lambdaS   = 0.02;
vars.lambdaB   = 0.03;
vars.lambdaB0  = 0.03;
vars.dr        = 1;
vars.nIters    = 50;

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
opts.tau      = 1;
opts.flag     = 1;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% k-t SLR initialization
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Knobs
vars0.p        = 1;
vars0.mu1      = 0.2;
vars0.mu2      = 0.003;
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
