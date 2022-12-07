function [vars, opts] = invivo_LASSI_R8_par1b()
% Syntax:   [vars, opts] = invivo_LASSI_R8_par1b();

% Knobs
vars.inpath   = 'invivo_R8.mat';
vars.initpath = 'invivo_R8_lps_nrmse.mat';
vars.rawpath  = 'invivo_LASSI_R8_par1b/data.mat';
vars.outpath  = 'invivo_LASSI_R8_par1b.mat';
vars.r        = nan;
vars.lambdaL  = [0.005, 0.01, 0.025, 0.05];
vars.lambdaS  = [0.005, 0.01, 0.02];
vars.lambdaB  = [0.0025, 0.005, 0.01];
vars.lambdaB0 = [0.0025, 0.005, 0.01];
vars.dr       = [1, 5];
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
opts.tau      = 1;
opts.flag     = 1;
