function [vars, opts] = invivo_LASSI_R8_par1f()
% Syntax:   [vars, opts] = invivo_LASSI_R8_par1f();

% Knobs
vars.inpath   = 'invivo_R8.mat';
vars.initpath = 'invivo_R8_ktslr_nrmse_pca2.mat'; % L = PCA(X,2), S = X - L
vars.rawpath  = 'invivo_LASSI_R8_par1f/data.mat';
vars.outpath  = 'invivo_LASSI_R8_par1f.mat';
vars.r        = nan;
vars.lambdaL  = [1.0, 2.0, 4.0, 6.0];
vars.lambdaS  = [0.05, 0.10];
vars.lambdaB  = [0.04, 0.08];
vars.lambdaB0 = [0.04, 0.08];
vars.dr       = 5;
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
