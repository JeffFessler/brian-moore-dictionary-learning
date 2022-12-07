function [vars, opts] = LASSI_R8_par1e()
% Syntax:   [vars, opts] = LASSI_R8_par1e();

% Knobs
vars.inpath   = 'otazo_R8.mat';
vars.initpath = 'otazo_R8_ktslr_nrmse_0X.mat'; % L = 0, S = X
vars.rawpath  = 'LASSI_R8_par1e/data.mat';
vars.outpath  = 'LASSI_R8_par1e.mat';
vars.r        = nan;
vars.lambdaL  = [0.1, 0.3, 0.5];
vars.lambdaS  = [0.01, 0.02, 0.03];
vars.lambdaB  = [0.02, 0.03, 0.04];
vars.lambdaB0 = [0.02, 0.03, 0.04];
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
