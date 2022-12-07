function [vars, opts] = LASSI_R8_par1g()
% Syntax:   [vars, opts] = LASSI_R8_par1g();

% Knobs
vars.inpath   = 'otazo_R8.mat';
vars.initpath = 'otazo_R8_lps_mse_0X.mat'; % L = 0, S = X
vars.rawpath  = 'LASSI_R8_par1g/data.mat';
vars.outpath  = 'LASSI_R8_par1g.mat';
vars.r        = nan;
vars.lambdaL  = 1e100;
vars.lambdaS  = [0.010, 0.015, 0.020];
vars.lambdaB  = [0.015, 0.025, 0.035];
vars.lambdaB0 = [0.015, 0.025, 0.035];
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
opts.fixedD   = true;                           % Fixed dictionary
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
