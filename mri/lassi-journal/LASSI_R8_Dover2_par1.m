function [vars, opts] = LASSI_R8_Dover2_par1()
% Syntax:   [vars, opts] = LASSI_R8_Dover2_par1();

% Knobs
vars.inpath   = 'otazo_R8.mat';
vars.initpath = 'otazo_R8_lps_mse_0X.mat'; % L = 0, S = X
vars.rawpath  = 'LASSI_R8_Dover2_par1/data.mat';
vars.outpath  = 'LASSI_R8_Dover2_par1.mat';
vars.r        = nan;
vars.lambdaL  = 0.5;
vars.lambdaS  = 0.02;
vars.lambdaB  = [0.005, 0.01, 0.02, 0.03, 0.05, 0.07, 0.09];
vars.lambdaB0 = [0.005, 0.01, 0.02, 0.03, 0.05, 0.07, 0.09];
vars.dr       = 1;
vars.np       = [-200, -100, 0, 100, 200];  % # dictionary atoms
vars.nIters   = 50;

% DRPCA opts
%opts.A       = A;
%opts.sdim    = [ny, nx, nt];
opts.pdim     = [8, 8, 5];
opts.pgap     = [2, 2, 2];
%opts.np      = np;
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
