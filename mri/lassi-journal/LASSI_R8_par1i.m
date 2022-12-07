function [vars, opts] = LASSI_R8_par1i()
% Syntax:   [vars, opts] = LASSI_R8_par1i();

% Knobs
vars.inpath   = 'otazo_R8.mat';
vars.initpath = 'otazo_R8_lps_mse_0X.mat';
vars.rawpath  = 'LASSI_R8_par1i/data.mat';
vars.outpath  = 'LASSI_R8_par1i.mat';
vars.r        = nan;
vars.lambdaL  = 0.5;
vars.lambdaS  = 0.02;
vars.lambdaB  = [0.002, 0.004, 0.006, 0.008, ...
                 0.010, 0.020, 0.030, 0.040, ...
                 0.060, 0.100, 0.150, 0.200];
vars.lambdaB0 = [0.002, 0.004, 0.006, 0.008, ...
                 0.010, 0.020, 0.030, 0.040, ...
                 0.060, 0.100, 0.150, 0.200];
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
opts.tau      = 1;
opts.flag     = 1;
