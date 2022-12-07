function [vars, opts] = LASSI_R8_par3c()
% Syntax:   [vars, opts] = LASSI_R8_par3c();

% Knobs
vars.inpath   = 'otazo_R8.mat';
vars.initpath = 'otazo_R8_lps_visual.mat';
vars.rawpath  = 'LASSI_R8_par3c/data.mat';
vars.outpath  = 'LASSI_R8_par3c.mat';
vars.r        = 1;
vars.lambdaL  = nan;
vars.lambdaS  = 0.02;
vars.lambdaB  = [0.020, 0.021, 0.022, 0.023, 0.024, 0.025];
vars.lambdaB0 = [0.020, 0.021, 0.022, 0.023, 0.024, 0.025];
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
