function [vars, opts] = invivo_LASSI_R8_Sonly_par2()
% Syntax:   [vars, opts] = invivo_LASSI_R8_Sonly_par2();

% Knobs
vars.inpath   = 'invivo_R8.mat';
vars.initpath = 'invivo_R8_ktslr_nrmse.mat';
vars.rawpath  = 'invivo_LASSI_R8_Sonly_par2/data.mat';
vars.outpath  = 'invivo_LASSI_R8_Sonly_par2.mat';
vars.r        = nan;
vars.lambdaL  = 1e100;
vars.lambdaS  = 0.01 * logspace(-1.5,1.5,11);
vars.lambdaB  = [0.01, 0.02, 0.04];
vars.lambdaB0 = [0.01, 0.02, 0.04];
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
opts.tau      = 2; % Since lambdaL = inf
opts.flag     = 1;
