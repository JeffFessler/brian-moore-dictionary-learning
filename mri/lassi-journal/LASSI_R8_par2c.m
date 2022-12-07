function [vars, opts] = LASSI_R8_par2c()
% Syntax:   [vars, opts] = LASSI_R8_par2c();

% Knobs
vars.inpath   = 'otazo_R8.mat';
vars.initpath = 'otazo_R8_lps_visual.mat';
vars.rawpath  = 'LASSI_R8_par2c/data.mat';
vars.outpath  = 'LASSI_R8_par2c.mat';
vars.r        = nan;
vars.lambdaL  = [0.1, 0.2, 0.3, 0.4];
vars.lambdaS  = [0.02, 0.03, 0.04];
vars.lambdaB  = [0.001, 0.003, 0.005];
vars.lambdaB0 = [0.001, 0.003, 0.005];
vars.dr       = 1;
vars.nIters   = 50;

% DRPCA opts
%opts.A       = A;
%opts.sdim    = [ny, nx, nt];
opts.pdim     = [8, 8, 5];
opts.pgap     = [2, 2, 2];
opts.type     = 'soft';
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
