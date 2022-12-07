function [vars, opts] = LASSI_R8_Dover_par2()
% Syntax:   [vars, opts] = LASSI_R8_Dover_par2();

% Knobs
vars.inpath   = 'otazo_R8.mat';
vars.initpath = 'otazo_R8_lps_visual.mat';
vars.rawpath  = 'LASSI_R8_Dover_par2/data.mat';
vars.outpath  = 'LASSI_R8_Dover_par2.mat';
vars.r        = nan;
vars.lambdaL  = 0.5;
vars.lambdaS  = 0.02;
vars.lambdaB  = 0.03;
vars.lambdaB0 = 0.03;
vars.dr       = 1;
vars.nIters   = 50;
vars.np       = [-240, -160, -80, 0];

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
