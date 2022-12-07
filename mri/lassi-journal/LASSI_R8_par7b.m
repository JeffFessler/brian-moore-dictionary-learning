function [vars, opts] = LASSI_R8_par7b()
% Syntax:   [vars, opts] = LASSI_R8_par7b();

% Knobs
vars.inpath   = 'otazo_R8.mat';
vars.initpath = 'otazo_R8_lps_visual.mat';
vars.rawpath  = 'LASSI_R8_par7b/data.mat';
vars.outpath  = 'LASSI_R8_par7b.mat';
vars.p        = 0.5; % "Firm" thresholding
vars.r        = nan;
vars.lambdaL  = 0.50 * logspace(-1,1,15);
vars.lambdaS  = 0.02;
vars.lambdaB  = 0.03;
vars.lambdaB0 = 0.03;
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
