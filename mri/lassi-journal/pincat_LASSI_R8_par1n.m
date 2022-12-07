function [vars, opts] = pincat_LASSI_R8_par1n()
% Syntax:   [vars, opts] = pincat_LASSI_R8_par1n();

% Knobs
vars.inpath   = 'pincat_R8n.mat';
vars.initpath = 'pincat_R8n_lps_nrmse.mat';
vars.rawpath  = 'pincat_LASSI_R8_par1n/data.mat';
vars.outpath  = 'pincat_LASSI_R8_par1n.mat';
vars.r        = nan;
vars.lambdaL  = [0.05, 0.1, 0.5, 1.0];
vars.lambdaS  = [0.01, 0.05, 0.1];
vars.lambdaB  = [0.005, 0.01, 0.02];
vars.lambdaB0 = [0.005, 0.01, 0.02];
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
