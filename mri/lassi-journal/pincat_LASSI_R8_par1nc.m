function [vars, opts] = pincat_LASSI_R8_par1nc()
% Syntax:   [vars, opts] = pincat_LASSI_R8_par1nc();

% Knobs
vars.inpath   = 'pincat_R8n.mat';
vars.initpath = 'pincat_R8n_lps_nrmse.mat';
vars.rawpath  = 'pincat_LASSI_R8_par1nc/data.mat';
vars.outpath  = 'pincat_LASSI_R8_par1nc.mat';
vars.r        = nan;
vars.lambdaL  = [2, 4, 6, 8, 10];
vars.lambdaS  = 0.1;
vars.lambdaB  = [0.02, 0.04, 0.06];
vars.lambdaB0 = [0.02, 0.04, 0.06];
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
opts.tau      = 1;
opts.flag     = 1;
