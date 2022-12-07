function [vars, opts] = pincat_LASSI_R8_par1nh2()
% Syntax:   [vars, opts] = pincat_LASSI_R8_par1nh2();

% Knobs
vars.inpath   = 'pincat_R8n.mat';
vars.initpath = 'pincat_R8n_ktslr_nrmse_0X.mat'; % L = 0, S = X
vars.rawpath  = 'pincat_LASSI_R8_par1nh2/data.mat';
vars.outpath  = 'pincat_LASSI_R8_par1nh2.mat';
vars.r        = nan;
vars.lambdaL  = [2.0, 3.0];
vars.lambdaS  = [0.1, 0.15, 0.2];
vars.lambdaB  = 0.06;
vars.lambdaB0 = 0.06;
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
