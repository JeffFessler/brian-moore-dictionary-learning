function [vars, opts] = pincat_LASSI_R8_Sonly_par2n()
% Syntax:   [vars, opts] = pincat_LASSI_R8_Sonly_par2n();

% Knobs
vars.inpath   = 'pincat_R8n.mat';
vars.initpath = 'pincat_R8n_ktslr_nrmse.mat';
vars.rawpath  = 'pincat_LASSI_R8_Sonly_par2n/data.mat';
vars.outpath  = 'pincat_LASSI_R8_Sonly_par2n.mat';
vars.r        = nan;
vars.lambdaL  = 1e100;
vars.lambdaS  = 0.158;
vars.lambdaB  = [0.04, 0.06, 0.08];
vars.lambdaB0 = [0.04, 0.06, 0.08];
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
opts.tau      = 2; % Since lambdaL = inf
opts.flag     = 1;
