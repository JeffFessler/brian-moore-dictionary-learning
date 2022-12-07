function [vars, opts] = otazo_R8_lassi_Sonly_par5()
% Syntax:   [vars, opts] = otazo_R8_lassi_Sonly_par5();

% Knobs
vars.inpath   = 'otazo_R8.mat';
vars.initpath = 'otazo_R8_lps_visual.mat';
vars.rawpath  = 'otazo_R8_lassi_Sonly_par5/data.mat';
vars.outpath  = 'otazo_R8_lassi_Sonly_par5.mat';
vars.r        = nan;
vars.lambdaL  = 1e100;
vars.lambdaS  = 0.474;
vars.lambdaB  = 0.008;
vars.lambdaB0 = 0.008;
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
opts.fixedD   = true; % Fixed dictionary
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
