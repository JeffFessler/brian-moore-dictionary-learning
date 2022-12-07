function [vars, opts] = invivo_LASSI_R8_Diter_par1()
% Syntax:   [vars, opts] = invivo_LASSI_R8_Diter_par1();

% Knobs
vars.inpath   = 'invivo_R8.mat';
vars.initpath = 'invivo_R8_lps_nrmse.mat';
vars.rawpath  = 'invivo_LASSI_R8_Diter_par1/data.mat';
vars.outpath  = 'invivo_LASSI_R8_Diter_par1.mat';
vars.r        = nan;
vars.lambdaL  = 0.05;
vars.lambdaS  = 0.01;
vars.lambdaB  = 0.01;
vars.lambdaB0 = 0.01;
vars.dr       = 5;
vars.nIters   = [50, 25, 10,  5];
vars.nItersDB = [ 1,  2,  5, 10];
vars.nItersLS = [ 5, 10, 25, 50];

% DRPCA opts
%opts.A        = A;
%opts.sdim     = [ny, nx, nt];
opts.pdim      = [8, 8, 5];
opts.pgap      = [2, 2, 2];
opts.type      = 'hard';
%opts.dr       = dr;
opts.ddim      = [prod(opts.pdim(1:(end - 1))), opts.pdim(end)];
%opts.nIters   = nIters;
%opts.nItersDB = nItersDB;
%opts.nItersLS = nItersLS;
%opts.L0       = reshape(init.Lhat,[],nt);
%opts.S0       = reshape(init.Shat,[],nt);
opts.D0        = dctmtx(prod(opts.pdim));
%opts.Xtrue    = reshape(data.Xtrue,[],nt);
opts.accel     = false;
opts.tau       = 1;
opts.flag      = 1;
