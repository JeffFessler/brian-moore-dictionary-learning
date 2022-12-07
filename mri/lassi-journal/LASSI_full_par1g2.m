function [vars, opts, vars0, opts0] = LASSI_full_par1g2()
% Syntax:   [vars, opts, vars0, opts0] = LASSI_full_par1g2();

% Knobs
vars.inpath    = 'otazo_full.mat';
vars.rawpath   = 'LASSI_full_par1g2/data.mat';
vars.outpath   = 'LASSI_full_par1g2.mat';
vars.ps        = 1 ./ [4, 8, 12, 16, 20, 24];
vars.seed      = 1:5;
vars.r         = nan;
vars.lambdaL   = 0.50;
vars.lambdaS   = 0.02;
vars.lambdaB   = 0.06;
vars.lambdaB0  = 0.06;
vars.dr        = 1;
vars.nIters    = 50;

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

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% FFT initialization
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

vars0 = [];
opts0 = [];
