function [vars, opts] = road_rpca_par1b()
% Syntax:   [vars, opts] = road_rpca_par1b();

% Knobs
vars.inpath   = 'data/road396.mat';
vars.rawpath  = 'road_rpca_par1b/data.mat';
vars.outpath  = 'road_rpca_par1b.mat';
vars.idim     = [128, nan, 200];
vars.nIters   = 400;
vars.algo     = 'paint2';
vars.method   = 1;
vars.SNR      = inf;
% sweep
vars.p        = [0.7];
vars.seed     = [1];
vars.lambda   = [0, logspace(-3,0,6)];
vars.mu       = logspace(-1,1,5) / sqrt(128 * 171);

% robustPCA opts
opts.A        = 1;
%opts.M;
opts.T        = 1;
opts.r        = nan;
%opts.nIters
%opts.L0;
%opts.S0;
%opts.Xtrue
opts.accel    = true;
opts.tau      = 0.99 / (1 + opts.accel);
opts.flag     = 1;
